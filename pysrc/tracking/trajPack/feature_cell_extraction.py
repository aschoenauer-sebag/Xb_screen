import getpass, os,pdb, operator, sys
import numpy as np
import cPickle as pickle
from optparse import OptionParser
from operator import itemgetter
from collections import Counter, defaultdict
from scipy.stats import scoreatpercentile

from tracking.trajPack import featuresSaved, featuresHisto, featuresNumeriques
from util.listFileManagement import fromShareToCBIO, appendingControl
from tracking.trajPack.clustering import histConcatenation
from tracking.PyPack.fHacktrack2 import initXml, finirXml
from tracking.histograms.k_means_transportation import DIVERGENCES, _distances
from tracking.histograms.transportation import costMatrix, computingBins
from util.listFileManagement import expSi, strToTuple
from util.sandbox import concatCtrl
from tracking.plots import plotAlignedTraj

from util import settings
from scipy.stats.stats import ks_2samp


nb_exp_list=[100, 500, 1000, 2000]

def compareDistributions(bins_type, nb_exp, folder="../resultData/features_on_films"):
    '''
    Ici on va comparer la distribution d'une feature dans differents tirages de nb_exp films.
    Le but est de trouver in fine le nb_exp minimal pour lequel la distribution est stable.
    Ca va me dire a quel point je peux paralleliser le calcul de la distance d'un film a ses controles
    '''
    
    pvalues = defaultdict(list)
    for feature in featuresNumeriques:
        print feature
        l = filter(lambda x: feature in x and bins_type in x and int(x.split('_')[2])>0.75*nb_exp and int(x.split('_')[2])<1.25*nb_exp, os.listdir(folder))
        for file_ in l:
            iter_ = int(file_.split('_')[-1].split('.')[0])
            print '----', iter_
            for other_file in filter(lambda x: int(x.split('_')[-1].split('.')[0])>iter_, l):
                iter2 = int(other_file.split('_')[-1].split('.')[0])
                print iter2
                f=open(os.path.join(folder, file_))
                l1=pickle.load(f); f.close()
    
                f=open(os.path.join(folder, other_file))
                l2=pickle.load(f); f.close()
                if bins_type =='quantile':
                    ks, pval = ks_2samp(l1[1], l2[1])
                else:
                    ks, pval = ks_2samp(l1, l2)
                    
                pvalues[feature].append(pval)
    for feature in featuresNumeriques:
        print feature, np.mean(pvalues[feature]), scoreatpercentile(pvalues[feature], 90)
        
    f=open(os.path.join(folder, 'pvalues_{}_{}'.format(bins_type, nb_exp)), 'w')
    pickle.dump(pvalues, f); f.close()
    return

def computeBins(bins_type,bin_size, folder="../resultData/features_on_films"):
    if type(bin_size)!=list:
        bin_size=[bin_size for k in range(len(featuresNumeriques))]
        
    if bins_type=="minmax":
        minMax=np.empty(shape=(len(featuresNumeriques), 2))
    else:
        quantiles = np.empty(shape=(len(featuresNumeriques),), dtype=object)
    
    for i,feature in enumerate(featuresNumeriques):
        print feature
        file_ = "quantiles_{}_11989_0.pkl".format(feature)
        f=open(os.path.join(folder, file_), 'r')
        distribution = pickle.load(f)[1]
        f.close()
        
        if bins_type=="minmax":
            minMax[i,0]=np.min(distribution)
            minMax[i,1]=np.max(distribution)
        elif bins_type=="quantile":
            quantiles[i] = [scoreatpercentile(distribution, per) for per in [k*100/float(bin_size[i]) for k in range(bin_size[i]+1)]]
            
    f=open(os.path.join(folder, 'distExp_123etctrl_{}_{}.pkl'.format(bins_type, bin_size[0])), 'w')
    if bins_type=="minmax":
        pickle.dump(minMax, f)
    else:
        pickle.dump(quantiles, f)
    f.close()
    return


def _writeXml(plate, well, resultCour, num=1):
    tmp=''
    for traj in resultCour:
        tmp+="\n      <Marker_Type>\n        <Type>{0}</Type>".format(num)#for the color of the cells in the xml file
        
        for t, x, y in zip(traj[0], traj[1], traj[2]):
            numFramePourXml = t+1
            
            tmp+="\n        <Marker>\n          <MarkerX>{1}</MarkerX>\n          <MarkerY>{2}</MarkerY>\n          <MarkerZ>{0}</MarkerZ>\n        </Marker>".format(numFramePourXml, int(x), int(y))
            
        tmp+="\n      </Marker_Type>"
    print len(tmp)
    assert len(tmp)>0
    return tmp

def compFeatures(data, length, who, ctrlStatus, percentile, *args):
    how_many=5
    outputFolder='../prediction/annotations'
    plwList=defaultdict(int)
    whoList={}
    extractors=[]
    for feature_name in args:
        extractor = cellExtractor(feature_name)
        extractors.append(extractor)
        bestMovies = extractor.dataToBestMovies(data, length, who, ctrlStatus, percentile, how_many=how_many)
        #plwList.extend(bestMovies.keys())
        whoList[feature_name]=bestMovies
        for el in bestMovies:
            plwList[el]+=len(bestMovies[el])
    plwList={el:plwList[el] for el in filter(lambda x: np.all([x in whoList[feature_name] for feature_name in args]), plwList.keys())}
    sort= sorted(plwList, key=operator.itemgetter(1))
    chosen=sort[-how_many:]
    num=1
    for el in chosen:
        plate, well=el
        filename = "PL{}___P{}___T00000.xml".format(plate, well)
        file_ = open(os.path.join(outputFolder, filename), 'w')
        file_.write(initXml())
        for i,feature_name in enumerate(args):
            resultCour = whoList[feature_name][el]
            print len(resultCour)#donc attention ici on n'a que les index des trajectoires dans le fichier coord
            file_.write(extractors[i].wellToXml(length,who, plate, well, resultCour, outputFolder, num)) # _writeXml(plate, well, resultCour, num))
            num+=1
        file_.write(finirXml())
        file_.close()
    return chosen
        
    
class cellExtractor():
    def __init__(self, siRNA, settings_file,div_name='total_variation', bin_type = 'quantile',
                 bin_size=50,lambda_=10,M=None, verbose=0):

        self.siRNA = siRNA
        self.settings = settings.Settings(settings_file, globals())
        
        self.div_name =div_name
        self.bin_type = bin_type
        self.bin_size = bin_size        
        self.verbose=verbose
        self.currInterestFeatures = featuresNumeriques
        if self.settings.histDataAsWell:
            raise AttributeError
    
    def _findExperiment(self):
        '''
        Purpose : find experiments that are related to this siRNA
        '''
        yqualDict=expSi(self.settings.quality_control_file, sens=0)
        l= strToTuple(yqualDict[self.siRNA], os.listdir(self.settings.data_folder))
        if self.verbose>0:
            print l
        return l
    
    def getData(self, histDataAsWell):
        '''
        Ici sur toutes les experiences liees au siRNA self.siRNA on construit l'histogramme des features numeriques
        R contains all the information for all experiments, and length is the list of the length of all experiments,
        meaning the number of trajectories in each
        '''
        histDict = defaultdict(list)

        r, _, _,_, length, _, _, _ = histConcatenation(self.settings.data_folder, self.expList, self.settings.mitocheck_file,
                                        self.settings.quality_control_file, verbose=self.verbose)
                    
        for i in range(len(length)):
            for k,feature in enumerate(self.currInterestFeatures):
                histDict[feature].append(r[np.sum(length[:i]):np.sum(length[:i+1]),k])
        
        f=open(os.path.join(self.settings.result_folder, 'distExp_123etctrl_{}_{}.pkl'.format(self.bin_type, self.bin_size)))
        bins = pickle.load(f); f.close()
        
        histogrammes, bins = computingBins(histDict, [self.bin_size for k in range(len(self.currInterestFeatures))], self.bin_type, previous_binning=bins)
                    
        return length, histogrammes, bins
    
    def ctrlHistograms(self, bins):
        '''
        Since we want to calculate the length of an experiment to its controls, we only need to keep one table per
        plate, and not one table per control experiment
        '''
        plates = Counter(np.array(self.expList)[:,0]).keys()

        histDict = defaultdict(list)
        length=[]
        for plate in plates:
            ctrlExpList = appendingControl([plate])
            try:
                curr_r, _, _,_, curr_length, _, _, _ = histConcatenation(self.settings.data_folder, ctrlExpList, self.settings.mitocheck_file,
                                            self.settings.quality_control_file, verbose=self.verbose)
            except:
                print "Problem with controls from plate {}".format(plate)
            else:    
                length.append(np.sum(curr_length))
                assert(np.sum(curr_length)==curr_r.shape[0])                        
                for k, feature in enumerate(self.currInterestFeatures):
                    histDict[feature].append(curr_r[:,k])
        assert(len(length)==len(plates))
        assert(len(histDict["entropy1"])==len(plates))
        pdb.set_trace()
                    
        histogrammes, bins = computingBins(histDict, [self.bin_size for k in range(len(self.currInterestFeatures))], self.bin_type, previous_binning=bins)
        
        return plates, length, histogrammes
    
    def parameters(self):
        r={
           'bins_type': self.bins_type,
           'bin_size':self.bin_size,
           'div_name':self.div_name,
           'lambda':self.lambda_,
           'cost_type':self.cost_type
           }
        return tuple(sorted(zip(r.keys(), r.values()), key=itemgetter(0)))
    
    def calculateDistances(self, plates, histogrammes, ctrl_histogrammes, length, ctrl_length):
        distances = np.zeros(shape=(len(self.expList),len(self.currInterestFeatures)))
        for i,experiment in enumerate(self.expList):
            corresponding_ctrl = plates.index(experiment[0])
            ctrl_hist = ctrl_histogrammes[corresponding_ctrl]

            for k, feature in enumerate(self.currInterestFeatures):
                dist_weights = np.zeros(shape=(len(self.currInterestFeatures)+1,))
                dist_weights[k+1]=1
                distances[i, k]=_distances(histogrammes[np.sum(length[:i]):np.sum(length[:i+1])],
                                           ctrl_hist[np.newaxis,:],div_name=self.div_name, lambda_=self.lambda_, M=None,
                                     mat_hist_sizes=[[self.bin_size for j in range(len(self.currInterestFeatures))]], nb_feat_num=0, dist_weights=dist_weights)
                
        return distances
    
    def saveResults(self, distances):
        if self.settings.outputFile.format(self.siRNA) in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, self.settings.outputFile.format(self.siRNA)))
            d=pickle.load(f)
            f.close()
            try:
                l=d[self.parameters()]
                sys.stderr.write("Bizarre, on dirait que cela a deja ete calcule")
            except KeyError:
                l=None
        else:
            d={}
            l=None
        l=distances if l is None else np.vstack((l, distances))
        d.update({self.parameters():l})
        f=open(os.path.join(self.settings.result_folder, self.settings.outputFile), 'w')
        pickle.dump(d, f)
        f.close()
        return
    
    def __call__(self):
        #i. getting experiments corresponding to siRNA
        self.expList=self._findExperiment()
        if self.expList==[]:
            sys.stderr.write("No experiments for siRNA {}".format(self.siRNA))
            return
        
        #iii.calculate the histograms and binnings of experiments, using pre-computed binnings, ON THE EXPERIMENTS ONLY
        length, histogrammes, bins = self.getData(self.settings.histDataAsWell)
        
        #iv. calculate the histograms of corresponding controls. Should be per plate then, a dictionary
        plates, ctrl_length, ctrl_histogrammes = self.ctrlHistograms(bins)
        
        #iii. calculate the distance from each experiment to its control
        distances = self.calculateDistances(plates, histogrammes, ctrl_histogrammes, length, ctrl_length)
        
        #iv. save results
        self.saveResults(distances)
        
        return
        
#    def trajToWells(self, length, who, trajectories):
#        cum=np.cumsum(length)
#        cum=np.insert(cum, 0, 0)
#        #l=[0]; l.extend(cum); cum=np.array(l)
#        result=defaultdict(list)
#        for traj_index in trajectories:
#            indexCour = np.where(cum<=traj_index)[0][-1]
#            result[who[indexCour]].append(traj_index-cum[indexCour])
#        return result
#        
#    def wellToPlot(self, length,who, all_, result_folder, median=False):
#        assert getpass.getuser()!='lalil0u', 'You are not running this on the right computer'
#        index=featuresSaved.index(self.feature_name)
#        folder='/share/data20T/mitocheck/tracking_results'
#        basename = 'traj'
#        resultCour=[]; sizes=[0]
#        valCour=[]
#        #if median: basename+='_median'
#        for el in all_:
#            trajectories=all_[el]
#            print el
#            for plate, well in trajectories:
#                print plate, well
#                if self.verbose:
#                    print "Taking care of plate {}, well {}".format(plate, well)
#        
#                f=open(os.path.join(folder,plate, 'hist_tabFeatures_{}.pkl'.format(well)))
#                tab, coord, _=pickle.load(f); f.close()
#                resultCour.extend(np.array(coord)[trajectories[(plate, well)]])
#                valCour.extend(tab[trajectories[(plate, well)], index])
#            sizes.append(len(resultCour))
#
#        XX=[]; YY=[]
#        for k in range(len(resultCour)):
#            X=np.array(resultCour[k][1]); X-=X[0]; XX.append(X)
#            Y=np.array( resultCour[k][2]); Y-=Y[0]; YY.append(Y)
#        minx,maxx = min([min(X) for X in XX]), max([max(X) for X in XX])
#        miny,maxy = min([min(Y) for Y in YY]), max([max(Y) for Y in YY])
#        
#        #saving trajectories 
#        f=open('trajectories_{}.pkl'.format(self.feature_name), 'w')
#        pickle.dump([all_, resultCour, valCour, sizes],f)
#        f.close()
#        
#        i=0
#        for el in all_:
#            for k in range(sizes[i], min(sizes[i]+50,sizes[i+1])):
#                plotAlignedTraj(XX[k], YY[k], minx, maxx, miny, maxy, show=False, 
#                                name=os.path.join(result_folder, '{}_{}_{}_{}.png'.format(basename,el, self.feature_name, k)), val=valCour[k])
#            i+=1
#        return 1#_writeXml(plate, well, resultCour)
#    
#    def wellToXml(self, length,who, plate, well, trajectories, result_folder, num):
#
#        assert getpass.getuser()!='lalil0u', 'You are not running this on the right computer'
#        folder='/share/data/mitocheck/tracking_results'
#
#        if self.verbose:
#            print "Taking care of plate {}, well {}".format(plate, well)
#
#        f=open(os.path.join(folder,plate, 'hist_tabFeatures_{}.pkl'.format(well)))
#        _, coord, _, _, _, _=pickle.load(f); f.close()
#        resultCour=np.array(coord)[trajectories]
#        return _writeXml(plate, well, resultCour, num)
#        
#    def find(self, data,ctrlStatus,length, percentile, how_many=5, div_name='total_variation', nb_feat_num=12, sample_size=20,
#                lambda_=10,M=None):
#        if self.numerical:
#            return self._findNumerical(data, percentile, how_many)
#        return self._findHistogram(data,ctrlStatus,length, percentile, how_many, div_name,nb_feat_num, sample_size,lambda_, M)
#    
#    def _findNumerical(self, data, percentile, how_many):
#        index=featuresSaved.index(self.feature_name)
#        scoreD = scoreatpercentile(data[:,index], percentile)
#        down=np.where(data[:,index]<=scoreD)
#        scoreU = scoreatpercentile(data[:,index], 100-percentile)
#        up=np.where(data[:,index]>=scoreU)
#            
#        med0, med1=scoreatpercentile(data[:,index], 48), scoreatpercentile(data[:,index], 52)
#        ctrl0 = np.where(data[:,index]>=med0)
#        ctrl1=np.where(data[:, index]<=med1)
#        ctrl=filter(lambda x: x in ctrl1[0], ctrl0[0])
#        
#        if self.verbose:
#            print "Feature {}".format(self.feature_name)
#            print "{} retrieved, percentile {}, value {}".format(up[0], 100-percentile, scoreU)
#            print "{} retrieved, percentile {}, value {}".format(down[0],percentile, scoreD)
#            print "{} retrieved, percentile {}, value {}".format(np.array(ctrl), percentile, med0)
#        return down[0], up[0], ctrl
#    
#    def _findHistogram(self, data,ctrlStatus,length, percentile, how_many, div_name,nb_feat_num, sample_size, lambda_, M):
#        #not implemented yet
#        assert(div_name in DIVERGENCES)
#        if div_name=='transportation' and M is None:
#            M={}
#            for m in range(self.mat_hist_sizes.shape[1]):
#                for n in range(self.mat_hist_sizes.shape[0]):
#                    M[(n,m)]=costMatrix(self.mat_hist_sizes[n,m], 1)
#        
#        print "CURRENT WEIGHTS", self.dist_weights
#        
#        random_state= np.random.mtrand._rand
#        ctrlData=concatCtrl(data, ctrlStatus, length, len(ctrlStatus)-np.sum(ctrlStatus))
#        seeds = random_state.permutation(ctrlData.shape[0])[:sample_size]
#        sample=ctrlData[seeds]
#        
#        #taille sample_size, data.shape[0]
#        distances = _distances(data, sample,div_name, lambda_, M, self.mat_hist_sizes, nb_feat_num, self.dist_weights)
#        distances=np.argmin(distances, 0)#donc doit etre de taille data.shape[0]
#        
#        #then we have a numerical feature so we can get the percentile the same way
#        score = scoreatpercentile(distances, percentile)
#        if percentile>50:
#            r=np.where(distances>=score)[0]
#        else:
#            r=np.where(distances<=score)[0]
#        if self.verbose:
#            print "{} examples retrieved for feature {}, percentile {}".format(len(r), self.feature_name, percentile)
#        return r
    
if __name__ == '__main__':
    '''
    Idea here is to calculate the distance of an experiment histogram to its control, for each numerical feature
    nb_exp should be the minimum number of experiments to get a stable distribution of the feature among
    considered experiments, so that the calculation is parallelizable
    '''
    parser = OptionParser(usage="usage: %prog [options]")    
    
    parser.add_option('--siRNA', type=str, dest='siRNA', default=None)

    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
    parser.add_option('--bin_size', type=int, dest="bin_size", default=50)
    parser.add_option('--nb_exp', type=int, dest='size', default=1000)
    
    parser.add_option("-l",type=int, dest="lambda_", default=10)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    (options, args) = parser.parse_args()
    
    settings_file = 'tracking/settings/settings_feature_extraction.py'

    extractor=cellExtractor(options.siRNA, settings_file,options.div_name, options.bins_type,
                 bin_size=options.bin_size,lambda_=options.lambda_, verbose=options.verbose)
    extractor()
