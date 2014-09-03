import getpass, os,pdb, operator, sys
import numpy as np
import cPickle as pickle

if getpass.getuser()=='aschoenauer':
    import matplotlib as mpl
    mpl.use('Agg')
    show = False
elif getpass.getuser()=='lalil0u':
    import matplotlib as mpl
    show=True

from optparse import OptionParser
from operator import itemgetter
from collections import Counter, defaultdict
from scipy.stats import pearsonr
import matplotlib.pyplot as p
import brewer2mpl
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
from scipy.stats.stats import ks_2samp, spearmanr, scoreatpercentile


nb_exp_list=[100, 500, 1000, 2000]

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

def Spearman(folder = '../resultData/features_on_films',feature=None, outputFile='pval_Spearman', sh=None, extreme=False, 
             percentile=95, saveExtreme=False, savePlots=True):
    '''
    Objective is to calculate Spearman's rank correlation coefficient between the distances given under
    different distance parameters
    '''
    if feature == None:
        featureL = featuresNumeriques
    elif type(feature)!=list:
        featureL=[feature]
    else:
        featureL = feature
    if sh==None:
        sh=show
           
    parameters = None
    fileL = filter(lambda x: 'distances_' in x and '.pkl' in x, os.listdir(folder))
    n_siRNA = len(fileL)
    set_number = 8
    result=[np.zeros(shape=(set_number,), dtype=object) for k in range(len(featureL))]
    
    for i in range(len(featureL)):
        for k in range(set_number):
            result[i][k]=[]

    siRNAL=[]
    count=0
    
    for j, file_ in enumerate(fileL):
        siRNA = file_.split('_')[-1][:-4]
        
            
        f=open(os.path.join(folder, file_))
        distances = pickle.load(f); f.close()
        
        if parameters == None:
            parameters = distances.keys()
        elif len(distances.keys())!=set_number:
            print file_
            pdb.set_trace()
        elif not np.all([parameter in distances.keys() for parameter in parameters]):
            print file_
            pdb.set_trace()
        elif np.any([len(distances[parameters[0]][:,featuresNumeriques.index(feature)])!=len(distances[parameter][:,featuresNumeriques.index(feature)]) for parameter in parameters]):
            print 'z', file_
            continue
        else:
            for param in filter(lambda x: x not in parameters, distances.keys()):
                parameters.append(param)
                
        siRNAL.extend([siRNA for k in range(len(distances[distances.keys()[0]]))])
        parameters = sorted(parameters, key=itemgetter(3, 1, 0))
        for i,feature in enumerate(featureL):
            for k,param in enumerate(parameters):
    
                l=list(distances[param][:,featuresNumeriques.index(feature)])
                count+=len(l)
                
                result[i][k].extend(l)
    pvalues_mat = np.zeros(shape=(len(featureL),set_number,set_number))
    if saveExtreme:
        targetExp={}; yeSiExp=expSi('../data/mapping/qc_export.txt', 0)
    for i, feature in enumerate(featureL):
        for k, param in enumerate(parameters):
            if i==0:
                print k, param
            pvalues_mat[i,k,k]=1
            result[i][k]=np.array(result[i][k])
            ext1 = np.where(result[i][k]>scoreatpercentile(result[i][k], percentile))[0]
            
            if saveExtreme and k==2:
                subExtreme1 = np.where(result[i][k]<scoreatpercentile(result[i][k], 100))[0]
                commonIndices = filter(lambda x: x in ext1, subExtreme1)
                orderedIndices = np.array(sorted(zip(commonIndices, result[i][k][commonIndices]), key=itemgetter(1), reverse=True))[:,0]

                targetExp[feature]=[]
                #np.random.shuffle(commonIndices)
                for el in np.array(orderedIndices, dtype=int):
                    L = np.where(np.array(siRNAL)==siRNAL[el])
                    currExp=yeSiExp[siRNAL[el]]
                    try:
                        assert len(currExp)==len(L[0])
                    except AssertionError:
                        pdb.set_trace()
                    targetExp[feature].append(currExp[np.where(L[0]==el)[0]])
                    
            for j in range(k+1, len(parameters)):
                if extreme==False:
                    pvalues_mat[i, k,j]=spearmanr(result[i][k], result[i][j])[0]
                else:
                    result[i][j]=np.array(result[i][j])
                    ext2=np.where(result[i][j]>scoreatpercentile(result[i][j], percentile))[0]
                    
        #here there's no reason that this calculus makes sense because no reason that the extremes for different parameter sets
        #are the same.
                    #pvalues_mat[i, k,j]=spearmanr(extreme1[-min(len(extreme1), len(extreme2)):], extreme2[-min(len(extreme1), len(extreme2)):])[0]
                    pvalues_mat[i, k,j]=(spearmanr(result[i][k][ext1], result[i][j][ext1])[0]+spearmanr(result[i][k][ext2], result[i][j][ext2])[0])/2
                pvalues_mat[i,j,k]=pvalues_mat[i, k,j]
    if saveExtreme:
        f=open(os.path.join(folder, 'targetExperiments_{}.pkl'.format(percentile)), 'w')
        pickle.dump(targetExp, f); f.close()
        
    if savePlots or sh:
        cmap = brewer2mpl.get_map('Blues', 'sequential', 3).mpl_colormap
        for i,feature in enumerate(featureL):
            f, ax = p.subplots(1, sharex=True, sharey=True)
            zou=ax.pcolormesh(pvalues_mat[i],
                          cmap=cmap,
                          norm=mpl.colors.Normalize(vmin=0, vmax=1),
                          edgecolors = 'None')
            ax.set_title("P-values, Spearman's rank correlation for feature {}, extreme experiments? {}".format(feature, extreme))
            f.colorbar(zou)
            if sh:
                p.show()
            else:
                p.savefig(os.path.join(folder, outputFile+'{}_e{}.png'.format(feature, extreme)))
    return

def plotDistanceDist(folder = '../resultData/features_on_films',siRNA_folder = '../data', feature=None, sh=None):
    if feature == None:
        featureL = featuresNumeriques
    else:
        featureL=[feature]
    f=open(os.path.join(siRNA_folder, 'siRNA_MAR_genes.pkl'))
    MAR_siRNA=pickle.load(f); f.close()
    
    if sh==None:
        sh=show
    
    parameters = None
    fileL = filter(lambda x: 'distances_' in x, os.listdir(folder))
    n_siRNA = len(fileL)
    result_MAX=[np.zeros(shape=(4,n_siRNA), dtype=float) for k in range(len(featureL))]
    result=[np.zeros(shape=(4,), dtype=object) for k in range(len(featureL))]
    for i in range(len(featureL)):
        for k in range(4):
            result[i][k]=[]

    siRNAL=[]; colors=[]
    count=0
    for j, file_ in enumerate(fileL):
        siRNA = file_.split('_')[-1][:-4]
        siRNAL.append(siRNA)
        if siRNA in MAR_siRNA:
            colors.append('red')
        else:
            colors.append('green')
            
        f=open(os.path.join(folder, file_))
        distances = pickle.load(f); f.close()
        if parameters == None:
            parameters = distances.keys()
        else:
            for param in filter(lambda x: x not in parameters, distances.keys()):
                parameters.append(param)
        
        for i,feature in enumerate(featureL):
            for k,param in enumerate(parameters):
                result_MAX[i][k,j]=np.max(distances[param][:,featuresNumeriques.index(feature)])
                l=list(distances[param][:,featuresNumeriques.index(feature)])
                count+=len(l)
                result[i][k].extend(l)
                
    print count, len(result[i][0])
    for i,feature in enumerate(featureL):
        f, axes = p.subplots(len(parameters), sharex=True, sharey=True,figsize=(24,17))
        for k,param in enumerate(parameters):
            axes[k].scatter(range(len(siRNAL)), result_MAX[i][k,:], color=colors)
            axes[k].text(3,3, np.std(result_MAX[i][k]))
            try:
                axes[k].text(5,5, pearsonr(result_MAX[i][0],result_MAX[i][k])[0])
            except ValueError:
                pass
            axes[k].set_title('Distances from exp to control {}'.format(param))
        f.suptitle('{}'.format(feature))
        if sh:
            p.show()
        else:
            p.savefig(os.path.join(folder, 'distribution_MAX_{}.png'.format(feature)))
            
    for i,feature in enumerate(featureL):
        f2, axes2 = p.subplots(len(parameters), sharex=True, sharey=True,figsize=(24,17))
        for k,param in enumerate(parameters):
            n,bins,patches=axes2[k].hist(result[i][k], bins=1000, normed=False)
            axes2[k].text(3,3, np.std(result[i][k]))
            try:
                axes2[k].text(5,5, pearsonr(result[i][0],result[i][k])[0])
            except ValueError:
                pass
            axes2[k].set_title('Distribution of distances from exp to control{}'.format(param))
        f2.suptitle('{}'.format(feature))
        if sh:
            p.show()
        else:
            p.savefig(os.path.join(folder, 'distribution_{}.png'.format(feature)))
    return
    
class cellExtractor():
    def __init__(self, siRNA, settings_file,
                 testCtrl = False,
                 div_name='total_variation', bin_type = 'quantile', cost_type = 'number',
                 bin_size=50,lambda_=10,M=None, verbose=0):

        self.siRNA = siRNA
        self.settings = settings.Settings(settings_file, globals())
        self.verbose=verbose
        if not testCtrl:
            self.plate = None
        else:
            if self.verbose:
                print "Testing controls for plate {}".format(testCtrl)
            assert (testCtrl in os.listdir(self.settings.data_folder))
            self.plate = testCtrl
            self.siRNA = 'CTRL_{}'.format(self.plate[:9])
            
        self.div_name =div_name
        self.cost_type=cost_type
        self.bins_type = bin_type
        self.bin_size = bin_size
        self.lambda_=lambda_        
        
        self.currInterestFeatures = featuresNumeriques
        self.currInterestFeatures.extend(['mean persistence',  'mean straight'])
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

        _,r, _, _,_, length, _, _, _ = histConcatenation(self.settings.data_folder, self.expList, self.settings.mitocheck_file,
                                        self.settings.quality_control_file, verbose=self.verbose)
                    
        for i in range(len(length)):
            for k,feature in enumerate(self.currInterestFeatures):
                histDict[feature].append(r[np.sum(length[:i]):np.sum(length[:i+1]),featuresSaved.index(feature)])
        
        f=open(os.path.join(self.settings.result_folder, 'distExp_ctrl_{}_{}.pkl'.format(self.bins_type, self.bin_size)))
        bins = pickle.load(f); f.close()
        
        histogrammes, bins = computingBins(histDict, [self.bin_size for k in range(len(self.currInterestFeatures))], self.bins_type, previous_binning=bins)
                    
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
                _,curr_r, _, _,_, curr_length, _, _, _ = histConcatenation(self.settings.data_folder, ctrlExpList, self.settings.mitocheck_file,
                                            self.settings.quality_control_file, verbose=self.verbose)
            except:
                print "Problem with controls from plate {}".format(plate)
            else:    
                length.append(np.sum(curr_length))
                assert(np.sum(curr_length)==curr_r.shape[0])                        
                for k, feature in enumerate(self.currInterestFeatures):
                    histDict[feature].append(curr_r[:,featuresSaved.index(feature)])
        assert(len(length)==len(plates))
        assert(len(histDict[feature])==len(plates))
                    
        histogrammes, bins = computingBins(histDict, [self.bin_size for k in range(len(self.currInterestFeatures))], self.bins_type, previous_binning=bins)
        
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
            if "transportation" in self.div_name:
                if self.cost_type!='number':
                    print 'Need to set up cost matrix to match with cost_type value'
                    raise
            for k, feature in enumerate(self.currInterestFeatures):
                dist_weights = np.zeros(shape=(len(self.currInterestFeatures)+1,))
                dist_weights[k+1]=1
                mat_hist_sizes = np.array([[self.bin_size for j in range(len(self.currInterestFeatures))]])
                distances[i, k]=_distances(histogrammes[np.newaxis, i], ctrl_hist[np.newaxis,:],div_name=self.div_name, lambda_=self.lambda_, M=None,
                                     mat_hist_sizes=mat_hist_sizes, nb_feat_num=0, dist_weights=dist_weights)[0][0]
                
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
        f=open(os.path.join(self.settings.result_folder, self.settings.outputFile.format(self.siRNA)), 'w')
        pickle.dump(d, f)
        f.close()
        return
    
    def __call__(self):
        
    #i. getting experiments corresponding to siRNA if we're not looking for control experiments only
        if self.plate is None:
            self.expList=self._findExperiment()
            if self.expList==[]:
                sys.stderr.write("No experiments for siRNA {}".format(self.siRNA))
                return
            #ii.calculate the histograms and binnings of experiments, using pre-computed binnings, ON THE EXPERIMENTS ONLY
            length, histogrammes, bins = self.getData(self.settings.histDataAsWell)
        
            #iii. calculate the histograms of corresponding controls. Should be per plate then, a dictionary
            plates, ctrl_length, ctrl_histogrammes = self.ctrlHistograms(bins)
            
        else:
            
            self.expList = appendingControl([self.plate])
            if self.verbose:
                print "all ctrl", self.expList
                
            #ii.calculate the histograms and binnings of experiments, using pre-computed binnings, ON all control experiments 
            #for the plate
            length, histogrammes, bins = self.getData(self.settings.histDataAsWell)
            
            plates = [self.plate]
            #randomly selecting two wells of the plate that will be used to be compared to the others
            false_exp = np.random.randint(len(length), size=2)
            true_ctrl = filter(lambda x: x not in false_exp, range(len(length)))
            print false_exp, true_ctrl
            pdb.set_trace()
            
            ctrl_histogrammes = histogrammes[true_ctrl]
            histogrammes = histogrammes[false_exp]
        
            ctrl_length = length[true_ctrl]
            length = length[false_exp]
        
        #iv. calculate the distance from each experiment to its control
        distances = self.calculateDistances(plates, histogrammes, ctrl_histogrammes, length, ctrl_length)
        print distances
        #v. save results
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
    parser.add_option('--testCtrl', type=str, dest='testCtrl', default=0)

    parser.add_option('--div_name', type=str, dest='div_name', default='total_variation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)
    
    parser.add_option("-l",type=int, dest="lambda_", default=10)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    (options, args) = parser.parse_args()
    
    settings_file = 'tracking/settings/settings_feature_extraction.py'

    extractor=cellExtractor(options.siRNA, settings_file,options.testCtrl, options.div_name, options.bins_type,
                 bin_size=options.bin_size,lambda_=options.lambda_, verbose=options.verbose)
    extractor()
