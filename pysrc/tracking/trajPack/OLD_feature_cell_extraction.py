import getpass, os,pdb, operator
import numpy as np
import cPickle as pickle

from collections import Counter, defaultdict
from scipy.stats import scoreatpercentile

from tracking.trajPack import featuresSaved, featuresHisto
from util.listFileManagement import fromShareToCBIO
from tracking.PyPack.fHacktrack2 import initXml, finirXml
from tracking.histograms.k_means_transportation import DIVERGENCES, _distances
from tracking.histograms.transportation import costMatrix
from util.sandbox import concatCtrl
from tracking.plots import plotAlignedTraj

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
    def __init__(self, feature_name, mat_hist_sizes=np.array([[50, 50, 50, 50, 50]]), verbose=0):
        assert((feature_name in featuresSaved) or
               (feature_name in featuresHisto) )
        self.feature_name=feature_name
        self.numerical=feature_name in featuresSaved
        
        if not self.numerical:
            self.mat_hist_sizes=mat_hist_sizes
            self.dist_weights = np.zeros(shape=(mat_hist_sizes.shape[1]+1))
            self.dist_weights[1+featuresHisto.index(self.feature_name)]=1
        
        self.verbose=verbose
        
    def __call__(self, data, length, who,ctrlStatus, percentile, outputFolder='../resultData/images', how_many=5,
                div_name='total_variation', nb_feat_num=12, sample_size=20,lambda_=10,M=None):
        
        
        down, up, ctrl=self.find(data, ctrlStatus,length, percentile, how_many, div_name, nb_feat_num, sample_size,
                lambda_,M)
        all_={}
    #DEALING WITH PERCENTILES 
        wells=self.trajToWells(length, who, up)
        
        tab=[(id_, len(wells[id_])) for id_ in wells]
        sort=sorted(tab, key=operator.itemgetter(1))
        bestMovies={a[0]:wells[a[0]] for a in sort[-how_many:]}

        all_['up']=dict(bestMovies)
        
    #dealing with the other half
        wells=self.trajToWells(length, who, down)
        
        tab=[(id_, len(wells[id_])) for id_ in wells]
        sort=sorted(tab, key=operator.itemgetter(1))
        bestMovies={a[0]:wells[a[0]] for a in sort[-how_many:]}
        all_['do']=dict(bestMovies)

    #DEALING WITH MEDIAN ONES
        ctrlW=self.trajToWells(length, who, ctrl)
        tab=[(id_, len(ctrlW[id_])) for id_ in ctrlW]
        sort=sorted(tab, key=operator.itemgetter(1))
        bestMovies={a[0]:ctrlW[a[0]] for a in sort[-how_many:]}
        
        all_['med']=dict(bestMovies)
#        #fromShareToCBIO({self.feature_name:bestMovies.keys()}, wellFolder=True)
        self.wellToPlot(length, who,all_, outputFolder, median=True)
        return 1
        
    def trajToWells(self, length, who, trajectories):
        cum=np.cumsum(length)
        cum=np.insert(cum, 0, 0)
        #l=[0]; l.extend(cum); cum=np.array(l)
        result=defaultdict(list)
        for traj_index in trajectories:
            indexCour = np.where(cum<=traj_index)[0][-1]
            result[who[indexCour]].append(traj_index-cum[indexCour])
        return result
        
    def wellToPlot(self, length,who, all_, outputFolder, median=False):
        assert getpass.getuser()!='lalil0u', 'You are not running this on the right computer'
        index=featuresSaved.index(self.feature_name)
        folder='/share/data20T/mitocheck/tracking_results'
        basename = 'traj'
        resultCour=[]; sizes=[0]
        valCour=[]
        #if median: basename+='_median'
        for el in all_:
            trajectories=all_[el]
            print el
            for plate, well in trajectories:
                print plate, well
                if self.verbose:
                    print "Taking care of plate {}, well {}".format(plate, well)
        
                f=open(os.path.join(folder,plate, 'hist_tabFeatures_{}.pkl'.format(well)))
                tab, coord, _=pickle.load(f); f.close()
                resultCour.extend(np.array(coord)[trajectories[(plate, well)]])
                valCour.extend(tab[trajectories[(plate, well)], index])
            sizes.append(len(resultCour))

        XX=[]; YY=[]
        for k in range(len(resultCour)):
            X=np.array(resultCour[k][1]); X-=X[0]; XX.append(X)
            Y=np.array( resultCour[k][2]); Y-=Y[0]; YY.append(Y)
        minx,maxx = min([min(X) for X in XX]), max([max(X) for X in XX])
        miny,maxy = min([min(Y) for Y in YY]), max([max(Y) for Y in YY])
        
        #saving trajectories 
        f=open('trajectories_{}.pkl'.format(self.feature_name), 'w')
        pickle.dump([all_, resultCour, valCour, sizes],f)
        f.close()
        
        i=0
        for el in all_:
            for k in range(sizes[i], min(sizes[i]+50,sizes[i+1])):
                plotAlignedTraj(XX[k], YY[k], minx, maxx, miny, maxy, show=False, 
                                name=os.path.join(outputFolder, '{}_{}_{}_{}.png'.format(basename,el, self.feature_name, k)), val=valCour[k])
            i+=1
        return 1#_writeXml(plate, well, resultCour)
    
    def wellToXml(self, length,who, plate, well, trajectories, outputFolder, num):

        assert getpass.getuser()!='lalil0u', 'You are not running this on the right computer'
        folder='/share/data/mitocheck/tracking_results'

        if self.verbose:
            print "Taking care of plate {}, well {}".format(plate, well)

        f=open(os.path.join(folder,plate, 'hist_tabFeatures_{}.pkl'.format(well)))
        _, coord, _, _, _, _=pickle.load(f); f.close()
        resultCour=np.array(coord)[trajectories]
        return _writeXml(plate, well, resultCour, num)
        
    def find(self, data,ctrlStatus,length, percentile, how_many=5, div_name='total_variation', nb_feat_num=12, sample_size=20,
                lambda_=10,M=None):
        if self.numerical:
            return self._findNumerical(data, percentile, how_many)
        return self._findHistogram(data,ctrlStatus,length, percentile, how_many, div_name,nb_feat_num, sample_size,lambda_, M)
    
    def _findNumerical(self, data, percentile, how_many):
        index=featuresSaved.index(self.feature_name)
        scoreD = scoreatpercentile(data[:,index], percentile)
        down=np.where(data[:,index]<=scoreD)
        scoreU = scoreatpercentile(data[:,index], 100-percentile)
        up=np.where(data[:,index]>=scoreU)
            
        med0, med1=scoreatpercentile(data[:,index], 48), scoreatpercentile(data[:,index], 52)
        ctrl0 = np.where(data[:,index]>=med0)
        ctrl1=np.where(data[:, index]<=med1)
        ctrl=filter(lambda x: x in ctrl1[0], ctrl0[0])
        
        if self.verbose:
            print "Feature {}".format(self.feature_name)
            print "{} retrieved, percentile {}, value {}".format(up[0], 100-percentile, scoreU)
            print "{} retrieved, percentile {}, value {}".format(down[0],percentile, scoreD)
            print "{} retrieved, percentile {}, value {}".format(np.array(ctrl), percentile, med0)
        return down[0], up[0], ctrl
    
    def _findHistogram(self, data,ctrlStatus,length, percentile, how_many, div_name,nb_feat_num, sample_size, lambda_, M):
        #not implemented yet
        assert(div_name in DIVERGENCES)
        if div_name=='transportation' and M is None:
            M={}
            for m in range(self.mat_hist_sizes.shape[1]):
                for n in range(self.mat_hist_sizes.shape[0]):
                    M[(n,m)]=costMatrix(self.mat_hist_sizes[n,m], 1)
        
        print "CURRENT WEIGHTS", self.dist_weights
        
        random_state= np.random.mtrand._rand
        ctrlData=concatCtrl(data, ctrlStatus, length, len(ctrlStatus)-np.sum(ctrlStatus))
        seeds = random_state.permutation(ctrlData.shape[0])[:sample_size]
        sample=ctrlData[seeds]
        
        #taille sample_size, data.shape[0]
        distances = _distances(data, sample,div_name, lambda_, M, self.mat_hist_sizes, nb_feat_num, self.dist_weights)
        distances=np.argmin(distances, 0)#donc doit etre de taille data.shape[0]
        
        #then we have a numerical feature so we can get the percentile the same way
        score = scoreatpercentile(distances, percentile)
        if percentile>50:
            r=np.where(distances>=score)[0]
        else:
            r=np.where(distances<=score)[0]
        if self.verbose:
            print "{} examples retrieved for feature {}, percentile {}".format(len(r), self.feature_name, percentile)
        return r
 
