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

def computeBins(bin_size, folder="../resultData/features_on_films"):
    '''
    WARNING DON'T FORGET TO SORT LISTS OF FEATURES
    '''
    print 'Using sorted numerical features'
    
    if type(bin_size)!=list:
        bin_size=[bin_size for k in range(len(featuresNumeriques))]
    minMax=np.empty(shape=(len(featuresNumeriques), 2))
    quantiles = np.empty(shape=(len(featuresNumeriques),), dtype=object)
    
    for i,feature in enumerate(sorted(featuresNumeriques)):
        print feature
        file_ = "quantiles_{}_11718_0.pkl".format(feature)
        f=open(os.path.join(folder, file_), 'r')
        distribution = pickle.load(f)[1]
        f.close()
        
        minMax[i,0]=np.min(distribution)
        minMax[i,1]=np.max(distribution)
        quantiles[i] = [scoreatpercentile(distribution, per) for per in [k*100/float(bin_size[i]) for k in range(bin_size[i]+1)]]
            

    f=open(os.path.join(folder, 'distExp_123etctrl_minmax_{}.pkl'.format(bin_size[0])), 'w')
    pickle.dump(minMax, f); f.close()
    f=open(os.path.join(folder, 'distExp_123etctrl_quantile_{}.pkl'.format(bin_size[0])), 'w')
    pickle.dump(quantiles, f); f.close()
    return
        
    
class distributionExtractor():
    def __init__(self, expList, settings_file,slot,iter_=None,bin_type = 'quantile', bin_size=50, verbose=0):

        self.expList = np.array(expList)
        self.settings = settings.Settings(settings_file, globals())
        self.iter_=iter_
        self.slot=slot
        
        self.bin_type = bin_type
        self.bin_size = bin_size        
        self.verbose=verbose
        self.currInterestFeatures = featuresNumeriques
        if self.settings.histDataAsWell:
            raise AttributeError
    
    def getData(self, histDataAsWell):
        '''
        Ici sur toutes les experiences dans self.expList on construit l'histogramme de toutes les features numeriques
        '''
        histDict = defaultdict(list)
        r, _, _,_, length, _, _, _ = histConcatenation(self.settings.data_folder, self.expList, self.settings.mitocheck_file,
                                        self.settings.quality_control_file, verbose=self.verbose)
        for k, feature in enumerate(self.currInterestFeatures):
            for i in range(len(length)):
                histDict[feature].append(r[np.sum(length[:i]):np.sum(length[:i+1]),k])
                    
        histogrammes, bins = computingBins(histDict, [self.bin_size for k in range(len(self.currInterestFeatures))], self.bin_type, iter_=self.iter_ )
                    
        return histogrammes, bins
    
    
    def __call__(self, distribution_only):
        #i.calculate the histograms and binnings of experiments, 20% of controls
        histogrammes, bins = self.getData(self.settings.histDataAsWell)
        if distribution_only:
            return
        
        return

    
if __name__ == '__main__':
    '''
    Idea here is to calculate the distance of an experiment histogram to its control, for each numerical feature
    nb_exp should be the minimum number of experiments to get a stable distribution of the feature among
    considered experiments, so that the calculation is parallelizable
    '''
    parser = OptionParser(usage="usage: %prog [options]")    
    
    parser.add_option('--experimentFile', type=str, dest='experimentFile', default = None)
    parser.add_option('--first_exp', type=int, dest='first_exp', default =0)
    parser.add_option('--iter', type=int, default =0, dest='iter_')
    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
    parser.add_option('--bin_size', type=int, dest="bin_size", default=50)
    parser.add_option('--nb_exp', type=int, dest='size', default=1000)
    
    parser.add_option("-l",type=int, dest="lambda_", default=10)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    (options, args) = parser.parse_args()
    
    settings_file = 'tracking/settings/settings_feature_extraction.py'
    
    distribution_only=True
    
    f=open(options.experimentFile)
    expList = pickle.load(f)
    f.close()
    
    if distribution_only:
        np.random.shuffle(expList)
    
    extractor=distributionExtractor(expList[options.first_exp:options.first_exp+options.size], settings_file,options.first_exp, options.iter_,options.bins_type,
                 bin_size=options.bin_size,verbose=options.verbose)
    extractor(distribution_only)
