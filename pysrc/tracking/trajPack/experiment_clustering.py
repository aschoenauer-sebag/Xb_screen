import os, pdb
import cPickle as pickle

import numpy as np

from sklearn.utils import check_random_state
from operator import itemgetter
from util import settings
from tracking.histograms import *
from tracking.trajPack import featuresSaved, featuresNumeriques
from tracking.trajPack.clustering import histConcatenation
from collections import defaultdict
from tracking.histograms.k_means_transportation import _costMatrix, histogramMiniKMeans
from tracking.histograms.transportation import computingBins
from tracking.histograms.summarization_clustering import stabilityCalculation
from util.listFileManagement import countingDone, appendingControl
from optparse import OptionParser

interestFeatures= featuresNumeriques
interestFeatures.extend(['mean persistence',  'mean straight'])

'''
Here we are going to cluster experiments according to the features of the trajectories inside them, put into histograms
'''

def collectingData(iter_, expList, debut, fin):
    folder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/resultData/experiment_clustering/'
    
    histDict=defaultdict(list)
    
    _,r, _, who,ctrlStatus, length, genes, siRNAs, _ = histConcatenation('/share/data20T/mitocheck/tracking_results', expList[debut:fin], 
                    '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/mitocheck_siRNAs_target_genes_Ens72.txt', '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/qc_export.txt')
    
    for i in range(len(length)):
        for k,feature in enumerate(interestFeatures):
            histDict[feature].append(r[np.sum(length[:i]):np.sum(length[:i+1]),featuresSaved.index(feature)])
    
    f=open('../resultData/experiment_clustering/distExp_ctrl_quantile_10.pkl')
    bins = pickle.load(f); f.close()
    
    histogrammes, bins = computingBins(histDict, [10 for k in range(16)], 'quantile', previous_binning=bins)
    f=open(os.path.join(folder, 'data_{}.pkl'.format(iter_)), 'w')
    pickle.dump((histogrammes, who, ctrlStatus, genes, siRNAs),f)
    f.close()



class clusteringExperiments():
    def __init__(self, settings_file,div_name, bins_type, cost_type, bin_size,
                 init='k-means++',lambda_=10, M=None, dist_weights=None,
                 n_init=10, verbose=0, iter_=0):
        
        assert(div_name in DIVERGENCES)
            
        self.settings = settings.Settings(settings_file, globals())
        
        self.bins_type=bins_type
        self.mat_hist_sizes=np.array([[bin_size for k in range(self.settings.num_features)]]) 
        self.bin_size=bin_size
        self.div_name=div_name
        self.lambda_=lambda_
        
#Cost matrix : give the power for function transportation.costMatrix if that is what is desired as ground metric. Otherwise
#if cost_type==values, the cost matrix will contain the distance between bin centers (since it's not uniform).
        self.cost_type = cost_type
        self.cost_matrix = M 
        self.dist_weights=[1 for k in range(self.settings.num_features+1)]; self.dist_weights[0]=0
        self.verbose=verbose
        self.random_state=None
        
        self.currInterestFeatures = interestFeatures
        
        self.iter_=iter_
        
    def parameters(self, pcaParameter, with_n_cluster=False):
        r={
           'bins_type': self.bins_type,
           'bin_size':self.bin_size,
           'div_name':self.div_name,
           'lambda':self.lambda_,
           'cost_type':self.cost_type,
           'dist_weights':self.dist_weights,
           }
        if with_n_cluster:
            r.update({'n_cluster':self.n_cluster})
            
        l = sorted(zip(r.keys(), r.values()), key=itemgetter(0))
        return tuple(l)
    
    def _dataPrep(self, pcaParameter):
        histDict = defaultdict(list)
        
        ctrlExp = appendingControl(self.expList)
        ctrlExp = countingDone(ctrlExp)
        np.random.shuffle(ctrlExp)
        ctrlExp=ctrlExp[:int(0.2*len(self.expList))]
        if self.verbose:
            print ctrlExp
        self.expList.extend(ctrlExp)
        
        _,r, _, _,_, length, _, _, _ = histConcatenation(self.settings.data_folder, self.expList, self.settings.mitocheck_file,
                                       self.settings.quality_control_file, verbose=self.verbose)
        for i in range(len(length)):
            for k,feature in enumerate(self.currInterestFeatures):
                histDict[feature].append(r[np.sum(length[:i]):np.sum(length[:i+1]),featuresSaved.index(feature)])
        
        f=open(os.path.join(self.settings.result_folder, 'distExp_ctrl_{}_{}.pkl'.format(self.bins_type, self.bin_size)))
        bins = pickle.load(f); f.close()
        
        histogrammes, bins = computingBins(histDict, [self.bin_size for k in range(len(self.currInterestFeatures))], self.bins_type, previous_binning=bins)
        print histogrammes.shape
        return histogrammes, bins
    
    def _saveResults(self, centers, bins, pcaParameter):
        
        if self.settings.clustering_filename.format(self.iter_) in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, self.settings.clustering_filename.format(self.iter_)), 'r')
            d=pickle.load(f); f.close()
        else:
            d={}
        d.update({self.parameters(pcaParameter,with_n_cluster=True):[centers, self.mean, self.std, bins, self.actual_expList]})
        f=open(os.path.join(self.settings.result_folder, self.settings.clustering_filename.format(self.iter_)), 'w')
        pickle.dump(d, f); f.close()

        return
        
    def __call__(self, histogrammes=None):
        filename=self.settings.filename.format(self.iter_)
        fraction = self.settings.fraction
        num_iterations_stability = self.settings.num_iterations_stability
        self.random_state = check_random_state(self.random_state)
#i.Data preparation: putting histogram data in bins and normalizing. We do that for each of the three possible cases:
#no PCA, PCA no whitening, PCA and whitening.

        for pcaParameter in self.settings.pcaParameters:
            if histogrammes == None:
                data, bins=self._dataPrep(pcaParameter)
            else:
                data=histogrammes
            if self.cost_type=='value':
                raise ValueError
#                print 'cost type value'
#                self.cost_matrix=bins
            self.cost_matrix=_costMatrix(self.mat_hist_sizes, self.cost_matrix)
                
    #ii.Looking for the number of clusters in the data
            stability=defaultdict(list)
            for n_clusters in range(self.settings.k_min, self.settings.k_max):
                if self.verbose:
                    print 'Launching minibatch k-means, n_clusters {}, divergence {}, bin_type {}, bin_size {}, cost_type {}, init size {}, batch size {}'\
                            .format(n_clusters, self.div_name, self.bins_type, self.bin_size, self.cost_type,\
                                                    self.settings.init_size, self.settings.batch_size) 
                else:
                    print n_clusters   
                
                if num_iterations_stability==0:
                    raise ValueError        
    
                elif num_iterations_stability>0:
            #IMPLEMENTATION DU TEST POUR LA STABILITE DU CLUSTERING TEL QU'IL FIGURE DANS BEN-HUR ET AL 2002
                    for it_ in range(num_iterations_stability):
                        if self.verbose>0:
                            print 'stability iteration ', it_
                        
                            print 'ONE'
                        model1 = histogramMiniKMeans(n_clusters,self.lambda_, self.mat_hist_sizes,
                             div_name=self.div_name, M=self.cost_matrix, 
                             dist_weights=self.dist_weights, nb_feat_num = self.settings.nb_feat_num,
                             init= self.settings.init,
                             batch_size= self.settings.batch_size,
                             init_size = self.settings.init_size, verbose=self.verbose,n_init= self.settings.n_init)
                        set1= self.random_state.permutation(data.shape[0])[:int(fraction*data.shape[0])]; set1.sort()
                        model1.fit(data[set1])
                        
                        if self.verbose>0:print 'TWO'
                        model2 = histogramMiniKMeans(n_clusters,self.lambda_, self.mat_hist_sizes,
                             div_name=self.div_name, M=self.cost_matrix, 
                             dist_weights=self.dist_weights, nb_feat_num = self.settings.nb_feat_num,
                             init= self.settings.init,
                             batch_size= self.settings.batch_size,
                             init_size = self.settings.init_size, verbose=self.verbose,n_init= self.settings.n_init)
                        set2= self.random_state.permutation(data.shape[0])[:int(fraction*data.shape[0])]; set2.sort()
                        model2.fit(data[set2])
                        
                        stability[n_clusters].append(stabilityCalculation(np.array(model1.labels_), np.array(model2.labels_), set1, set2))
                        if self.verbose>0:print 'stability calculation over'
            stab_array = np.array([stability[n_clusters] for n_clusters in range(self.settings.k_min, self.settings.k_max)])
            print stab_array
            print np.array([(np.mean(stability[n_clusters]), np.std(stability[n_clusters])) for n_clusters in range(self.settings.k_min, self.settings.k_max)])
            
            arr=np.array([np.mean(stability[n_clusters])+np.std(stability[n_clusters]) for n_clusters in range(self.settings.k_min, self.settings.k_max)])
            try:
                k = self.settings.k_min+np.where(arr>0.6)[0][-1]
            except IndexError:
                k=self.settings.k_min+np.argmax(arr)
            print 'Chosen k', k
            model = histogramMiniKMeans(k,self.lambda_, self.mat_hist_sizes,
                             div_name=self.div_name, M=self.cost_matrix, 
                             dist_weights=self.dist_weights, nb_feat_num = self.settings.nb_feat_num,
                             init=self.settings.init,
                             batch_size=self.settings.batch_size,
                             init_size =self.settings.init_size, verbose=self.verbose,n_init=self.settings.n_init)
            model.fit(data)
            representatives = model.find_representatives(N=self.settings.n_representatives) 
            
        #saving results
            if not os.path.isdir(self.settings.result_folder):
                os.mkdir(self.settings.result_folder)
        
            if filename in os.listdir(self.settings.result_folder):
                f=open(os.path.join(self.settings.result_folder, filename), 'r')
                d = pickle.load(f); f.close()
                d.append([self.parameters(pcaParameter),np.array(self.expList), stab_array, representatives])
            else:
                d=[(self.parameters(pcaParameter),stab_array, representatives)]
                
    #So in the file summary_experiment.pkl, we have a list of indices in the original feature array, of representatives
    #from the experiment
            f=open(os.path.join(self.settings.result_folder, filename), 'w')
            pickle.dump(d, f)
            f.close()
        return
    
    
if __name__ == '__main__':
    
    parser = OptionParser(usage="usage: %prog [options]")    
    
    #parser.add_option('--sim', type=int, dest='simulated', default=0)

    parser.add_option('--experimentFile', type=str, dest='experimentFile', default = None)
    parser.add_option('--iter', type=int, dest='iter_', default=0)
    
    parser.add_option('--div_name', type=str, dest='div_name', default='etransportation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)
    parser.add_option("--init", dest="init", type=str,default='k-means++')
    parser.add_option("--n_init", dest="n_init", type=int,default=10)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    (options, args) = parser.parse_args()
    settings_file = 'tracking/settings/settings_exp_clustering.py'
    
    file_ = open(options.experimentFile, 'r')
    histogrammes, who, ctrlStatus, genes, siRNAs = pickle.load(file_); file_.close()
    
    model = clusteringExperiments(settings_file, 
            options.div_name, options.bins_type, options.cost_type, options.bin_size, 
                 init=options.init,
                 n_init=options.n_init, 
                 verbose=options.verbose, 
                 iter_=options.iter_)
        
    model(histogrammes)




#TO COLLECT DATA
#
#
#    parser.add_option('-d', type=int, dest='debut')
#    parser.add_option('-f', type=int, dest='fin')#all_Simpson_data.pkl
#    
#    (options, args) = parser.parse_args()
#    file_ = open(options.experimentFile, 'r')
#    experimentList, _ = pickle.load(file_); file_.close()
#    
#    collectingData(options.iter_, experimentList, options.debut, options.fin)
    