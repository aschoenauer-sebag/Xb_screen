import warnings, pdb, time, sys, os
import getpass

sys.path.insert(0, '/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc/') #CHAUD
import cPickle as pickle
import numpy as np

from scipy.stats import chi2
from sklearn.utils import check_random_state
from optparse import OptionParser
from collections import defaultdict, Counter
from operator import itemgetter

from transportation import multSinkhorn,multEMD1d, costMatrix, ddcostMatrix, computingBins, computingddBins
from k_means_transportation import _costMatrix, histogramMiniKMeans
from tracking.trajPack import featuresNumeriques
from tracking.trajPack.clustering import histConcatenation
from tracking.histograms import *
from tracking.histograms.k_means_transportation import _labels_inertia_precompute_dense
from util.sandbox import histLogTrsforming
from util.listDealing import expSi, appendingControl
from util.fileManagement import strToTuple
from util import settings

from rpy2 import robjects
from rpy2.robjects.packages import importr
rStats=importr("stats")

def stabilityCalculation(labels1, labels2, set1, set2):
    intersection = filter(lambda x: x in set2, set1)
    indices1=[np.where(set1==i)[0][0] for i in intersection]
    indices2=[np.where(set2==i)[0][0] for i in intersection]
    labels1=labels1[indices1]; labels2=labels2[indices2]
    
    corr1=np.array([labels1==label for label in labels1])
    corr2=np.array([labels2==label for label in labels2])
    for i in range(corr1.shape[0]):
        corr1[i,i]=0; corr2[i,i]=0
    
    return np.sum(corr1*corr2)/np.sqrt(np.sum(corr1*corr1)*np.sum(corr2*corr2))

class hitFinder():
    def __init__(self, settings_file, siRNA):            
        self.settings = settings.Settings(settings_file, globals())
        self.siRNA = siRNA
        
    def _findExperiment(self):
        '''
        Purpose : find experiments that are related to this siRNA
        '''
        yqualDict=expSi(self.settings.quality_control_file, sens=0)
        l= strToTuple(yqualDict[self.siRNA], os.listdir(self.settings.data_folder))
        assert( np.all([self.settings.summary_filename.format(el[0], el[1]) in os.listdir(self.settings.result_folder) for el in l]))
            
        return l
    
    def _findControls(self, expList):
        '''
        For a given list of experiments, returns a list of their control experiments.
        There is NOT a filter on control experiments which have been used to do the clustering for the moment
        '''
        plates = Counter(np.array(expList)[:,0]).keys()
#        result = filter(lambda x: self.settings.summary_filename.format(x[0], x[1]) in os.listdir(self.settings.result_folder), appendingControl(plates))
        result = appendingControl(plates)
        return result
    
    def _dataPrep(self, expList, ctrlList):
        '''
        Data should not be normalized as the normalization needs to be done following the original center normalization
        '''
        expList.extend(ctrlList)
        
        r, histNtot,  who,ctrlStatus, length, _, sirna, _ = histConcatenation(self.settings.data_folder, expList, self.settings.mitocheck_file,
                                        self.settings.quality_control_file)
        assert(np.all(np.array(sirna)[np.where(np.array(ctrlStatus)==1)]==self.siRNA))
        
        return r[:,:self.settings.nb_feat_num], histNtot, who, ctrlStatus, length
    
    def _attribute_labels(self, r, hist_features):
        result={}
        try:
            f=open(os.path.join(self.settings.result_folder, self.settings.clustering_filename), 'r')
            d=pickle.load(f); f.close()
        except IOError:
            raise
        else:
            for param_tuple in d:
                centers, mean, std, bins = d[param_tuple]
                curr_parameters = dict(param_tuple)
                
                mat_hist_sizes = np.array([[curr_parameters['bin_size'] for k in range(5)]])
                cost_matrix = None
                if curr_parameters['ddimensional']:
                    #histogrammeMatrix, bins = computingddBins(hist_features, self.mat_hist_sizes[0], bin_type=self.bins_type, previous_binning = bins)
                    pass
                else:
                    histogrammeMatrix, bins = computingBins(hist_features, nb_bins_list = mat_hist_sizes[0], bin_type=curr_parameters["bins_type"], previous_binning = bins)
                    
                r=(r-mean)/std
                data=np.hstack((r, histogrammeMatrix))
                if curr_parameters['ddimensional']:
                    pass
#                    self.cost_matrix={(0,0):ddcostMatrix(np.product(self.mat_hist_sizes[0]), self.mat_hist_sizes[0])}
#                    self.mat_hist_sizes=np.array([[np.product(self.mat_hist_sizes)]])
                    #filename = 'dd_'+self.settings.clustering_filename
                else:   
                    if curr_parameters["cost_type"]=='value':
                        print 'cost type value'
                        cost_matrix=bins
                    cost_matrix=_costMatrix(mat_hist_sizes, cost_matrix)
                
                labels, _, _ = _labels_inertia_precompute_dense(data, 
                                                                curr_parameters['div_name'], curr_parameters["lambda_"], 
                                                                M=cost_matrix, 
                                                                mat_hist_sizes=mat_hist_sizes, 
                                                                nb_feat_num=self.settings.nb_feat_num, 
                                                                dist_weights=curr_parameters["dist_weights"], centers=centers)
                result[param_tuple] = labels
            
            return result
        
    def _computePValues(self, labelDict, ctrlList, who, ctrlStatus, length):
        info = np.array( zip((who, ctrlStatus, length)))
        p_vals = {}
        for param_tuple in labelDict:
            labels = labelDict[param_tuple]
            p_vals[param_tuple] = defaultdict(list)
            
            for experiment in info[np.where(info[1]==1), 0]:
                i=who.index(experiment)
                pLabelsCour = labels[np.sum(length[:i]):np.sum(length[:i+1])]
                ctrlPl = filter(lambda x: x[0]==experiment[0], ctrlList)
                cLabelsCour =[]
                for ctrlel in ctrlPl:
                    try:
                        index = who.index(ctrlel)
                    except ValueError:
                        continue
                    else:
                        cLabelsCour.extend(labels[np.sum(length[:index]):np.sum(length[:index+1])])
                
                vecLongueurs = [0 for k in range(len(cLabelsCour))]; vecLongueurs.extend([1 for k in range(len(pLabelsCour))])
                cLabelsCour.extend(pLabelsCour)
                try:
                    p_vals[param_tuple].append(np.float64(rStats.fisher_test(robjects.IntVector(cLabelsCour), robjects.IntVector(vecLongueurs), 
                                                                             workspace=200000000, simulate_p_value=True)[0])[0])
                except:
                    pdb.set_trace()
        #statistical test according to Fisher's method http://en.wikipedia.org/wiki/Fisher%27s_method
            stat = -2*np.sum(np.log(p_vals[param_tuple]))
            chi2_dist = chi2(df = 2*len(p_vals[param_tuple]))
        #Given that what we're checking for is small p-values = big stat values, we can limit ourselves to computing
            #right-tail p-values
            p_vals[param_tuple] = chi2_dist.sf(stat)
        return p_vals
    
    def _saveResults(self, p_values):
        if self.settings.hit_filename in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, self.settings.hit_filename), 'r')
            hits = pickle.load(f); f.close()
        else:
            hits={}
            
        if self.siRNA in hits:
            hits[self.siRNA].update(p_values)
        else:
            hits[self.siRNA]=p_values
            
        f=open(os.path.join(self.settings.result_folder, self.settings.hit_filename), 'w')
        pickle.dump(hits, f); f.close()
    
    def __call__(self):
        #i. getting experiments corresponding to siRNA
        expList=self._findExperiment()
        
        #ii.getting controls corresponding to experiments
        ctrlList = self._findControls(expList)
        
        #iii. get data
        r, histNtot, who, ctrlStatus, length = self._dataPrep(expList, ctrlList) 
        
        #iv. find labels for data
        labelDict = self._attribute_labels(r, histNtot)
        
        #v. compute p-values: are experiment trajectories significantly differently clustered than control trajectories
        pvalues = self._computePValues(labelDict, ctrlList, who, ctrlStatus, length)
        
        #vi. save results
        self._saveResults(pvalues)
        
        return


class clusteringExperiments():
    def __init__(self, settings_file, experimentList,n_cluster, div_name, bins_type, cost_type, bin_size, 
                 init='k-means++',lambda_=10, M=None, dist_weights=None, batch_size=100, init_size=1000,
                 n_init=5, verbose=0, ddim = 0):
        
        assert(div_name in DIVERGENCES)
            
        self.settings = settings.Settings(settings_file, globals())
        self.expList = experimentList
        
        self.n_cluster = n_cluster
        self.bins_type=bins_type
        self.mat_hist_sizes=np.array([[bin_size for k in range(5)]]) if ddim==0 else np.array([[bin_size for k in range(2)]])
        self.bin_size=bin_size
        self.init=init
        self.div_name=div_name
        self.lambda_=lambda_
        
#Cost matrix : give the power for function transportation.costMatrix if that is what is desired as ground metric. Otherwise
#if cost_type==values, the cost matrix will contain the distance between bin centers (since it's not uniform).
        self.cost_type = cost_type
        self.cost_matrix = M 
        self.dist_weights=dist_weights
        self.batch_size = batch_size
        self.init_size=init_size
        self.n_init=n_init
        self.verbose=verbose
        self.random_state=None
        
        self.ddimensional = ddim
        
    def parameters(self, with_n_cluster=False):
        r={
           'bins_type': self.bins_type,
           'bin_size':self.bin_size,
           'div_name':self.div_name,
           'lambda':self.lambda_,
           'cost_type':self.cost_type,
           'dist_weights':self.dist_weights,
           'ddimensional':self.ddimensional
           }
        if with_n_cluster:
            r.update({'n_cluster':self.n_cluster})
            
        l = sorted(zip(r.keys(), r.values()), key=itemgetter(0))
        return tuple(l)
    
    def _dataPrep(self):
        num_features=None
        hist_features = {hist_feature : [] for hist_feature in featuresHisto}

        count=0
        for i, experiment in enumerate(self.expList):
            print '------{}/{} experiments'.format(i, len(self.expList))
            try:
                f=open(os.path.join(self.settings.result_folder, self.settings.summary_filename.format(experiment[0], experiment[1])))
                representatives = pickle.load(f); f.close()
            except IOError:
                sys.stderr.write("Representative trajectories for experiment {} have not been calculated at all.".format(experiment))

            else:
                
                r, histNtot,  _,_, _, _, _, _ = histConcatenation(self.settings.data_folder, [experiment], self.settings.mitocheck_file,
                                        self.settings.quality_control_file, verbose=self.verbose)
                
                curr_parameters = self.parameters(with_n_cluster=False)
                try:
                    num_features=r[representatives[curr_parameters], :self.settings.nb_feat_num] if num_features is None else \
                    np.vstack((num_features, r[representatives[curr_parameters], :self.settings.nb_feat_num]))
                except KeyError:
                    sys.stderr.write("Representative trajectories for experiment {}, parameters {} have not yet been calculated.".format(experiment, curr_parameters))
                else:                    
                    for hist_feature in hist_features:
                        hist_features[hist_feature].extend([histNtot[hist_feature][k] for k in representatives[curr_parameters]])
                    count+=len(representatives[curr_parameters]) 
                    
        if self.verbose:     
            print num_features.shape, 'not normalized'
            print 'computing bins, bin type {}, d-dimensional? {}'.format(self.bins_type, bool(self.ddimensional))
        assert(count == num_features.shape[0])
            
        if self.ddimensional:
            histogrammeMatrix, bins = computingddBins(hist_features, self.mat_hist_sizes[0], bin_type=self.bins_type)
        else:
            histogrammeMatrix, bins = computingBins(hist_features, self.mat_hist_sizes[0], bin_type=self.bins_type)
            
        self.mean = np.mean(num_features,0)
        self.std = np.std(num_features,0)
        num_features=(num_features-self.mean)/self.std
        data=np.hstack((num_features, histogrammeMatrix))
        
        return data, bins     
    
    def _saveResults(self, centers, bins):
        
        if self.settings.clustering_filename in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, self.settings.clustering_filename), 'r')
            d=pickle.load(f); f.close()
        else:
            d={}
        d.update({self.parameters():[centers, self.mean, self.std, bins]})
        f=open(os.path.join(self.settings.result_folder, self.settings.clustering_filename), 'w')
        pickle.dump(d, f); f.close()

        return
        
    def __call__(self):
        dist_weights = WEIGHTS[self.dist_weights] if not self.ddimensional else ddWEIGHTS[self.dist_weights]
        
        #i. assembling all data 
        data, bins=self._dataPrep()
        if self.ddimensional:
            self.cost_matrix={(0,0):ddcostMatrix(np.product(self.mat_hist_sizes[0]), self.mat_hist_sizes[0])}
            self.mat_hist_sizes=np.array([[np.product(self.mat_hist_sizes)]])
            #filename = 'dd_'+self.settings.clustering_filename
        else:   
            if self.cost_type=='value':
                print 'cost type value'
                self.cost_matrix=bins
            self.cost_matrix=_costMatrix(self.mat_hist_sizes, self.cost_matrix)
        
        #ii.clustering the data
        model = histogramMiniKMeans(self.n_cluster,self.lambda_, self.mat_hist_sizes,
                         div_name=self.div_name, M=self.cost_matrix, 
                         dist_weights=dist_weights, nb_feat_num = self.settings.nb_feat_num,
                         init=self.init,
                         batch_size=self.batch_size,
                         init_size =self.init_size, verbose=self.verbose,n_init=self.n_init)
        model.fit(data)
        
        #iii.saving the centers, mean and std to be further used movie per movie
        self._saveResults(model.cluster_centers_, bins)        
        return

class summarizingExperiment():
    '''
    The purpose of this class is to take an experiment, find the right number of clusters in its trajectory
    set, and summarize it with N trajectories per cluster.
    '''
    
    def __init__(self,settings_file, experiment,div_name, bins_type, cost_type, bin_size, 
                 init='k-means++',lambda_=10, M=None, dist_weights=None, batch_size=1000, init_size=1000,
                 n_init=5, verbose=0, ddim = 0, random_state=None,max_no_improvement=None):
        
        assert(div_name in DIVERGENCES)
            
        self.settings = settings.Settings(settings_file, globals())
        
        self.experiment = experiment
        self.bins_type=bins_type
        self.mat_hist_sizes=np.array([[bin_size for k in range(5)]]) if ddim==0 else np.array([[bin_size for k in range(2)]])
        self.bin_size=bin_size
        self.init=init
        self.div_name=div_name
        self.lambda_=lambda_
        
#Cost matrix : give the power for function transportation.costMatrix if that is what is desired as ground metric. Otherwise
#if cost_type==values, the cost matrix will contain the distance between bin centers (since it's not uniform).
        self.cost_type = cost_type
        self.cost_matrix = M 
        self.dist_weights=dist_weights
        self.batch_size = batch_size
        self.init_size=init_size
        self.n_init=n_init
        self.verbose=verbose
        self.random_state=random_state
        
        self.ddimensional = ddim
        
    def _dataPrep(self):
        r, hist, _,_, _, _,_, _ = histConcatenation(self.settings.data_folder, [self.experiment], self.settings.mitocheck_file,
                                        self.settings.quality_control_file, verbose=self.verbose) 
        r=r[:,:len(featuresNumeriques)]
        self.nb_feat_num = len(featuresNumeriques)
        
        if self.verbose:     
            print r.shape, 'not normalized'
            print 'computing bins, bin type {}, d-dimensional? {}'.format(self.bins_type, bool(self.ddimensional))
            
        if self.ddimensional:
            histogrammeMatrix, bins = computingddBins(hist, self.mat_hist_sizes[0], bin_type=self.bins_type)
        else:
            histogrammeMatrix, bins = computingBins(hist, self.mat_hist_sizes[0], bin_type=self.bins_type)

        r=(r-np.mean(r,0))/np.std(r,0)
        r=np.hstack((r, histogrammeMatrix))
        
        return r, bins
    
    def parameters(self):
        '''
        The number of clusters that is finally chosen for the experiment is not a parameter per se
        as it comes from the output of the clustering algorithm
        '''
        
        r={'bins_type': self.bins_type,
           'bin_size':self.bin_size,
           'div_name':self.div_name,
           'lambda':self.lambda_,
           'cost_type':self.cost_type,
           'dist_weights':self.dist_weights,
           'ddimensional':self.ddimensional
           }
        
        l = sorted(zip(r.keys(), r.values()), key=itemgetter(0))
        return tuple(l)
    
    def __call__(self):
        filename=self.settings.summary_filename.format(self.experiment[0], self.experiment[1])
        fraction = self.settings.fraction
        num_iterations_stability = self.settings.num_iterations_stability
        dist_weights = WEIGHTS[self.dist_weights] if not self.ddimensional else ddWEIGHTS[self.dist_weights]
        self.random_state = check_random_state(self.random_state)
#i.Data preparation: putting histogram data in bins and normalizing
        
        #if self.cost_type==value then returning the bins for further ground metric calculation
        data, bins=self._dataPrep()
        
        if self.ddimensional:
            raise ValueError
#            self.cost_matrix={(0,0):ddcostMatrix(np.product(self.mat_hist_sizes[0]), self.mat_hist_sizes[0])}
#            self.mat_hist_sizes=np.array([[np.product(self.mat_hist_sizes)]])
#            filename = 'dd_'+filename 
        else:   
            if self.cost_type=='value':
                print 'cost type value'
                self.cost_matrix=bins
            self.cost_matrix=_costMatrix(self.mat_hist_sizes, self.cost_matrix)
            
#ii.Looking for the number of clusters in the data
        stability=defaultdict(list)
        for n_clusters in range(self.settings.k_min, self.settings.k_max):
            if self.verbose:
                print 'Launching minibatch k-means, n_clusters {}, divergence {}, bin_type {}, bin_size {}, cost_type {}, init size {}, batch size {}'\
                        .format(n_clusters, self.div_name, self.bins_type, self.bin_size, self.cost_type,\
                                                self.init_size, self.batch_size) 
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
                         dist_weights=dist_weights, nb_feat_num = self.nb_feat_num,
                         init=self.init,
                         batch_size=self.batch_size,
                         init_size =self.init_size, verbose=self.verbose,n_init=self.n_init)
                    set1= self.random_state.permutation(data.shape[0])[:int(fraction*data.shape[0])]; set1.sort()
                    model1.fit(data[set1])
                    
                    if self.verbose>0:print 'TWO'
                    model2 = histogramMiniKMeans(n_clusters,self.lambda_, self.mat_hist_sizes,
                         div_name=self.div_name, M=self.cost_matrix, 
                         dist_weights=dist_weights, nb_feat_num = self.nb_feat_num,
                         init=self.init,
                         batch_size=self.batch_size,
                         init_size =self.init_size, verbose=self.verbose,n_init=self.n_init)
                    set2= self.random_state.permutation(data.shape[0])[:int(fraction*data.shape[0])]; set2.sort()
                    model2.fit(data[set2])
                    
                    stability[n_clusters].append(stabilityCalculation(np.array(model1.labels_), np.array(model2.labels_), set1, set2))
                    if self.verbose>0:print 'stability calculation over'
        print np.array([stability[n_clusters] for n_clusters in range(self.settings.k_min, self.settings.k_max)])
        print np.array([(np.mean(stability[n_clusters]), np.std(stability[n_clusters])) for n_clusters in range(self.settings.k_min, self.settings.k_max)])
        
        arr=np.array([np.mean(stability[n_clusters])+np.std(stability[n_clusters]) for n_clusters in range(self.settings.k_min, self.settings.k_max)])
        try:
            k = self.settings.k_min+np.where(arr>0.6)[0][-1]
        except IndexError:
            k=self.settings.k_min+np.argmax(arr)
        print 'Chosen k', k
        model = histogramMiniKMeans(k,self.lambda_, self.mat_hist_sizes,
                         div_name=self.div_name, M=self.cost_matrix, 
                         dist_weights=dist_weights, nb_feat_num = self.nb_feat_num,
                         init=self.init,
                         batch_size=self.batch_size,
                         init_size =self.init_size, verbose=self.verbose,n_init=self.n_init)
        model.fit(data)
        representatives = model.find_representatives(N=self.settings.n_representatives) 
        
    #saving results
        if not os.path.isdir(self.settings.result_folder):
            os.mkdir(self.settings.result_folder)
    
        if filename in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, filename), 'r')
            d = pickle.load(f); f.close()
            d.update({self.parameters():representatives})
        else:
            d={self.parameters():representatives}
            
#So in the file summary_experiment.pkl, we have a list of indices in the original feature array, of representatives
#from the experiment
        f=open(os.path.join(self.settings.result_folder, filename), 'w')
        pickle.dump(d, f)
        f.close()
        return
    
if __name__ == '__main__':
    
    parser = OptionParser(usage="usage: %prog [options]")    
    
    #parser.add_option('--sim', type=int, dest='simulated', default=0)
    
    parser.add_option('--plate', type=str, dest='plate', default=None)
    parser.add_option('--well', type=str, dest='well', default=None)
    parser.add_option('--experimentFile', type=str, dest='experimentFile', default = None)
    parser.add_option('--siRNA', type=str, dest='siRNA', default=None)
    
    parser.add_option('-a', type=str, dest='action',default = 'summary')
    
    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)
    parser.add_option('--ddimensional', type=int, dest='ddimensional', default=0)
    parser.add_option("-k", type=int, dest="n_cluster")
    parser.add_option("-w", type=int, dest="weights", default=0)
    parser.add_option("-l",type=int, dest="lambda_", default=10)
    parser.add_option("--batch_size", dest="batch_size", type=int,default=1000)
    parser.add_option("--init", dest="init", type=str,default='k-means++')
    parser.add_option("--n_init", dest="n_init", type=int,default=5)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    (options, args) = parser.parse_args()
    if getpass.getuser()=='lalil0u':
        settings_file = 'tracking/histograms/settings/settings_Trulove.py'
    else:
        settings_file = 'tracking/histograms/settings/settings_Thalassa.py'
#    if options.simulated:
#        datafile = 'hist_tabFeatures_WO.pkl'
#        folder = '../resultData/simulated_traj'
#    else:
#        datafile = 'MAR_histdata_PCAtrsf.pkl'
#        folder='../resultData/sinkhornKMeans'
#    sys.path.append('/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc/')
    if options.action == 'summary':
        model = summarizingExperiment(settings_file, (options.plate, options.well), 
                options.div_name, options.bins_type, options.cost_type, options.bin_size, 
                  init=options.init,lambda_=options.lambda_, dist_weights=options.weights, batch_size=options.batch_size,
                 n_init=options.n_init, verbose=options.verbose, ddim = options.ddimensional)
        model()
        
    elif options.action=='clustering':
        file_ = open(options.experimentFile, 'r')
        experimentList = pickle.load(file_); file_.close()
        
        model = clusteringExperiments(settings_file, experimentList,options.n_cluster, 
                options.div_name, options.bins_type, options.cost_type, options.bin_size, 
                 init=options.init,lambda_=options.lambda_, dist_weights=options.weights, batch_size=options.batch_size,
                 n_init=options.n_init, verbose=options.verbose, ddim = options.ddimensional)
        
        model()
    
    elif options.action =='hitFinder':
        pass
        
    
                 
                 
                 
                 
                 
                 
