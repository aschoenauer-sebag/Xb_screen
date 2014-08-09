import warnings, pdb, time, sys, os
import getpass
from scipy.stats.stats import ks_2samp, pearsonr, scoreatpercentile, spearmanr

if getpass.getuser()=='aschoenauer':
    import matplotlib as mpl
    mpl.use('Agg')
    show = False
elif getpass.getuser()=='lalil0u':
    import matplotlib as mpl
    show=True
from matplotlib import pyplot as p
import brewer2mpl

import cPickle as pickle
import numpy as np

from scipy.stats import chi2
from sklearn.utils import check_random_state
from optparse import OptionParser
from collections import defaultdict, Counter
from operator import itemgetter
from itertools import product

from transportation import multSinkhorn,multEMD1d, costMatrix, ddcostMatrix, computingBins, computingddBins
from k_means_transportation import _costMatrix, histogramMiniKMeans
from tracking.trajPack import featuresNumeriques
from tracking.trajPack.clustering import histConcatenation
from tracking.plots import plotAlignedTraj
from tracking.histograms import *
from tracking.histograms.k_means_transportation import _labels_inertia_precompute_dense
from util.listFileManagement import strToTuple, expSi, appendingControl, is_ctrl, siEntrez
from util import settings
from util.plots import couleurs

#code pour plotter la dist des p-val ctrl
#import pylab as pl
#for k in range(result.shape[1]):
#    pl.hist(result[0,k, :], bins=np.logspace(-12, 0,num=100)); pl.ylim(0,100)
#    pl.gca().set_xscale("log")
#    pl.savefig('dist_pvalCtrl_{}.png'.format(k))
#    pl.close()


def _ctrlPca(r, whiten,nb_features_num,nb_composantes, folder, filename='pca_ctrlSimpson_whiten{}.pkl'):
    print "PCA, whiten ", whiten
    f=open(os.path.join(folder, filename.format(whiten)))
    pca, mean, std = pickle.load(f)
    f.close()
    
    r=(r-mean[:nb_features_num])/std[:nb_features_num]
    return pca.transform(r)[:,:nb_composantes]

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
    def __init__(self, settings_file, siRNA, verbose=0,iter_=0, testCtrl = False):
        '''
        This class is supposed to look at experiments for a given siRNA, and compare
        the repartition of its trajectories in clusters with the repartition of control trajectories
        in the same clusters. However, if testCtrl is a plate, then this class will select two CONTROL
        wells at random on the plate and do the exact same procedure but for the control wells.
        '''            
        self.settings = settings.Settings(settings_file, globals())
        self.siRNA = siRNA
        self.verbose=verbose
        self.iter_=iter_
        if not testCtrl:
            self.plate = None
        else:
            if self.verbose:
                print "Testing controls for plate {}".format(testCtrl)
            assert (testCtrl in os.listdir(self.settings.data_folder))
            self.plate = testCtrl
            self.siRNA = 'CTRL_{}'.format(self.plate[:9])
        
    def _findExperiment(self):
        '''
        Purpose : find experiments that are related to this siRNA
        '''
        yqualDict=expSi(self.settings.quality_control_file, sens=0)
        l= strToTuple(yqualDict[self.siRNA], os.listdir(self.settings.data_folder))
        if self.verbose:
            print l
        assert( np.all([self.settings.summary_filename.format(self.iter_, el[0], el[1]) in os.listdir(self.settings.result_folder) for el in l]))
            
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
        
        _, r, histNtot,  who,ctrlStatus, length, _, sirna, _ = histConcatenation(self.settings.data_folder, expList, self.settings.mitocheck_file,
                                        self.settings.quality_control_file)
        if self.plate is None:
            assert(np.all(np.array(sirna)[np.where(np.array(ctrlStatus)==1)]==self.siRNA))
        else:
            assert np.all(np.array(ctrlStatus)==0)
            ctrlStatus[0]=1
        
        return r[:,:self.settings.nb_feat_num], histNtot, who, ctrlStatus, length
    
    def _attribute_labels(self, r, hist_features,ctrlIter, param_tuple=None):
        result={}
        try:
            f=open(os.path.join(self.settings.result_folder, self.settings.clustering_filename.format(self.iter_)), 'r')
            d=pickle.load(f); f.close()
        except IOError:
            raise
        else:
            if param_tuple is not None:
                small_param_set = [param_tuple]                
            elif not self.settings.redo:
                small_param_set=self._cleanParameterSetsFromDoneWork(d, ctrlIter)
            
            for param_tuple in small_param_set:
                centers, mean, std, bins, _ = d[param_tuple]
                curr_parameters = dict(param_tuple)
                if self.verbose:
                    print curr_parameters
                mat_hist_sizes = np.array([[curr_parameters['bin_size'] for k in range(5)]])
                cost_matrix = None
                if curr_parameters['ddimensional']:
                    #histogrammeMatrix, bins = computingddBins(hist_features, self.mat_hist_sizes[0], bin_type=self.bins_type, previous_binning = bins)
                    pass
                else:
                    histogrammeMatrix, bins = computingBins(hist_features, nb_bins_list = mat_hist_sizes[0], bin_type=curr_parameters["bins_type"], previous_binning = bins)
                if curr_parameters['PCA']==0:
                    norm_r=(r-mean)/std
                else:
                    norm_r = _ctrlPca(r, curr_parameters['whiten'], self.settings.nb_feat_num, self.settings.nb_composantes, self.settings.result_folder)
                data=np.hstack((norm_r, histogrammeMatrix))
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
                    weights = WEIGHTS[curr_parameters["dist_weights"]]
                
                labels, _, _ = _labels_inertia_precompute_dense(data, 
                                                                curr_parameters['div_name'], curr_parameters["lambda"], 
                                                                M=cost_matrix, 
                                                                mat_hist_sizes=mat_hist_sizes, 
                                                                nb_feat_num=self.settings.nb_feat_num, 
                                                                dist_weights=weights, centers=centers)
                result[param_tuple] = labels

            return result
        
    def _computePValues(self, labelDict, ctrlList, who, ctrlStatus, length):
        info = np.array( zip((who, ctrlStatus, length)))[:,0]
        p_vals = defaultdict(list); labelsTot = {}
#        r=[]
        for param_tuple in labelDict:
            labels = labelDict[param_tuple]
            d=dict(param_tuple)
            labelsTot[param_tuple]={}
            
            for experiment in info[0,np.where(info[1]==1)[0]]:
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
                min_ = d['n_cluster']
                if self.verbose:
                    print np.bincount(cLabelsCour, minlength=min_)
                    print np.bincount(pLabelsCour, minlength=min_)
                    
                dist_ = np.sum( np.absolute(
                                            np.array(np.bincount(cLabelsCour, minlength=min_)/float(len(cLabelsCour)), dtype=float) 
                                                     -
                                            np.array(np.bincount(pLabelsCour,minlength=min_)/float(len(pLabelsCour)), dtype=float)
                                            ))
                if self.verbose:
                    print dist_
                labelsTot[param_tuple][experiment] = [cLabelsCour, pLabelsCour, dist_]
#
#                llength = len(cLabelsCour)
#                vecLongueurs = [0 for k in range(len(cLabelsCour))]; vecLongueurs.extend([1 for k in range(len(pLabelsCour))])
#                cLabelsCour.extend(pLabelsCour)
                #r.append([cLabelsCour, vecLongueurs])
#                if np.all(np.array(cLabelsCour)==np.zeros(len(cLabelsCour))):
#                    p_vals[param_tuple].append(1.0)
#                    continue
#                
#                p_vals[param_tuple].append(np.float64(rStats.fisher_test(IntVector(cLabelsCour), IntVector(vecLongueurs), 
#                                                                         simulate_p_value=True, B=2000000)[0][0]))
                #
#                if self.verbose:
#                    print '---------------weights', d['dist_weights'], 'k', d['n_cluster'], 'bins_type', d['bins_type']
#                    print len(pLabelsCour), np.bincount(pLabelsCour)/float(len(pLabelsCour))
#                    print llength, np.bincount(cLabelsCour[:llength])/float(llength)
#                    print p_vals[param_tuple][-1]
        #statistical test according to Fisher's method http://en.wikipedia.org/wiki/Fisher%27s_method
#            if len(p_vals[param_tuple])>1:
#                stat = -2*np.sum(np.log(p_vals[param_tuple]))
#                p_vals[param_tuple] = chi2.sf(stat, 2*len(p_vals[param_tuple]))
#
        return labelsTot
    
    def _cleanParameterSetsFromDoneWork(self, parameter_set_list, ctrlIter):
        '''
        The purpose of this function is to return a parameter_set_list that is cleaned from
        parameter sets for which the p-values have already be calculated, if settings.redo is not True
        '''
        if self.plate is not None:
            siRNA = self._giveFilename(ctrlIter)
        else:
            siRNA = self.siRNA
            
        if self.settings.pval_filename.format(self.iter_,siRNA) in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, self.settings.pval_filename.format(self.iter_, siRNA)), 'r')
            hits = pickle.load(f); f.close()
            
            return filter(lambda x: x not in hits, parameter_set_list)
        else:
            return parameter_set_list
            
        
    
    def _giveFilename(self, ctrlIter):
        if ctrlIter==0:
            return '{}{}'.format(ctrlIter, self.siRNA)
        else:
            return '{}{}'.format(ctrlIter, self.siRNA[1:])
    
    def _saveResults(self, p_values, ctrlIter):
        if self.plate is not None:
            self.siRNA = self._giveFilename(ctrlIter)
        if self.settings.pval_filename.format(self.iter_, self.siRNA) in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, self.settings.pval_filename.format(self.iter_, self.siRNA)), 'r')
            hits = pickle.load(f); f.close()
        else:
            hits={}
            
        hits.update(p_values)
            
        f=open(os.path.join(self.settings.result_folder, self.settings.pval_filename.format(self.iter_, self.siRNA)), 'w')
        pickle.dump(hits, f); f.close()
        
    def _saveLabels(self, labelDict, ctrlIter):
        if self.plate is not None:
            self.siRNA = self._giveFilename(ctrlIter)
        if self.settings.label_filename.format(self.iter_, self.siRNA) in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, self.settings.label_filename.format(self.iter_, self.siRNA)), 'r')
            labels = pickle.load(f); f.close()
        else:
            labels={}
            
        labels.update(labelDict)
            
        f=open(os.path.join(self.settings.result_folder, self.settings.label_filename.format(self.iter_, self.siRNA)), 'w')
        pickle.dump(labels, f); f.close()
        
    def plotTrajectories(self, expList, labels):
        resultCour=[]
        clusters = Counter(labels).keys()
        for plate, well in expList:
            if self.verbose:
                print "Taking care of plate {}, well {}".format(plate, well)
        
            f=open(os.path.join(self.settings.data_folder,plate, 'hist2_tabFeatures_{}.pkl'.format(well)))
            _, coord, _=pickle.load(f); f.close()
            resultCour.extend(coord)

        XX=[]; YY=[]
        for k in range(len(resultCour)):
            X=np.array(resultCour[k][1]); X-=X[0]; XX.append(X)
            Y=np.array( resultCour[k][2]); Y-=Y[0]; YY.append(Y)
        minx,maxx = min([min(X) for X in XX]), max([max(X) for X in XX])
        miny,maxy = min([min(Y) for Y in YY]), max([max(Y) for Y in YY])
        
        for cluster in clusters:
            for i in np.where(np.array(labels)==cluster)[0][:10]:
                if not os.path.isdir(os.path.join(self.settings.result_folder, 'images')): 
                    os.mkdir(os.path.join(self.settings.result_folder, 'images'))
                plotAlignedTraj(XX[i], YY[i], minx, maxx, miny, maxy, show=False, 
                     name=os.path.join(self.settings.result_folder, 'images', 'cluster{}_{}_{}.png'.format(cluster, self.siRNA,i)))

        return 1
    
    def labelOnly(self, param_tuple,plot=True, ctrlIter=0):
        
        #i. getting experiments corresponding to siRNA if we're not looking for control experiments only
        if self.plate is None:
            expList=self._findExperiment()
            ctrlList = []
            if expList==[]:
                sys.stderr.write("No experiments for siRNA {}".format(self.siRNA))
                return
        else:
            #not returning any control labels if it's not only control labels
            expList = [(self.plate, '00001_01')]
            if self.verbose:
                print "ctrl", expList
            ctrlList = self._findControls(expList)
            if self.plate is not None:
                np.random.shuffle(ctrlList)
                expList = ctrlList[:2]
                ctrlList = ctrlList[2:]
                print 'Picked ctrl experiments', expList
        
        #iii. get data
        r, histNtot, who, ctrlStatus, length = self._dataPrep(expList, ctrlList)
        
        #iv. find labels for data
        labelDict = self._attribute_labels(r, histNtot,ctrlIter, param_tuple)
        
        #vi. save labels
        self._saveLabels(labelDict, ctrlIter)
        print self.siRNA, np.bincount(labelDict[param_tuple])
        if self.plate is not None and ctrlIter==0:
            print 'going for another round'
            self.__call__(ctrlIter=1)
        if plot:
            self.plotTrajectories(expList, labelDict[param_tuple])
        return labelDict[param_tuple]
    
    def __call__(self, ctrlIter=0):
        
        #i. getting experiments corresponding to siRNA if we're not looking for control experiments only
        if self.plate is None:
            expList=self._findExperiment()
            if expList==[]:
                sys.stderr.write("No experiments for siRNA {}".format(self.siRNA))
                return
        else:
            expList = [(self.plate, '00001_01')]
            if self.verbose:
                print "ctrl", expList
        
        #ii.getting controls corresponding to experiments if not looking for ctrl experiments only,
            #else getting all controls on self.plate
        ctrlList = self._findControls(expList)
        if self.plate is not None:
            np.random.shuffle(ctrlList)
            expList = ctrlList[:1]
            ctrlList = ctrlList[1:]
            print 'Picked ctrl experiments', expList
        
        #iii. get data
        r, histNtot, who, ctrlStatus, length = self._dataPrep(expList, ctrlList)
        
        #iv. find labels for data
        labelDict = self._attribute_labels(r, histNtot, ctrlIter)
        
        #v. compute p-values: are experiment trajectories significantly differently clustered than control trajectories
        result = self._computePValues(labelDict, ctrlList, who, ctrlStatus, length)
        
        #vi. save results
        self._saveResults(result, ctrlIter)
        if self.plate is not None and ctrlIter==0:
            print 'going for another round'
            self.__call__(ctrlIter=1)
        return
    
    def _correlationParamParam(self, comparisons, parameters, result, testCtrl,iterations, sh,test=spearmanr, save=True):
        '''
        One should be aware that scoreatpercentile is influenced by nan values. 
        Indeed scoreatpercentile([0,1,0,1,np.nan]) is 1 when scoreatpercentile([0,1,0,1]) is 0.5
        '''

        parameter_corr = np.zeros(shape=(len(comparisons), result.shape[1], result.shape[1]))
        
        for k, comparison in enumerate(comparisons):
            msg=''; thresholds=[]
            msg+='Value of ctrl p-values\n'
            for j in range(result.shape[1]):
                parameter_corr[k,j,j]=1    
                try:
                    thresholds.append(scoreatpercentile(result[k,j], 10))
                except:
                    pdb.set_trace()
                msg+='{} {} \n'.format(j, thresholds[-1])

                for l in range(j+1, result.shape[1]):
                    parameter_corr[k,j,l]=test(result[k,j], result[k,l])[0]
                    parameter_corr[k,l,j]=parameter_corr[k,j,l]
            if testCtrl and save:
                print msg
                for iter_ in range(len(iterations)):
                    dict_thresholds = dict(zip(parameters, thresholds[iter_::len(iterations)]))    
                    f=open(os.path.join(self.settings.result_folder, self.settings.threshold_filename.format(comparison, iter_)), 'w')
                    pickle.dump((dict_thresholds, result[k,:]), f); f.close()
            elif not testCtrl:
                thresholds=np.zeros(shape=(len(parameters)*len(iterations)))
                for iter_ in range(len(iterations)):
                    f=open(os.path.join(self.settings.result_folder, self.settings.threshold_filename.format(comparison, iter_)), 'r')
                    dict_thresholds=pickle.load(f)[0]; f.close()
                    thresholds[iter_::len(iterations)]=np.array([dict_thresholds[el] for el in parameters])
                for j in range(result.shape[1]):
                    print '{} {}'.format(j, thresholds[j])
                if k==0:
                    mean_corr=[]
                    for z, param in enumerate(parameters):
                        mean_corr.append(np.mean(parameter_corr[k, 3*z:3*z+3, 3*z:3*z+3]))
                    f=open(os.path.join(self.settings.result_folder, '{}_mean_itercorr.pkl'.format(test.func_name)), 'w')
                    pickle.dump(mean_corr, f); f.close()
        #plotting this correlation
        cmap = brewer2mpl.get_map('RdBu', 'diverging', 3).mpl_colormap
        for k, comparison in enumerate(comparisons):
            f, ax = p.subplots(1, sharex=True, sharey=True)
            zou=ax.pcolormesh(parameter_corr[k],
                          cmap=cmap,
                          norm=mpl.colors.Normalize(vmin=0, vmax=1),
                          edgecolors = 'None')
            ax.set_xlim(0,result.shape[1])
            ax.set_ylim(0,result.shape[1])
            ax.set_title("P-values, {}'s correlation for {}".format(test.func_name, comparison))
            f.colorbar(zou)
            if sh:
                p.show()
            else:
                p.savefig(os.path.join(self.settings.result_folder, 'Pvalues{}_CTRL{}.png'.format(comparison, testCtrl)))
                
        return np.array(thresholds)
    
    def _correlationExpParam(self, comparisons, parameters, result, testCtrl, sh,thresholds=None):
        cmap = brewer2mpl.get_map('RdBu', 'diverging', 11).mpl_colormap
        
        if thresholds is not None:
            result=np.array(result)
            min_=0; max_=1
            title = "Hits"  
            for k in range(result.shape[0]):
                for j in range(result.shape[1]):
                    threshold = thresholds[j]
                    result[k,j]=result[k,j]<threshold
        else:
            title = 'P-values'
            min_=10**(-5)
            max_=0.1
        
        axes= p.subplots(len(comparisons), sharex=True)
        k=0
        if len(comparisons)==1:
            #have to do this because if this is the case then axes[1] is not a list
            zou=axes[1].pcolormesh(result[k],
                      cmap=cmap,
                      norm = mpl.colors.Normalize(vmin= min_, vmax=max_),
                      edgecolors = 'None')
            axes[1].set_ylim(0, result[0].shape[0])
            axes[1].set_xlim(0, result[0].shape[1])
            axes[1].set_title("{} for file {}, CTRL {}".format(title, comparisons[k], testCtrl))

        else:
            for k in range(len(comparisons)):
                zou=axes[1][k].pcolormesh(result[k],
                      cmap=cmap,
                      norm = mpl.colors.Normalize(vmin=min_ , vmax=max_),
                      edgecolors = 'None')
                axes[1][k].set_ylim(0, result[0].shape[0])
                axes[1][k].set_xlim(0, result[0].shape[1])
                axes[1][k].set_title("{} for file {}, CTRL {}".format(title, comparisons[k], testCtrl))
        axes[0].colorbar(zou)        
        
        if sh:
            p.show()
        else:
            p.savefig(os.path.join(self.settings.result_folder, self.settings.figure_name.format(title, testCtrl)))
        return result
    
    def plot_heatmaps(self, testCtrl=False, iterations=[0,1,2], sh=None, saveOnly=False, loadOnly=False):
        if sh==None:
            sh=show
        parameters = None
        if not loadOnly:
        #getting list of all parameter sets, from all iterations and taking parameter sets that are in all iterations
            for iter_ in iterations:            
                try:
                    f=open(os.path.join(self.settings.result_folder, self.settings.clustering_filename.format(iter_)), 'r')
                    d=pickle.load(f); f.close()
                except IOError:
                    raise
                if parameters == None:
                    parameters = sorted(d.keys(), key=itemgetter(0, 9,5, 2, 7, 8))
                else:
                    paramCourants = sorted(d.keys(), key=itemgetter(0, 9,5, 2, 7, 8))
                    parameters = filter(lambda x: x in paramCourants, parameters)
        #getting list of all siRNAs for which p-values have been calculated in at least one case
    
            if not testCtrl:
                pvalL = filter(lambda x: 'pval' in x and 'CTRL' not in x and 'png' not in x, os.listdir(self.settings.result_folder))
                siRNAL=[el.split('_')[-1][:-4] for el in pvalL]
            else:
                #in this case we're working with ctrl experiments so siRNAL will contain plate names
                pvalL = filter(lambda x: 'pval' in x and 'CTRL' in x and 'png' not in x, os.listdir(self.settings.result_folder))
                siRNAL=['{}_{}'.format(el.split('_')[-2],el.split('_')[-1][:-4]) for el in pvalL]
            
            siRNAL = sorted(Counter(siRNAL).keys())
            platesL=[]
            comparisons = [el.split('_')[1] for el in pvalL] 
            comparisons = Counter(comparisons).keys()
            print comparisons
    
    
    #idea: get iterations of same parameter set next to one another
            result = np.empty(shape=(len(comparisons), len(parameters),len(iterations)), dtype=object)
            result.fill(None)
            
            for k, comparison in enumerate(comparisons):
                iterations=['iter{}'.format(m) for m in iterations]
                if testCtrl:
                    ctrl_iterations=['_0CTRL', '_1CTRL']#[comparison+'_0CTRL', comparison+'_1CTRL']
                else:
                    ctrl_iterations=['']
                for i, iter_ in enumerate(iterations):
                    for ctrl_iter in ctrl_iterations:
                        for siRNA in siRNAL:
                            print '---', siRNA
                            platesL.append(siRNA)
                            try:
                                f=open(os.path.join(self.settings.result_folder, 'pval_{}_{}{}_{}.pkl'.format(comparison, iter_, ctrl_iter, siRNA)))
                                d=pickle.load(f); f.close()
                            except IOError:
                                print "zou", siRNA, iter_
                                d={}
                            for j,parameter_set in enumerate(parameters):
                                
                                if parameter_set in d:
                                    if type(d[parameter_set])==list and d[parameter_set]!=[]:
                                        try:
                                            result[k,j,i].append(d[parameter_set][0])
                                        except AttributeError:
                                            result[k,j,i]=d[parameter_set]
                                    elif d[parameter_set]==[]:
                                        print siRNA, 'liste vide de pvalues'
                                        result[k,j,i].append(np.nan)
                                    else:
                                        try:
                                            result[k,j,i].append(d[parameter_set])
                                        except AttributeError:
                                            result[k,j,i]=[d[parameter_set]]
                                else:
                                    try:
                                        result[k,j,i].append(np.nan)
                                    except AttributeError:
                                        result[k,j,i] = [np.nan]
        
#        if testCtrl == True:
#            for j in range(result.shape[1]):
#                f=p.figure()
#                ax=f.add_subplot(111)
#                for i in range(result.shape[2]):
#                    nb,bins,patches=ax.hist(result[0,j,i], bins=100, range=(0,0.001), color = couleurs[i],alpha=0.2)
#                f.savefig( 'histParam{}.png'.format(j) )
#
        else:
            f=open(os.path.join(self.settings.result_folder, 'pickledPval_{}.pkl'.format(testCtrl)))
            result, siRNAL, platesL = pickle.load(f); f.close()
                    
        if saveOnly:
            f=open(os.path.join(self.settings.result_folder, 'pickledPval{}.pkl'.format(testCtrl)), 'w')
            pickle.dump((result, siRNAL, platesL), f); f.close()
            return
            
        new=None
        for k in range(result.shape[0]):
            N=None
            for i in range(result.shape[2]):
                for j in range(result.shape[1]):
                    try:
                        N=np.vstack((N, np.array(result[k,j,i]) )) if N is not None else np.array(result[k,j,i])
                    except:
                        pdb.set_trace()
            if new is None:
                new = N[np.newaxis, :]
            else:
                new=np.vstack((new, N[np.newaxis, :]))
        result=new
        print result.shape
        if len(np.where(np.isnan(result))[2])>0:
            pdb.set_trace()
        genesToDel = np.where(np.isnan(result))[2]
        result = np.delete(result, genesToDel, 2)
        print "And after deleting lines with NaN", result.shape           
     
        #calculating correlations between differents parameter sets and computing p-values thresholds for different parameters
        thresholds = self._correlationParamParam(comparisons, parameters, result, testCtrl,iterations, sh)
        
        #plotting results: experiments vs parameters. Just another way to see the correlations between different parameters
        _ = self._correlationExpParam(comparisons,parameters, result, testCtrl, sh)
        
        #plotting results in terms of hits, and returning the matrix with 1 where the hits were found
        hit_matrix = self._correlationExpParam(comparisons, parameters, result, testCtrl, sh, thresholds=thresholds)
        
        #plotting correlations for different iterations and parameters between the hit lists
        _ = self._correlationParamParam(comparisons, parameters, hit_matrix, testCtrl,iterations, sh, test=pearsonr, save=False)
        #giving gene list
        genes=platesL
        if not testCtrl:
            genes=[]
            dictSiEntrez=siEntrez(self.settings.mitocheck_file)
            for siRNA in siRNAL:
                try:
                    genes.append(dictSiEntrez[siRNA])
                except KeyError:
                    pdb.set_trace()
        genes = np.delete(np.array(genes), genesToDel)
        siRNAL = np.delete(np.array(siRNAL), genesToDel)
        return hit_matrix, result, parameters, genes, siRNAL


class clusteringExperiments():
    def __init__(self, settings_file, experimentList,n_cluster, div_name, bins_type, cost_type, bin_size, 
                 init='k-means++',lambda_=10, M=None, dist_weights=None, batch_size=100, init_size=1000,
                 n_init=10, verbose=0, ddim = 0, iter_=0):
        
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
        self.iter_=iter_
        
    def parameters(self, pcaParameter, with_n_cluster=False):
        r={
           'bins_type': self.bins_type,
           'bin_size':self.bin_size,
           'div_name':self.div_name,
           'lambda':self.lambda_,
           'cost_type':self.cost_type,
           'dist_weights':self.dist_weights,
           'ddimensional':self.ddimensional,
           'PCA': pcaParameter[0],
           'whiten': pcaParameter[1]
           }
        if with_n_cluster:
            r.update({'n_cluster':self.n_cluster})
            
        l = sorted(zip(r.keys(), r.values()), key=itemgetter(0))
        return tuple(l)
    
    def _dataPrep(self):
        '''
        This function gets the data from each experiment, and then selects the lines corresponding to the representatives of
        the experiment which have been selected in the previous step (summarizingExperiment)
        '''
        num_features=[None, None, None]; dataL=[]
        hist_features = [{hist_feature : [] for hist_feature in featuresHisto} for k in range(3)]

        count=0; ctrl_count=0; ctrl_exp_count=0
        self.actual_expList = []
        for i, experiment in enumerate(self.expList):
            print '------{}/{} experiments'.format(i, len(self.expList))
            try:
                f=open(os.path.join(self.settings.result_folder, self.settings.summary_filename.format(self.iter_, experiment[0], experiment[1])))
                representatives = pickle.load(f); f.close()
            except IOError:
                sys.stderr.write("Representative trajectories for experiment {} have not been calculated at all.".format(experiment))

            else:
                try:
                    _, r, histNtot,  _,_, _, _, _, _ = histConcatenation(self.settings.data_folder, [experiment], self.settings.mitocheck_file,
                                        self.settings.quality_control_file, verbose=self.verbose)
                except AttributeError:
                    sys.stderr.write("No data for experiment {}".format(experiment))
                else:        
                    for i, pcaParameter in enumerate(self.settings.pcaParameters):
                        curr_parameters = self.parameters(pcaParameter, with_n_cluster=False)
                        try:
                            num_features[i]=r[representatives[curr_parameters], :self.settings.nb_feat_num] if num_features[i] is None else \
                            np.vstack((num_features[i], r[representatives[curr_parameters], :self.settings.nb_feat_num]))
            
                        except KeyError:
                            sys.stderr.write("Representative trajectories for experiment {}, parameters {} have not yet been calculated.".format(experiment, curr_parameters))
                        else:                    
                            
                            for hist_feature in hist_features[i]:
                                hist_features[i][hist_feature].extend([histNtot[hist_feature][k] for k in representatives[curr_parameters]])
                            if i==0:
                                self.actual_expList.append(experiment)
                                count+=len(representatives[curr_parameters]) 
                                if is_ctrl(experiment):
                                    ctrl_exp_count+=1
                                    ctrl_count +=len(representatives[curr_parameters]) 
        if self.verbose:   
            print "{} representatives from {} experiments, among which {} ctrl representatives from {} ctrl experiments".format(count,\
                                         len(self.actual_expList), ctrl_count, ctrl_exp_count)  
            print num_features[0].shape, 'not normalized'
            print 'computing bins, bin type {}, d-dimensional? {}'.format(self.bins_type, bool(self.ddimensional))
        assert(count == num_features[0].shape[0])
        
        for i, pcaParameter in enumerate(self.settings.pcaParameters):            
            if self.ddimensional:
                histogrammeMatrix, bins = computingddBins(hist_features[i], self.mat_hist_sizes[0], bin_type=self.bins_type)
            else:
    
                histogrammeMatrix, bins = computingBins(hist_features[i], self.mat_hist_sizes[0], bin_type=self.bins_type)
            
            if pcaParameter[0]==0:
    #no PCA at all
                self.mean = np.mean(num_features[i],0)
                self.std = np.std(num_features[i],0)
                nnum_features=(num_features[i]-self.mean)/self.std
            
            else:
    #PCA
                nnum_features=_ctrlPca(num_features[i], pcaParameter[1], self.settings.nb_feat_num, self.settings.nb_composantes, self.settings.result_folder)
            dataL.append(np.hstack((nnum_features, histogrammeMatrix)))
        
        return dataL, bins     
    
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
        
    def __call__(self):
        dist_weights = WEIGHTS[self.dist_weights] if not self.ddimensional else ddWEIGHTS[self.dist_weights]
        
        #i. assembling all data 
        dataL, bins=self._dataPrep()
        
        if self.ddimensional:
            self.cost_matrix={(0,0):ddcostMatrix(np.product(self.mat_hist_sizes[0]), self.mat_hist_sizes[0])}
            self.mat_hist_sizes=np.array([[np.product(self.mat_hist_sizes)]])
            #filename = 'dd_'+self.settings.clustering_filename
        else:   
            if self.cost_type=='value':
                print 'cost type value'
                self.cost_matrix=bins
            self.cost_matrix=_costMatrix(self.mat_hist_sizes, self.cost_matrix)
#here we deal with different pca parameters after loading data since it's so long
        for i,pcaParameter in enumerate(self.settings.pcaParameters):
            #ii.clustering the data
            model = histogramMiniKMeans(self.n_cluster,self.lambda_, self.mat_hist_sizes,
                             div_name=self.div_name, M=self.cost_matrix, 
                             dist_weights=dist_weights, nb_feat_num = self.settings.nb_feat_num,
                             init=self.init,
                             batch_size=self.batch_size,
                             init_size =self.init_size, verbose=self.verbose,n_init=self.n_init)
            model.fit(dataL[i])
            
            #iii.saving the centers, mean and std to be further used movie per movie
            self._saveResults(model.cluster_centers_, bins, pcaParameter)        
        return

class summarizingExperiment():
    '''
    The purpose of this class is to take an experiment, find the right number of clusters in its trajectory
    set, and summarize it with N trajectories per cluster.
    '''
    
    def __init__(self,settings_file, experiment,div_name, bins_type, cost_type, bin_size, 
                 init='k-means++',lambda_=10, M=None, dist_weights=None, batch_size=1000, init_size=1000,
                 n_init=5, verbose=0, ddim = 0,iter_=0, random_state=None,max_no_improvement=None):
        
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
        
        self.iter_=iter_
        
    def _dataPrep(self, pcaParameter):
        _, r, hist, _,_, _, _,_, _ = histConcatenation(self.settings.data_folder, [self.experiment], self.settings.mitocheck_file,
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
        if pcaParameter[0]==0:
#no PCA at all
            r=(r-np.mean(r,0))/np.std(r,0)
        else:
#PCA
            r=_ctrlPca(r, pcaParameter[1], self.settings.nb_feat_num, self.settings.nb_composantes, self.settings.result_folder)
        r=np.hstack((r, histogrammeMatrix))
        
        return r, bins
    
    def parameters(self, pcaParameter):
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
           'ddimensional':self.ddimensional,
           'PCA': pcaParameter[0],
           'whiten': pcaParameter[1]
           }
        
        l = sorted(zip(r.keys(), r.values()), key=itemgetter(0))
        return tuple(l)
    
    def __call__(self):
        filename=self.settings.summary_filename.format(self.iter_, self.experiment[0], self.experiment[1])
        fraction = self.settings.fraction
        num_iterations_stability = self.settings.num_iterations_stability
        dist_weights = WEIGHTS[self.dist_weights] if not self.ddimensional else ddWEIGHTS[self.dist_weights]
        self.random_state = check_random_state(self.random_state)
#i.Data preparation: putting histogram data in bins and normalizing. We do that for each of the three possible cases:
#no PCA, PCA no whitening, PCA and whitening.

        for pcaParameter in self.settings.pcaParameters:
            #if self.cost_type==value then returning the bins for further ground metric calculation
            try:
                data, bins=self._dataPrep(pcaParameter)
            except AttributeError:
                return
            
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
                d.update({self.parameters(pcaParameter):representatives})
            else:
                d={self.parameters(pcaParameter):representatives}
                
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
    parser.add_option('--testCtrl', type=str, dest='testCtrl', default=0)
    parser.add_option('--iter', type=int, dest='iter_', default=0)
    
    parser.add_option('-a', type=str, dest='action',default = 'summary')
    
    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)
    parser.add_option('--ddimensional', type=int, dest='ddimensional', default=0)
    parser.add_option("-k", type=int, dest="n_cluster", default=4)
    parser.add_option("-w", type=int, dest="weights", default=0)
    parser.add_option("-l",type=int, dest="lambda_", default=10)
    parser.add_option("--batch_size", dest="batch_size", type=int,default=1000)
    parser.add_option("--init", dest="init", type=str,default='k-means++')
    parser.add_option("--n_init", dest="n_init", type=int,default=10)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    (options, args) = parser.parse_args()
    if getpass.getuser()=='lalil0u':
        settings_file = 'tracking/settings/settings_summaries_Trulove.py'
    else:
        settings_file = 'tracking/settings/settings_summaries_Thalassa.py'
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
                 n_init=options.n_init, verbose=options.verbose, ddim = options.ddimensional, iter_=options.iter_)
        model()
        
    elif options.action=='clustering':
        file_ = open(options.experimentFile, 'r')
        experimentList = pickle.load(file_); file_.close()
        
        model = clusteringExperiments(settings_file, experimentList,options.n_cluster, 
                options.div_name, options.bins_type, options.cost_type, options.bin_size, 
                 init=options.init,lambda_=options.lambda_, dist_weights=options.weights, batch_size=options.batch_size,
                 n_init=options.n_init, verbose=options.verbose, ddim = options.ddimensional, iter_=options.iter_)
        
        model()
    
    elif options.action =='hitFinder':
        model = hitFinder(settings_file, options.siRNA, options.verbose, testCtrl = options.testCtrl, iter_=options.iter_)
        model()
    elif options.action == 'hitCollection':
        model = hitFinder(settings_file, '')
        model.plot_heatmaps(options.testCtrl, iterations=[1,2,3], sh=False, saveOnly=True)
        
    
                 
                 
                 
                 
                 
                 
