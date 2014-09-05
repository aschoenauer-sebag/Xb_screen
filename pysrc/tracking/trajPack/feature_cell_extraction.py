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
from util.listFileManagement import expSi, strToTuple, siEntrez, EnsemblEntrezTrad, geneListToFile
from util.sandbox import concatCtrl
from tracking.plots import plotAlignedTraj

from util import settings
from scipy.stats.stats import ks_2samp, spearmanr, scoreatpercentile, percentileofscore
from scipy.stats import chi2

nb_exp_list=[100, 500, 1000, 2000]
parameters=[(('bin_size', 10),
  ('bins_type', 'quantile'),
  ('cost_type', 'number'),
  ('div_name', 'etransportation'),
  ('lambda', 10)),
 (('bin_size', 10),
  ('bins_type', 'quantile'),
  ('cost_type', 'number'),
  ('div_name', 'total_variation'),
  ('lambda', 10))]

def hitDistances(folder, filename='all_distances2.pkl', ctrl_filename ="all_distances2_CTRL.pkl", threshold=0.05):
    if filename not in os.listdir(folder):
        exp = collectingDistances(folder, testCtrl=False)
        ctrl = collectingDistances(folder, testCtrl=True)
        f=open(os.path.join(folder, filename), 'w')
        pickle.dump(exp, f); f.close()
        
        f=open(os.path.join(folder, ctrl_filename), 'w')
        pickle.dump(ctrl,f); f.close()
        
    else:
        f=open(os.path.join(folder, filename))
        exp=pickle.load(f); f.close()
        f=open(os.path.join(folder, ctrl_filename), 'r')
        ctrl =pickle.load(f); f.close()
    r=[]
    combined_pval, combined_qval = empiricalDistributions({param:ctrl[param][-1] for param in parameters},
                                                           {param:exp[param][-1] for param in parameters}, folder, sup=True)
    
    for param in exp:
        siRNAL, expL, geneL, _ = exp[param]
        curr_pval=combined_pval[param]; curr_qval=combined_qval[param]
        
        siRNA_count=Counter(siRNAL)
        gene_count=Counter(geneL)
    
        exp_hit=[]
        siRNA_hit=[]
        gene_hit=[]
        for k in range(len(siRNAL)):
            if curr_pval[k]<threshold and curr_qval[k]<threshold:
                exp_hit.append(expL[k])
                siRNA_hit.append(siRNAL[k])
                gene_hit.append(geneL[k])
                
        siRNA_hit=Counter(siRNA_hit)
        siRNA_hit={siRNA:siRNA_hit[siRNA]/float(siRNA_count[siRNA]) for siRNA in siRNA_hit}
        
        siRNA_highconf = filter(lambda x: siRNA_hit[x]>0.5, siRNA_hit)
        gene_highconf = Counter([geneL[siRNAL.index(siRNA)] for siRNA in siRNA_highconf])
        
        trad = EnsemblEntrezTrad('../data/mapping/mitocheck_siRNAs_target_genes_Ens72.txt')
        gene_hits_Ensembl=[trad[el] for el in gene_hit]
        gene_highconf_Ensembl=[trad[el] for el in gene_highconf]
        gene_Ensembl = [trad[el] for el in gene_count]
        
        print 'Hits ', len(gene_hit)
        print 'High conf', len(gene_highconf)
        print 'Background', len(gene_count)
        geneListToFile(gene_hits_Ensembl, 'gene_hits_p{}.txt'.format(parameters.index(param)))
        geneListToFile(gene_highconf_Ensembl, 'gene_high_conf_p{}.txt'.format(parameters.index(param)))
        geneListToFile(gene_Ensembl, 'gene_list_p{}.txt'.format(parameters.index(param)))
        
        r.append([exp_hit, gene_hit, gene_highconf])
    return r

    
def collectingDistances(folder, testCtrl =False):
    if not testCtrl:
        files = filter(lambda x: 'distances2' in x and 'CTRL' not in x, os.listdir(folder))
        yqualDict=expSi('../data/mapping/qc_export.txt', sens=0)
        dictSiEntrez=siEntrez('../data/mapping/mitocheck_siRNAs_target_genes_Ens72.txt')
    else:
        files = filter(lambda x: 'distances2_CTRL' in x, os.listdir(folder))
    print len(files)

    result={param:[[], [], [], None] for param in parameters}
    
    for file_ in files:
        f=open(os.path.join(folder, file_))
        d=pickle.load(f)
        f.close()
        for param in parameters:
            if param not in d:
                print "NO parameter sets ", file_, 'param ', parameters.index(param)
            else:
                siRNA = file_.split('_')[1][:-4]
                result[param][0].extend([siRNA for k in range(d[param].shape[0])])
                if not testCtrl:
                    try:
                        gene = dictSiEntrez[siRNA]
                    except:
                        pdb.set_trace()
                    result[param][1].extend(yqualDict[siRNA])
                    result[param][2].extend([gene for k in range(d[param].shape[0])])
                
                result[param][3]=np.vstack((result[param][3], d[param])) if result[param][3] is not None else d[param]
            
    return result

def empiricalDistributions(dist_controls, dist_exp, folder, sup=False):
    '''
    This function computes empirical p-value distributions
    dist_controls et dist_exp are dictionaries with keys=parameter sets
    dist_controls[param] array of size nb experiments x nb of features
    '''

    empirical_pval = {param : np.zeros(shape = (dist_exp[param].shape[0], dist_exp[param].shape[1]), dtype=float) for param in parameters}
    empirical_qval = {param : np.zeros(shape = (dist_exp[param].shape[0], dist_exp[param].shape[1]), dtype=float) for param in parameters}

    if 'empirical_p_qval.pkl' not in os.listdir(folder):
        for param in parameters:
            #for all parameter sets
            for k in range(dist_exp[param].shape[1]):
                #for all features
                for j in range(dist_exp[param].shape[0]):
                    #for all experiments
                    if sup:
        #if we are dealing with distances to calculate empirical distributions,
        #we're looking for distances that are big, so a small p-value should be for big distances.
        #Hence our empirical p-values are 100 - the percentage of distances that are lower
                        empirical_pval[param][j,k] = (100 - percentileofscore(dist_controls[param][:,k], dist_exp[param][j,k]))/100.0
                    else:
                        empirical_pval[param][j,k] = percentileofscore(dist_controls[param][:,k], dist_exp[param][j,k])/100.0
            
                empirical_qval[param][:,k] = empirical_pval[param][:,k]*empirical_pval[param].shape[0]/(np.argsort(empirical_pval[param][:,k])+1)
    
            #plotting empirical p-value distribution
            f=p.figure()
            ax=f.add_subplot(111)
            nb,bins,patches=ax.hist(empirical_pval[param], bins=50, alpha=0.2)
            p.savefig(os.path.join(folder, 'pvalParam{}.png'.format(parameters.index(param))))
            p.close('all')
            
        f = open(os.path.join(folder, 'empirical_p_qval.pkl'), 'w')
        pickle.dump((empirical_pval, empirical_qval),f); f.close()
    else:
        f=open(os.path.join(folder, 'empirical_p_qval.pkl'))
        empirical_pval, _ = pickle.load(f); f.close()
        
    
    # statistical test according to Fisher's method http://en.wikipedia.org/wiki/Fisher%27s_method
    combined_pval = {param : np.zeros(shape = (empirical_pval[param].shape[0],), dtype=float) for param in parameters}
    combined_qval = {param : np.zeros(shape = (empirical_pval[param].shape[0],), dtype=float) for param in parameters}
    for param in parameters:
        stat = -2*np.sum(np.log(empirical_pval[param]),1)
        combined_pval[param] = chi2.sf(stat, 2*empirical_pval[param].shape[1])
        combined_qval[param]= combined_pval[param]*combined_pval[param].shape[0]/(1+np.argsort(combined_pval[param]))
    
    return combined_pval, combined_qval


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
                    
        return histogrammes, bins
    
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
        
        return plates, histogrammes
    
    def parameters(self):
        r={
           'bins_type': self.bins_type,
           'bin_size':self.bin_size,
           'div_name':self.div_name,
           'lambda':self.lambda_,
           'cost_type':self.cost_type
           }
        return tuple(sorted(zip(r.keys(), r.values()), key=itemgetter(0)))
    
    def calculateDistances(self, plates, histogrammes, ctrl_histogrammes):
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
                pass
        else:
            d={}

        d.update({self.parameters():distances})
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
            histogrammes, bins = self.getData(self.settings.histDataAsWell)
        
            #iii. calculate the histograms of corresponding controls. Should be per plate then, a dictionary
            plates, ctrl_histogrammes = self.ctrlHistograms(bins)
            
        else:
            
            self.expList = appendingControl([self.plate])
            if self.verbose:
                print "all ctrl", self.expList
                
            #ii.calculate the histograms and binnings of experiments, using pre-computed binnings, ON all control experiments 
            #for the plate
            histogrammes, bins = self.getData(self.settings.histDataAsWell)
            
            plates = [self.plate]
            #randomly selecting two wells of the plate that will be used to be compared to the others
            false_exp = np.random.randint(histogrammes.shape[0], size=2)
            while false_exp[0]==false_exp[1]:
                false_exp = np.random.randint(histogrammes.shape[0], size=2)
            assert(histogrammes.shape[1]==self.bin_size*len(self.currInterestFeatures))
            true_ctrl = filter(lambda x: x not in false_exp, range(histogrammes.shape[0]))
            
            ctrl_histogrammes = histogrammes[true_ctrl]
            histogrammes = histogrammes[false_exp]
            self.expList = list(np.array(self.expList)[false_exp])
            
            if self.verbose:
                print self.expList, true_ctrl, histogrammes.shape[0], ctrl_histogrammes.shape[0]
        
        #iv. calculate the distance from each experiment to its control
        distances = self.calculateDistances(plates, histogrammes, ctrl_histogrammes)
        print distances
        #v. save results
        self.saveResults(distances)
        
        return
    
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
