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
from random import sample
from scipy.stats import pearsonr
import matplotlib.pyplot as p
import brewer2mpl
from tracking.trajPack import featuresSaved, featuresHisto, featuresNumeriques
from util.listFileManagement import fromShareToCBIO, appendingControl, txtToList
from tracking.trajPack.clustering import histConcatenation,outputBin, usable, correct_from_Nan
from tracking.PyPack.fHacktrack2 import initXml, finirXml
from tracking.histograms.k_means_transportation import DIVERGENCES, _distances
from tracking.histograms.transportation import costMatrix, computingBins
from util.listFileManagement import expSi, strToTuple, siEntrez, EnsemblEntrezTrad, geneListToFile
from util.sandbox import concatCtrl
from tracking.plots import plotAlignedTraj
from sklearn.decomposition import PCA

from util import settings
from scipy.stats.stats import ks_2samp, spearmanr, scoreatpercentile, percentileofscore
from scipy.stats import chi2

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

stats = importr('stats')

nb_exp_list=[100, 500, 1000, 2000]
parameters=[
#(('bin_size', 10),
#  ('bins_type', 'quantile'),
#  ('cost_type', 'number'),
#  ('div_name', 'etransportation'),
#  ('lambda', 10)),
# (('bin_size', 10),
#  ('bins_type', 'quantile'),
#  ('cost_type', 'number'),
#  ('div_name', 'total_variation'),
#  ('lambda', 10)),
 (('bin_size', 10),
  ('bins_type', 'quantile'),
  ('cost_type', 'number'),
  ('div_name', 'KS'),
  ('lambda', 10))]

def plotTrajectories(verbose, expList, labels, data_folder, result_folder):
    resultCour=[]
    clusters = Counter(labels).keys()
    for plate, well in expList:
        if verbose:
            print "Taking care of plate {}, well {}".format(plate, well)
    
        f=open(os.path.join(data_folder,plate, 'hist_tabFeatures_{}.pkl'.format(well)))
        tab, coord, _=pickle.load(f); f.close()
        tab, toDel = correct_from_Nan(tab, perMovie=True)
        new_coord=[]
        for i,el in enumerate(coord):
            if i not in toDel:
                new_coord.append(el)
        resultCour.extend(new_coord)

    XX=[]; YY=[]
    for k in range(len(resultCour)):
        X=np.array(resultCour[k][1]); X-=X[0]; XX.append(X)
        Y=np.array( resultCour[k][2]); Y-=Y[0]; YY.append(Y)
    minx,maxx = min([min(X) for X in XX]), max([max(X) for X in XX])
    miny,maxy = min([min(Y) for Y in YY]), max([max(Y) for Y in YY])
    
    for cluster in clusters:
        for i in np.where(np.array(labels)==cluster)[0][:10]:
            if not os.path.isdir(os.path.join(result_folder, 'images')): 
                os.mkdir(os.path.join(result_folder, 'images'))
            plotAlignedTraj(XX[i], YY[i], minx, maxx, miny, maxy, show=False, val=cluster, 
                 name=os.path.join(result_folder, 'images', 'cluster{}_{}.png'.format(cluster, i)))

    return 1

def plotDistances(folder, filename='all_distances_whole.pkl', ctrl_filename ="all_distances_whole_CTRL.pkl", sigma=0.1, binSize=10,texts=None):
    f=open(os.path.join(folder, filename))
    distances = pickle.load(f);f.close()
    
    f=open(os.path.join(folder, ctrl_filename))
    ctrl_distances = pickle.load(f); f.close()
    
    for param in parameters:
        _,_,_, ctrl_arr=ctrl_distances[param]
        _,expList,_, exp_arr =distances[param]
        if texts is not None:
            assert(len(expList)==len(texts))
        
        print parameters.index(param), ctrl_arr.shape
        pca = PCA(n_components=16)
        moy=np.mean(ctrl_arr,0); std=np.std(ctrl_arr, 0)
        ctrl_arr=pca.fit_transform((ctrl_arr-moy)/std) #six premieres composantes pour garder le gros ~94% 
        exp_arr=pca.transform((exp_arr-moy)/std)[:,:6]
        f=p.figure(); ax=f.add_subplot(111)
        ax.scatter(exp_arr[:,0], exp_arr[:,1], color='red', alpha=0.5)
        if texts is not None:
            for k in range(exp_arr.shape[0]):
                ax.text(exp_arr[k,0], exp_arr[k,1], texts[k])
        ax.scatter(ctrl_arr[:,0], ctrl_arr[:,1], color='green', alpha=0.5)
        p.show()
        
        outputBin(np.vstack((ctrl_arr[:,:2], exp_arr[:,:2])), ctrl_arr.shape[0], 1, [exp_arr.shape[0]], binSize, sigma)
        
#        p.plot(np.cumsum(pca.explained_variance_ratio_), label=param[3][1])
#    p.grid(True)
#    p.legend()
#    p.show()

def hitDistances(folder, filename='all_distances_whole2.pkl', ctrl_filename ="all_distances_whole2_CTRL.pkl", threshold=0.05, sup=False, renorm_first_statistic=False,
                 renorm_second_statistic=True, redo=False, without_mitotic_hits=False):
    '''
    This function collects all distances or p-values from files, then with or without renormalizing them with respect to control distributions
    combines them with Fisher's method and finally with or without renormalizing them computes the hit list
    
    '''
    #i. Collecting distances
        #parameters : names of quality control and mapping files
    
    exp = collectingDistances(filename, folder, testCtrl=False, qc_filename='../data/mapping_2014/qc_export.txt',mapping_filename='../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt',
                              redo=redo)
    ctrl = collectingDistances(ctrl_filename, folder, testCtrl=True, qc_filename='../data/mapping_2014/qc_export.txt',mapping_filename='../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt',
                               redo=redo)        
    r=[]
    
    #ii. Computing renormalized distributions if renorm_first_statistic==True, else simply combining them using Fisher's method
        #and returning p-val anq q-val distributions, normalized or not
        #parameters: renorm_first_statistic, renorm_second_statistic
                    #sup: if True, indicates that outliers have big values (eg distances) else outliers have small values (eg p-values)

    ctrl_pval, ctrl_qval, combined_pval, combined_qval = empiricalDistributions({param:ctrl[param][-1] for param in parameters},
                                                           {param:exp[param][-1] for param in parameters}, folder, sup=sup,
                                                           renorm_first_statistic=renorm_first_statistic, renorm_second_statistic=renorm_second_statistic,
                                                           redo=redo)
    for param in parameters:
        siRNAL, expL, geneL, _ = exp[param]
        platesL,_,_,_=ctrl[param]
        curr_pval=combined_pval[param]; curr_qval=combined_qval[param]
        
        siRNA_count=Counter(siRNAL)
        gene_count=Counter(geneL)
        
        ctrl_hit=[]
        for k in range(len(platesL)):
            if ctrl_qval[param][k]<threshold:
                ctrl_hit.append(platesL[k])
    
        exp_hit=[]
        siRNA_hit=[]
        gene_hit=[]

        
        for k in range(len(siRNAL)):
            if curr_qval[k]<threshold:
                exp_hit.append(expL[k])
                siRNA_hit.append(siRNAL[k])
                gene_hit.append(geneL[k])
                
        siRNA_hit=Counter(siRNA_hit)
        siRNA_hit_perc={siRNA:siRNA_hit[siRNA]/float(siRNA_count[siRNA]) for siRNA in siRNA_hit}
        
        siRNA_highconf = filter(lambda x: siRNA_hit_perc[x]>0.5 and siRNA_hit[x]>=2, siRNA_hit_perc)
        
        if without_mitotic_hits:
            print "Filtering out mitotic hits from Mitocheck supp table 2"
            mito_hits = txtToList('../data/mitocheck_hits.txt')[:,0]
            siRNA_highconf=[siRNA for siRNA in siRNA_highconf if siRNA not in mito_hits]
        
        gene_highconf = Counter([geneL[siRNAL.index(siRNA)] for siRNA in siRNA_highconf])
        #gene_highconf={gene:gene_highconf[gene]/float(gene_count[gene]) for gene in gene_highconf}
        gene_highconf=filter(lambda x: gene_highconf[x]>=1, gene_highconf)
            
        trad = EnsemblEntrezTrad('../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt')
        gene_hits_Ensembl=[trad[el] for el in gene_hit]
        gene_highconf_Ensembl=[trad[el] for el in gene_highconf]
        gene_Ensembl = [trad[el] for el in gene_count]
        
        print 'Ctrl hits ', len(ctrl_hit), 'out of', len(platesL), 'ie ', len(ctrl_hit)/float(len(platesL))
        print 'Experiences ', len(exp_hit), 'out of',len(expL), 'ie ', len(exp_hit)/float(len(expL))
        print 'siRNA high conf ', len(siRNA_highconf), 'out of',len(siRNA_count), 'ie', len(siRNA_highconf)/float(len(siRNA_count))
        print 'Genes high conf', len(gene_highconf), 'out of', len(gene_count), 'ie ', len(gene_highconf)/float(len(gene_count))

        geneListToFile(gene_hits_Ensembl, 'gene_hits_p{}.txt'.format(parameters.index(param)))
        geneListToFile(gene_highconf_Ensembl, 'gene_high_conf_p{}.txt'.format(parameters.index(param)))
        geneListToFile(gene_Ensembl, 'gene_list_p{}.txt'.format(parameters.index(param)))
        
        r.append([exp_hit, gene_hit, gene_highconf])
        
        result = defaultdict(list)
        for i,gene in enumerate(geneL):
            if gene not in result: 
                result[gene]=defaultdict(list)
            result[gene][siRNAL[i]].append(curr_qval[i])

    return r, siRNA_highconf, siRNAL, curr_qval

    
def collectingDistances(filename, folder, qc_filename='../data/mapping_2014/qc_export.txt',mapping_filename='../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt', testCtrl =False,
                        redo=False):
    if filename not in os.listdir(folder) or redo:
        if not testCtrl:
            files = filter(lambda x: 'distances_whole_5Ctrl' in x and 'CTRL' not in x and 'all' not in x, os.listdir(folder))
            yqualDict=expSi(qc_filename, sens=0)
            dictSiEntrez=siEntrez(mapping_filename)
        else:
            files = filter(lambda x: 'distances_whole_5Ctrl' in x and 'CTRL' in x  and 'all' not in x, os.listdir(folder))
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
                    siRNA = file_.split('_')[-1][:-4]
                    l= Counter(np.where(~np.isnan(d[param]))[0]).keys()
                    
#                    if len(l)==1 and not testCtrl:
#                        #meaning it is an orphan siRNA with only one experiment. But do we really want to take them out? The experiment is real after all
#                        continue 
#                    
                    result[param][0].extend([siRNA for k in range(len(l))])
                    if not testCtrl:
                        try:
                            gene = dictSiEntrez[siRNA]
                        except:
                            pdb.set_trace()
                        expList = np.array(strToTuple(yqualDict[siRNA], os.listdir('/share/data20T/mitocheck/tracking_results')))
                       
                        if d[param].shape[0]==len(expList):
                            used_experiments = expList
                        else:
                            pdb.set_trace()
                            used_experiments = expList[np.where(usable('/share/data20T/mitocheck/tracking_results', expList,
                                                                       qc='../data/qc_export.txt',mitocheck='../data/mitocheck_siRNAs_target_genes_Ens75.txt'))]
                            
                        result[param][1].extend(list(used_experiments[l]))
                        result[param][2].extend([gene for k in range(len(l))])

                    
                    result[param][3]=np.vstack((result[param][3], d[param][l])) if result[param][3] is not None else d[param][l]

        f=open(os.path.join(folder, filename), 'w')
        pickle.dump(result, f); f.close()

    else:
        f=open(os.path.join(folder, filename))
        result=pickle.load(f); f.close()
    
    return result

    

def empiricalPvalues(dist_controls, dist_exp, folder, name, sup=False):
    empirical_pval = {param : np.zeros(shape = (dist_exp[param].shape[0], dist_exp[param].shape[1]), dtype=float) for param in parameters}
    empirical_qval = {param : np.zeros(shape = (dist_exp[param].shape[0], dist_exp[param].shape[1]), dtype=float) for param in parameters}

    for param in parameters:
        #for all parameter sets
        for k in range(dist_exp[param].shape[1]):
            #for all features
            print '---------------',k
            for j in range(dist_exp[param].shape[0]):
                #for all experiments
                
                if sup:
    #if we are dealing with distances to calculate empirical distributions,
    #we're looking for distances that are big, so a small p-value should be for big distances.
    #Hence our empirical p-values are 100 - the percentage of distances that are lower
                    empirical_pval[param][j,k] = (100 - percentileofscore(dist_controls[param][:,k], dist_exp[param][j,k]))/100.0
                else:
                    empirical_pval[param][j,k] = percentileofscore(dist_controls[param][:,k], dist_exp[param][j,k])/100.0
                    
    #The minimum precision is 1/#controls. It is not possible that there is p-value 0 here        
            empirical_pval[param][np.where(empirical_pval[param]==0)]=1/float(dist_controls[param].shape[0])
    #Calling R function p-adjust, applying Benjamini-Hochberg method
            empirical_qval[param] = stats.p_adjust(FloatVector(list(empirical_pval[param][:,k])), method = 'BH')
        #plotting empirical p-value distribution
        f=p.figure()
        ax=f.add_subplot(111)
        nb,bins,patches=ax.hist(empirical_pval[param], bins=50, alpha=0.2)
        nb,bins,patches=ax.hist(empirical_qval[param], bins=50, alpha=0.2, color='red')
        ax.set_xlim(0,1)
        p.savefig(os.path.join(folder, 'pvalParam{}_{}.png'.format(parameters.index(param), name)))
        p.close('all')
            
    return empirical_pval, empirical_qval

def empiricalDistributions(dist_controls, dist_exp, folder, sup=False, union=False, redo=False, renorm_first_statistic=False, renorm_second_statistic=True):
    '''
    This function computes empirical p-value distributions
    dist_controls et dist_exp are dictionaries with keys=parameter sets
    dist_controls[param] array of size nb experiments x nb of features
    '''
    param=parameters[0]
    print "First, univariate p-values"
    if renorm_first_statistic:
        if 'empirical_p_qval.pkl' not in os.listdir(folder) or redo:
            ctrl_pval, ctrl_qval = empiricalPvalues(dist_controls, dist_controls, folder,name='ctrlPval', sup=sup)
            empirical_pval, empirical_qval = empiricalPvalues(dist_controls, dist_exp, folder,name='expPval', sup=sup)
            
            
            f = open(os.path.join(folder, 'empirical_p_qval.pkl'), 'w')
            pickle.dump((ctrl_pval, ctrl_qval, empirical_pval, empirical_qval),f); f.close()
        else:
            f=open(os.path.join(folder, 'empirical_p_qval.pkl'))
            ctrl_pval, ctrl_qval, empirical_pval, empirical_qval = pickle.load(f); f.close()
    else:
        ctrl_pval = dist_controls; empirical_pval = dist_exp
    
    print 'Second, combining p-values'
    # statistical test according to Fisher's method http://en.wikipedia.org/wiki/Fisher%27s_method
    ctrl_combined_pval = {param : np.zeros(shape = (ctrl_pval[param].shape[0],1), dtype=float) for param in parameters}
    combined_pval = {param : np.zeros(shape = (empirical_pval[param].shape[0],1), dtype=float) for param in parameters}

    ctrl_combined_qval = {param : np.zeros(shape = (ctrl_pval[param].shape[0],), dtype=float) for param in parameters}
    combined_qval = {param : np.zeros(shape = (empirical_pval[param].shape[0],), dtype=float) for param in parameters}
    for param in parameters:
        stat = -2*np.sum(np.log(ctrl_pval[param]),1) 
        ctrl_combined_pval[param] = stat[:,np.newaxis]
        
        stat = -2*np.sum(np.log(empirical_pval[param]),1)
        combined_pval[param] = stat[:,np.newaxis]
    
    
    if renorm_second_statistic:
        if 'comb_empirical_p_qval.pkl' not in os.listdir(folder) or redo:
            ctrl_combined_pval2, ctrl_combined_qval = empiricalPvalues(ctrl_combined_pval, ctrl_combined_pval, folder,name='ctrlCombinedStat', sup=True)
            combined_pval, combined_qval = empiricalPvalues(ctrl_combined_pval, combined_pval, folder,name='expCombinedStat', sup=True)
            
            f = open(os.path.join(folder, 'comb_empirical_p_qval.pkl'), 'w')
            pickle.dump((ctrl_combined_pval2, ctrl_combined_qval, combined_pval, combined_qval),f); f.close()
        else:
            f = open(os.path.join(folder, 'comb_empirical_p_qval.pkl'), 'r')
            ctrl_combined_pval2, ctrl_combined_qval, combined_pval, combined_qval=pickle.load(f); f.close()

    else:
        for param in parameters:
            ctrl_combined_pval[param] = chi2.sf(ctrl_combined_pval[param], 2*ctrl_pval[param].shape[1])
            ctrl_combined_qval[param]= ctrl_combined_pval[param]*ctrl_combined_pval[param].shape[0]/(1+np.argsort(ctrl_combined_pval[param]))
            
            combined_pval[param] = chi2.sf(combined_pval[param], 2*empirical_pval[param].shape[1])
            combined_qval[param]= combined_pval[param]*combined_pval[param].shape[0]/(1+np.argsort(combined_pval[param]))
    print 'End'
    return ctrl_combined_pval2, ctrl_combined_qval, combined_pval, combined_qval


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
        targetExp={}; yeSiExp=expSi('../data/mapping_2014/qc_export.txt', 0)
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
        
        self.currInterestFeatures = list(featuresNumeriques)
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
        
        _,r, _, self.expList,_, length, _, _, _ = histConcatenation(self.settings.data_folder, self.expList, self.settings.mitocheck_file,
                                        self.settings.quality_control_file, verbose=self.verbose, perMovie=True)
                    
        for i in range(len(length)):
            for k,feature in enumerate(self.currInterestFeatures):
                histDict[feature].append(r[np.sum(length[:i]):np.sum(length[:i+1]),featuresSaved.index(feature)])
        
        f=open(os.path.join(self.settings.result_folder, 'distExp_ctrl_{}_{}.pkl'.format(self.bins_type, self.bin_size)))
        bins = pickle.load(f); f.close()
        
        if self.div_name !='KS':
            raise TypeError
#            histogrammes, bins = computingBins(histDict, [self.bin_size for k in range(len(self.currInterestFeatures))], self.bins_type, previous_binning=bins)
#            return histogrammes, bins
        else:
            return histDict, None
    
    def ctrlHistograms(self, bins, toDel=[]):
        '''
        Since we want to calculate the length of an experiment to its controls, we only need to keep one table per
        plate, and not one table per control experiment
        '''
        plates = Counter(np.array(self.expList)[:,0]).keys()
        if self.plate is not None:
            assert len(plates)==1
            assert plates[0]==self.plate
            
        histDict = defaultdict(list)
        length=[]
        for plate in plates:
            
            #meaning we are dealing with control control tests
            if len(toDel)>0:
                ctrlExpList = appendingControl([plate])
                if self.verbose:
                    print 'before cleaning', len(ctrlExpList)
                
                true_ctrl = filter(lambda x: x not in toDel, range(len(ctrlExpList)))
                ctrlExpList = list(np.array(ctrlExpList)[true_ctrl])
                if self.verbose:
                    print 'after cleaning', len(ctrlExpList)
            #saving the controls that were used to compare the 'false experiments' to use them to compare the real experiments of the same plate
             
                f=open(os.path.join(self.settings.result_folder, self.settings.ctrl_exp_filename.format(plate)), 'w')
                pickle.dump(ctrlExpList, f); f.close()
                
            else:
            #loading the controls that were used to compute control control p-values
                try:
                    f=open(os.path.join(self.settings.result_folder, self.settings.ctrl_exp_filename.format(plate)))
                    ctrlExpList = pickle.load(f); f.close()
                except IOError:
                    sys.stderr.write('No file registering used ctrls for plate {}'.format(plate))
                    ctrlExpList=None
                    
            try:
                _,curr_r, _, _,_, curr_length, _, _, _ = histConcatenation(self.settings.data_folder, ctrlExpList, self.settings.mitocheck_file,
                                            self.settings.quality_control_file, verbose=self.verbose)
            except:
                print "Problem with controls from plate {}".format(plate)
                length.append(np.nan)
                for k, feature in enumerate(self.currInterestFeatures):
                    histDict[feature].append(np.nan)
            else:    
                length.append(np.sum(curr_length))
                assert(np.sum(curr_length)==curr_r.shape[0])                        
                for k, feature in enumerate(self.currInterestFeatures):
                    histDict[feature].append(curr_r[:,featuresSaved.index(feature)])
        assert(len(length)==len(plates))
        assert(len(histDict[feature])==len(plates))
        if self.div_name !='KS':
            raise TypeError
#            histogrammes, _ = computingBins(histDict, [self.bin_size for k in range(len(self.currInterestFeatures))], self.bins_type, previous_binning=bins)
#            return plates, histogrammes
        else:
            return plates, histDict
    
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
            
            if "transportation" in self.div_name:
                if self.cost_type!='number':
                    print 'Need to set up cost matrix to match with cost_type value'
                    raise
            for k, feature in enumerate(self.currInterestFeatures):
                if self.div_name !='KS':
                    raise TypeError
#                    ctrl_hist = ctrl_histogrammes[corresponding_ctrl]
#                    dist_weights = np.zeros(shape=(len(self.currInterestFeatures)+1,))
#                    dist_weights[k+1]=1
#                    mat_hist_sizes = np.array([[self.bin_size for j in range(len(self.currInterestFeatures))]])
#                    distances[i, k]=_distances(histogrammes[np.newaxis, i], ctrl_hist[np.newaxis,:],div_name=self.div_name, lambda_=self.lambda_, M=None,
#                                         mat_hist_sizes=mat_hist_sizes, nb_feat_num=0, dist_weights=dist_weights)[0][0]
                else:
                    h1 = histogrammes[feature][i]
                    h2 = ctrl_histogrammes[feature][corresponding_ctrl]
                    if not (type(h2)==float and np.isnan(h2)):
                #taking the p-value because the distance does not take into account the sample sizes explicitly
                        distances[i,k]=ks_2samp(h1[~np.isnan(h1)], h2[~np.isnan(h2)])[1]
                    else:
                        distances[i,k]=np.nan
        return distances
    
    def alreadyDone(self):
        return self.settings.outputFile.format(self.siRNA) in os.listdir(self.settings.result_folder)
    
    def saveResults(self, distances):
#        if self.settings.outputFile.format(self.siRNA) in os.listdir(self.settings.result_folder):
#            f=open(os.path.join(self.settings.result_folder, self.settings.outputFile.format(self.siRNA)))
#            d=pickle.load(f)
#            f.close()
#            try:
#                l=d[self.parameters()]
#            except KeyError:
#                pass
#        else:
#            d={}

        d={self.parameters():distances}
        f=open(os.path.join(self.settings.result_folder, self.settings.outputFile.format(self.siRNA)), 'w')
        pickle.dump(d, f)
        f.close()
        return
    
    def __call__(self):
        
    #before anything, testing for existence if not redo anyway
        if not self.settings.redo and self.alreadyDone():
            print 'Already done'
            return
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
            ctrl = np.array(appendingControl([self.plate]))
            usable_ctrl = ctrl[np.where( usable(self.settings.data_folder, ctrl, qc=self.settings.quality_control_file, mitocheck=self.settings.mitocheck_file))]
            
            if len(usable_ctrl)<2:
                return
            elif len(usable_ctrl)<6:
                #randomly selecting two wells of the plate that will be used to be compared to the others
                false_exp = np.random.randint(len(usable_ctrl), size=1)
                chosen_ctrl = filter(lambda x: x not in false_exp, range(len(usable_ctrl)))
            else:
                chosen_ctrl = sample(xrange(len(usable_ctrl)), 5)
                false_exp = filter(lambda x: x not in chosen_ctrl, range(len(usable_ctrl)))
            
            self.expList = list(usable_ctrl[false_exp])
            histogrammes, bins = self.getData(self.settings.histDataAsWell)
            
            plates, ctrl_histogrammes = self.ctrlHistograms(bins, toDel=false_exp)
            
            if self.verbose:
                print "all ctrl", ctrl
                print 'all usable ctrl', usable_ctrl
                print 'false exp', self.expList
                
            #ii.calculate the histograms and binnings of experiments, using pre-computed binnings, ON all control experiments 
            #for the plate
    
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
