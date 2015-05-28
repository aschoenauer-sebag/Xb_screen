import getpass, os,pdb, operator, sys, datetime
import numpy as np
import cPickle as pickle
from analyzer.plateSetting import fromXBToWells

if getpass.getuser()=='aschoenauer':
    import matplotlib as mpl
    mpl.use('Agg')
    show = False
elif getpass.getuser()=='lalil0u':
    import matplotlib as mpl
    show=True
    
from sklearn.decomposition import PCA
from optparse import OptionParser
from operator import itemgetter
from itertools import combinations, product
from collections import Counter, defaultdict
from random import sample
from scipy.stats import pearsonr
import matplotlib.pyplot as p
import brewer2mpl
from tracking.trajPack import featuresSaved, featuresHisto, featuresNumeriques,\
    time_windows
from util.listFileManagement import fromShareToCBIO, appendingControl, txtToList
from tracking.trajPack.clustering import histConcatenation,outputBin, correct_from_Nan
from tracking.PyPack.fHacktrack2 import initXml, finirXml
#from tracking.histograms.k_means_transportation import DIVERGENCES, _distances
#from tracking.histograms.transportation import costMatrix, computingBins
from util.listFileManagement import expSi, strToTuple, siEntrez, EnsemblEntrezTrad, geneListToFile, usable_MITO
from util.sandbox import concatCtrl
from util.plots import couleurs, markers
from analyzer import CONTROLS, compoundL, xbL, plates
from analyzer.xb_analysis import xbConcatenation
from analyzer.quality_control import usable_XBSC, computingToRedo

from util import settings
from scipy.stats.stats import ks_2samp, spearmanr, scoreatpercentile, percentileofscore,\
    mannwhitneyu
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
#  ('div_name', 'MW'),
#  ('lambda', 10)),
  (('div_name', 'KS'),
  ('iter', z)) for z in range(5)]

def plotComparison(expDict, inDir,outputFile ="{}_length_distribution.png", filename="cell_cycle_info_{}.pkl", b=50):
    controls=defaultdict(list)

    for gene in expDict:
        f,axes=p.subplots(1, len(expDict[gene]),sharey=True)
        
        for i,exp in enumerate(expDict[gene]):
            folder = filter(lambda x: x[:9] == exp[:9], os.listdir(inDir))[0]
            if folder not in controls:
                for well in ["00074_01", "00315_01"]:
                    try:
                        f=open(os.path.join(inDir, folder, filename.format(well)))
                        d=pickle.load(f)
                        f.close()
                    except:
                        print "Pas ", os.path.join(inDir, folder, filename.format(well))
                        continue
                    else:
                        controls[folder].extend(d['length'].values())
                        
            try:
                f=open(os.path.join(inDir, folder, filename.format(exp[11:]+'_01')))
                d=pickle.load(f)
                f.close()
            except:
                print "Pas ", os.path.join(inDir, folder, filename.format(exp[11:]+'_01'))
                continue
            else:
                axes[i].hist(d['length'].values(), bins=b, color='red', normed=True, alpha=0.5, label=exp)
                axes[i].hist(controls[folder], bins=b, color='green', normed=True, alpha=0.5, label='Ctrl')
        p.legend()
        p.savefig(os.path.join('../resultData/cell_cycle/movies_median', outputFile.format(gene)))
        p.close('all')
        
    return

#def plotDistances(folder, filename='all_distances_whole.pkl', ctrl_filename ="all_distances_whole_CTRL.pkl", sigma=0.1, binSize=10,texts=None):
#    f=open(os.path.join(folder, filename))
#    distances = pickle.load(f);f.close()
#    
#    f=open(os.path.join(folder, ctrl_filename))
#    ctrl_distances = pickle.load(f); f.close()
#    
#    for param in parameters:
#        _,_,_, ctrl_arr=ctrl_distances[param]
#        _,expList,_, exp_arr =distances[param]
#        if texts is not None:
#            assert(len(expList)==len(texts))
#        
#        print parameters.index(param), ctrl_arr.shape
#        pca = PCA(n_components=16)
#        moy=np.mean(ctrl_arr,0); std=np.std(ctrl_arr, 0)
#        ctrl_arr=pca.fit_transform((ctrl_arr-moy)/std) #six premieres composantes pour garder le gros ~94% 
#        exp_arr=pca.transform((exp_arr-moy)/std)[:,:6]
#        f=p.figure(); ax=f.add_subplot(111)
#        ax.scatter(exp_arr[:,0], exp_arr[:,1], color='red', alpha=0.5)
#        if texts is not None:
#            for k in range(exp_arr.shape[0]):
#                ax.text(exp_arr[k,0], exp_arr[k,1], texts[k])
#        ax.scatter(ctrl_arr[:,0], ctrl_arr[:,1], color='green', alpha=0.5)
#        p.show()
#        
#        outputBin(np.vstack((ctrl_arr[:,:2], exp_arr[:,:2])), ctrl_arr.shape[0], 1, [exp_arr.shape[0]], binSize, sigma)
#        
##        p.plot(np.cumsum(pca.explained_variance_ratio_), label=param[3][1])
##    p.grid(True)
##    p.legend()
##    p.show()

def multipleHitDistances(folder, key_name,
                         threshold=0.05, 
                         qc_filename='../data/mapping_2014/qc_export.txt',
                         mapping_filename = '../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt',
                         filename='all_distances_whole_5Ctrl', 
                         combination='min', redo=False,
                         trad=True,
                         without_mitotic_hits=False,
                         without_mean_persistence=True,
                         save=False):
    features = list(featuresNumeriques)
    features.append('mean persistence')
    features.append('mean straight')
    
    expL=None
    if 'comb_empirical_p_qval{}{}.pkl'.format(combination,without_mean_persistence) not in os.listdir(folder) or redo:
        if 'all_distances_{}.pkl'.format(without_mean_persistence) not in os.listdir(folder) or redo:
            expL, geneL, siRNAL, global_result = hitDistances(folder,
                        key_name = key_name, 
                        filename='{}.pkl'.format(filename), 
                        qc_filename=qc_filename,
                        mapping_filename=mapping_filename,
                        redo=redo,
                         without_mean_persistence=True, features=features)
                
    #            curr_expL=curr_result[:,0]
    #            curr_pval=np.array(curr_result[:,-1],dtype=np.float_)
    #            result.append((np.array(curr_expL), np.array(curr_pval)))
    #            if expL is None:
    #                expL=sorted(curr_expL) 
    #                siRNAL=np.array(curr_result[:,2])
    #                geneL=np.array(curr_result[:,1])
    #                
    #            elif expL !=sorted(curr_expL):
    #                print '######################## different experiment list - need to recompute siRNA list'
    #                pdb.set_trace()
    #                expL=filter(lambda x: x in expL, curr_expL)
    #            else:
    #                print 'Same experiment list ', k
                    
            expL=np.array(expL)
            print 'combination of iterations done ', len(expL)           
            print 'Everything together'        
                    
            f=open(os.path.join(folder, 'all_distances_{}.pkl'.format(without_mean_persistence)), 'w')
            pickle.dump([expL, siRNAL, geneL, global_result], f);f.close()
        else:
            f=open(os.path.join(folder, 'all_distances_{}.pkl'.format(without_mean_persistence)))
            expL, siRNAL, geneL, global_result=pickle.load(f); f.close()


        if combination=='min':
            global_pval=np.min(global_result,1)
        elif combination=="median":
            global_pval=np.median(global_result,1)
        elif combination=="max":
            global_pval=np.max(global_result,1)
        elif combination=='mean':
            global_pval=np.mean(global_result,1)
        
        ctrl = collectingDistances("{}_CTRL.pkl".format(filename), folder,key_name =key_name,
                                   testCtrl=True, qc_filename=qc_filename,mapping_filename=mapping_filename,
                                   redo=redo)
        param = (('div_name', 'KS'),('iter', 0))
        ctrl_pval=ctrl[param][-1]
        if without_mean_persistence:
            ctrl_pval=np.delete(ctrl_pval, features.index('mean persistence'),1)
            #ctrl_pval=np.delete(ctrl_pval, features.index('diffusion coefficient'),1)
    
    
        stat = -2*np.sum(np.log(ctrl_pval),1) 
        ctrl_combined_pval = stat[:,np.newaxis]

        ctrl_qval = empiricalPvalues(ctrl_combined_pval, ctrl_combined_pval, folder, name='ctrlCombStatC', sup=True)
        empirical_qval = empiricalPvalues(ctrl_combined_pval,global_pval[:,np.newaxis], folder, name='expCombStat{}'.format(combination), sup=True)
        
        f = open(os.path.join(folder, 'comb_empirical_p_qval{}{}.pkl'.format(combination,without_mean_persistence)), 'w')
        pickle.dump((ctrl_qval, empirical_qval),f); f.close()
    else:
        f=open(os.path.join(folder, 'all_distances_{}.pkl'.format(without_mean_persistence)))
        expL, siRNAL, geneL, global_result=pickle.load(f); f.close()
        
        f = open(os.path.join(folder, 'comb_empirical_p_qval{}{}.pkl'.format(combination,without_mean_persistence)), 'r')
        ctrl_qval, empirical_qval=pickle.load(f); f.close()

    empirical_qval=np.array(empirical_qval)
    exp_hit, gene_hit, gene_highconf, exp_of_highconfsiRNAs, siRNA_highconf=finding_hit(empirical_qval, threshold=threshold, siRNAL=list(siRNAL), geneL=list(geneL), expL=list(expL),
                                                                    mapping_filename=mapping_filename,
                                                                    trad=trad, without_mitotic_hits=without_mitotic_hits)
    if save:
        f=open(os.path.join(folder, '{}_hit_exp.pkl'.format(filename)), 'w')
        pickle.dump(exp_of_highconfsiRNAs,f)
        f.close()
        
    return empirical_qval,siRNAL, exp_hit,siRNA_highconf, exp_of_highconfsiRNAs, gene_highconf

def finding_hit(curr_qval,threshold, siRNAL, geneL, expL,mapping_filename,trad=True, without_mitotic_hits=False):
    exp_hit=[]
    siRNA_hit=[]
    gene_hit=[]
    siRNA_count=Counter(siRNAL)
    gene_count=Counter(geneL)
    
    for k in range(len(siRNAL)):
        if curr_qval[k]<threshold:
            exp_hit.append(expL[k])
            siRNA_hit.append(siRNAL[k])
            gene_hit.append(geneL[k])
                
    siRNA_hit=Counter(siRNA_hit)
    siRNA_hit_perc={siRNA:siRNA_hit[siRNA]/float(siRNA_count[siRNA]) for siRNA in siRNA_hit}
    
    siRNA_highconf = filter(lambda x: siRNA_hit_perc[x]>=0.5 and siRNA_hit[x]>=2, siRNA_hit_perc)
    
    if without_mitotic_hits:
        print "Filtering out mitotic hits from Mitocheck supp table 2"
        mito_hits = txtToList('../data/mitocheck_hits.txt')[:,0]
        siRNA_highconf=[siRNA for siRNA in siRNA_highconf if siRNA not in mito_hits]
    
    gene_highconf = Counter([geneL[siRNAL.index(siRNA)] for siRNA in siRNA_highconf])
    #gene_highconf={gene:gene_highconf[gene]/float(gene_count[gene]) for gene in gene_highconf}
    gene_highconf=filter(lambda x: gene_highconf[x]>=1, gene_highconf)
    
    print 'Experiences ', len(exp_hit), 'out of',len(expL), 'ie ', len(exp_hit)/float(len(expL))
    print 'siRNA high conf ', len(siRNA_highconf), 'out of',len(siRNA_count), 'ie', len(siRNA_highconf)/float(len(siRNA_count))
    print 'Genes high conf', len(gene_highconf), 'out of', len(gene_count), 'ie ', len(gene_highconf)/float(len(gene_count))

    exp_of_highconfsiRNAs=[]
    for siRNA in siRNA_highconf:
        experiments = np.array(expL)[np.where(np.array(siRNAL)==siRNA)]
        exp_of_highconfsiRNAs.extend([experiment for experiment in experiments if experiment in exp_hit])
    
    if trad:
        trad = EnsemblEntrezTrad(mapping_filename)
        gene_highconf_Ensembl=[trad[el] for el in gene_highconf]
        gene_Ensembl = [trad[el] for el in gene_count]
        
        for geneList in [gene_highconf_Ensembl, gene_Ensembl]:
            for i,gene in enumerate(geneList):
                    if '/' in gene:
                        geneList[i]=gene.split('/')[0]
                        geneList.append(gene.split('/')[1])

        geneListToFile(gene_highconf_Ensembl, 'gene_high_conf_KS.txt')
        geneListToFile(gene_Ensembl, 'gene_list_KS.txt')

    return exp_hit, gene_hit, gene_highconf, exp_of_highconfsiRNAs, siRNA_highconf
        
        

def hitDistances(folder,key_name = 'distances_whole_5Ctrl{}', filename='all_distances_whole_5Ctrl{}.pkl',
                 mapping_filename = '../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt',
                 qc_filename='../data/mapping_2014/qc_export.txt',
                 redo=False,
                 without_mean_persistence=True, features=None):
    '''
    This function collects all distances or p-values from files, then with or without renormalizing them with respect to control distributions
    combines them with Fisher's method and finally with or without renormalizing them computes the hit list
    
    '''
    
    #i. Collecting distances
        #parameters : names of quality control and mapping files
    
    exp = collectingDistances(filename, folder,key_name, 
                              testCtrl=False, qc_filename=qc_filename,mapping_filename=mapping_filename,
                              redo=redo)
    
    combined_pval = np.zeros(shape = (exp[exp.keys()[0]][-1].shape[0],len(parameters)), dtype=float)
    
    for i,param in enumerate(parameters):
        empirical_pval=exp[param][-1]
        if without_mean_persistence:
            empirical_pval = np.delete(empirical_pval, features.index('mean persistence'),1)

#ii. for each parameter, combining pvalues from KS using Fisher's formula
        stat = -2*np.sum(np.log(empirical_pval),1)
        combined_pval[:,i] = stat
#        
        siRNAL, expL, geneL, _ = exp[param]

#    result = zip(expL, geneL, siRNAL, list(combined_pval[:,0]))
#    
#    result.sort(key=itemgetter(0))
    return expL, geneL, siRNAL, combined_pval#np.array(result) 

# def collectingDistances_XB(folder, filename='distances_tw_', iter_range=range(5), use_time_window=True, compL=None, norm='neg_ctrl'):
#     '''
#     This function aims at collecting p-values from the xenobiotic screen. No need to look again for the experiments that are involved since it's saved in
#     the file
#     '''
#     files=os.listdir(folder)
#     if use_time_window:
#         time_window_range=range(1)#len(time_windows)
#         result={(i,j):None for (i,j) in product(iter_range, time_window_range)}
#     else:
#         result=[None for k in iter_range]
#     who=defaultdict(list); compounds=defaultdict(list); doses=defaultdict(list)
#     
#     if compL==None:
#         compL=compoundL
#     
#     for compound in compL:
#     #we get the list of files from dose 1 to dose max
#         cFiles=sorted(filter(lambda x: filename in x and compound in x, files))
#         if compound in CONTROLS.values() and norm=='neg_ctrl':
#             tag=(lambda x:x.split('_')[-2])
#             for file_ in cFiles:
#                 try:
#                     f=open(os.path.join(folder, file_), 'r')
#                     d=pickle.load(f); f.close()
#                 except:
#                     pdb.set_trace()
#                 else:
#                     for param_set in d:
#                         currParams=dict(param_set)
#                         try:
#                             who[currParams['time_window']].extend(currParams['wells'])
#                         except KeyError:
#                             pdb.set_trace()
#                         compounds[currParams['time_window']].extend([compound for k in range(d[param_set].shape[0])])
#                         doses[currParams['time_window']].extend([0 for k in range(d[param_set].shape[0])])
#                     #deleting mean persistence if it is still there
#                         if d[param_set].shape[1]==16:
#                             d[param_set]=np.delete(d[param_set], 15, 1)
#                         for iter_ in iter_range:
#                             if use_time_window:
#                                 try:
#                                     assert(len(currParams['wells'])==d[param_set].shape[0] )
#                                 except AssertionError:
#                                     pdb.set_trace()
#                                 result[(iter_,currParams['time_window'])]= d[param_set] if result[(iter_,currParams['time_window'])] is None\
#                                      else np.vstack((result[(iter_,currParams['time_window'])], d[param_set]))
#                             else:
#                                 raise AttributeError
#         else:
#             for file_ in cFiles:
#                 #ATTENTION ca ne marche pas dans le cas ou les doses vont de 1 a 10 pcq ds ce cas l'ordre est 1 10 2 3 4 5 6 7 8 9
#                 #tag = lambda x: cFiles.index(x)+1
#                 tag = lambda x: x.split('_')[-1].split('.')[0]
#                 print file_,
#                 
#                 try:
#                     f=open(os.path.join(folder, file_), 'r')
#                     d=pickle.load(f); f.close()
#                 except:
#                     pdb.set_trace()
#                 else:
#                     print len(d.keys())
#                     for param_set in d:
#                         currParams=dict(param_set)
#                         iter_ = currParams['iter']
#                         if iter_==0:
#                             try:
#                                 who[currParams['time_window']].extend(currParams['wells'])
#                             except KeyError:
#                                 pdb.set_trace()
#                             compounds[currParams['time_window']].extend([compound for k in range(d[param_set].shape[0])])
#                             doses[currParams['time_window']].extend([tag(file_) for k in range(d[param_set].shape[0])])
#                         #deleting mean persistence if it is still there
#                         if currParams['time_window']==0:
#                             if d[param_set].shape[1]==16:
#                                 d[param_set]=np.delete(d[param_set], 15, 1)
#                             if use_time_window:
#     #                             if np.any(np.isnan(d[param_set])):
#     #                                 pdb.set_trace()
#                                 result[(iter_, currParams['time_window'])]= d[param_set] if result[(iter_, currParams['time_window'])] is None\
#                                          else np.vstack((result[(iter_, currParams['time_window'])], d[param_set]))
#                             else:
#                                 raise AttributeError
#     print result[(0,0)].shape
#     result={tw: np.vstack((-np.log(result[(el,tw)])[np.newaxis] for el in range(5))) for tw in time_window_range}
#     result2={tw:2*np.sum(result[tw],2) for tw in result}
#     
#     doses={el:np.array(doses[el], dtype=int) for el in who}
#     compounds={el:np.array(compounds[el]) for el in who}
#     who={el:np.array(who[el]) for el in who}
#     conditions = {el:["{}_{}".format(a,b) for a,b in zip(compounds[el],doses[el])] for el in doses}
#     
#     return result[0], result2[0], who[0], compounds[0], doses[0], conditions[0]
# 
# def finding_hit_XB(result, who, compounds, doses, combination=(lambda x: np.max(x,0)), 
#                    outputFolder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions__settings2/features_on_films',
#                    plot=True,
#                    saveTextFile=True):
#     '''
#     Working on min(stat)=working on max(pval) == CONSERVATIVE approach
#     Working on max(stat)=working on min(pval) == OPTIMISTIC APPROACH
#     '''
#     result = np.array(combination(result))[:,np.newaxis]
#     
#     #do the comparison with DMSO distribution
#     control_stat=result[np.where(compounds=='DMSO')]
#     stat=result[np.where(np.array([compound not in ['Nonane', 'DMSO', 'TCDD'] for compound in compounds]))]
#     pval_DMSO, qval_DMSO = empiricalPvalues(control_stat, stat, outputFolder, 'DMSO_compared', sup=True, also_pval=True)
#     compounds_DMSO=compounds[np.where(np.array([compound not in ['Nonane', 'DMSO', 'TCDD'] for compound in compounds]))]
#     who_DMSO=who[np.where(np.array([compound not in ['Nonane', 'DMSO', 'TCDD'] for compound in compounds]))]
#     doses_DMSO=doses[np.where(np.array([compound not in ['Nonane', 'DMSO', 'TCDD'] for compound in compounds]))]
#     
#     conditions_DMSO=["{} {} {}".format(a,b,c) for a,b,c in zip(compounds_DMSO, doses_DMSO, who_DMSO)]
#     
#     
#     #now with Nonane for TCDD
#     control_stat=result[np.where(compounds=='Nonane')]
#     stat=result[np.where(compounds=='TCDD')]
#     pval_Nonane, qval_Nonane = empiricalPvalues(control_stat, stat, outputFolder, 'Nonane_compared', sup=True, also_pval=True)
#     
#     conditions_Nonane=["{} {} {}".format(a,b,c) for a,b,c in zip(['TCDD' for k in range(len(pval_Nonane))], doses[np.where(compounds=='TCDD')],\
#                                                                  who[np.where(compounds=='TCDD')])]
#     l=len(pval_DMSO)
#     conditions_DMSO.extend(conditions_Nonane); pval_DMSO=np.vstack((pval_DMSO,pval_Nonane)); qval_DMSO=np.vstack((qval_DMSO[:,np.newaxis],qval_Nonane[:,np.newaxis])); 
#     result=zip(pval_DMSO, qval_DMSO, conditions_DMSO);result.sort(key=itemgetter(0))
#     
#     if plot:
#         qval_DMSO=np.array(qval_DMSO)
#         for compound in ['BPA', 'MeHg', 'PCB', 'Endo']:
#             zz=[]
#             f=p.figure(figsize=(24,16))
#             ax=f.add_subplot(111)
#             for i,el in enumerate(qval_DMSO[np.where(compounds_DMSO==compound)[0]]):
#                 if plates.index(who_DMSO[np.where(compounds_DMSO==compound)][i,0]) not in zz:
#                     ax.scatter(i,el, marker=markers[plates.index(who_DMSO[np.where(compounds_DMSO==compound)][i,0])], \
#                                color=couleurs[plates.index(who_DMSO[np.where(compounds_DMSO==compound)][i,0])], label=who_DMSO[np.where(compounds_DMSO==compound)][i,0])
#                     zz.append(plates.index(who_DMSO[np.where(compounds_DMSO==compound)][i,0]))
#                 else:
#                     ax.scatter(i,el, marker=markers[plates.index(who_DMSO[np.where(compounds_DMSO==compound)][i,0])],\
#                                color=couleurs[plates.index(who_DMSO[np.where(compounds_DMSO==compound)][i,0])])
#                 ax.text(i,el,"{} {}".format(compounds_DMSO[np.where(compounds_DMSO==compound)][i][0],doses_DMSO[np.where(compounds_DMSO==compound)][i]))
#             ax.set_ylim(0,1)
#             p.legend(fontsize=7)
#             p.savefig(os.path.join(outputFolder, "pval_{}_comb_med.png".format(compound)))
#             
#         who_Nonane=who[np.where(compounds=='TCDD')];zz=[]
#         f=p.figure(figsize=(24,16))
#         ax=f.add_subplot(111)
#         for i,el in enumerate(qval_Nonane):
#             if who_Nonane[i,0] not in zz:
#                 ax.scatter(i,el, marker=markers[plates.index(who_Nonane[i,0])], color=couleurs[plates.index(who_Nonane[i,0])], label=who_Nonane[i,0])
#                 zz.append(who_Nonane[i,0])
#             else:
#                 ax.scatter(i,el, marker=markers[plates.index(who_Nonane[i,0])], color=couleurs[plates.index(who_Nonane[i,0])])
#             ax.text(i,el,"T {}".format(doses[np.where(compounds=='TCDD')][i]))
#         ax.set_ylim(0,1)
#         p.legend(fontsize=7)
#         p.savefig(os.path.join(outputFolder, "pval_TCDD_comb_med.png"))
#     if saveTextFile:
#         f=open(os.path.join(outputFolder, 'ranking_pval_med{}.txt'.format(datetime.date.today())), 'w')
#         f.write('P-val\t Q-val\t Xb\t Dose\t Plate\t Well\n')
#         for el in result:
#             f.write("{}\t {} \t {}\t {} \t {} \t {}\n".format(el[0][0], el[1][0], *el[2].split(" ")))
#         f.close()
#     return result, pval_DMSO[:l], qval_DMSO[:l], pval_Nonane,  qval_Nonane
#     
# 
# def plotDistances_XB(result, compounds,who, doses, combination=(lambda x: np.max(x,0)), cmap=mpl.cm.bwr, norm=2, compL=None, sh=True,
#                      outputFolder = '../../../projects/Xb_screen/dry_lab_results/MOTILITY/track_predictions__settings2/features_on_films',
#                      filename='heatmaps_pval_comb_med_medstandardized_TW0.png', length=None, normal_z_score=False, ctrl="neg_plate"
#                      ):
#     '''
#     So just plot-wise we look at the different -log(pvalues) for all features for the different conditions,
#     normalizing with the mean,std of the controls from all plates.
#     
#     It seems though that some plates have a wide range of pvalues for control-control comparisons,
#     which is not the case on other plates... So it might be interesting to look at things plate-wise
#     
#     By the way, if the combination here is min, it means it's max(pval) because result is supposed to be -log(pval),
#     so it's just the other way around with all mitocheck results. We chose in the latter to use max(pval) ie most
#     conservative approach.
#     
#     Possible values for ctrl: neg_plate for negative controls on the same plate, neg_all fr negative controls on all plates
#     and all for all wells on the same plate
#     '''
#     if result.shape[0] >500:
#         #It means we're dealing with feature tables and not statistic or pvalue tables.
#         if length ==None:
#             raise ValueError
#         else:
#             new_result=np.zeros(shape=(len(length), result.shape[1]))
#             for i in range(len(length)):
#                 new_result[i]=np.median(result[np.sum(length[:i]):np.sum(length[:i+1])],0)
#             combination=lambda x:x
#             
#             result=np.hstack((new_result[:,:len(featuresNumeriques)],new_result[:, featuresSaved.index('mean straight'), np.newaxis]))
#             
#     if compL==None:
#         compL=compoundL
#     if compL[0] not in xbL:
# #second plot for the controls
#         f, axes=p.subplots(1,4,sharey=True, figsize=(15,16))
#         L=compL
#     else:
# #first plot for the xenobiotics
#         f, axes=p.subplots(5,10,sharex=True, sharey=True,figsize=(24,16))
#         L=filter(lambda x: x in xbL, compL)
#         result=np.array(combination(result))
#     print result.shape
#     k=0
#     for i,compound in enumerate(L):
#         print compound
#         subR = result[np.where(compounds==compound)]
#         where_compound=np.where(compounds==compound)[0]
#         if compound in CONTROLS.keys():
#             where_1=np.where(compounds==CONTROLS[compound])[0]
#             for el,j in enumerate(where_compound):
#                 pl=who[j][0]; print "-------------------------",pl
#                 if ctrl=='neg_plate':
#                     where_control=filter(lambda x:x in where_1, np.where(who[:,0]==pl)[0])
#                 elif ctrl=="neg_all":
#                     where_control=where_1
#                 elif ctrl=="plate":
#                     where_control=np.where(who[:,0]==pl)[0]
#                 else:
#                     raise AttributeError
#                 
#                 if normal_z_score:
#                     sig=np.std(result[where_control],0)
#                     sig[np.where(sig==0)]+=0.001
#                     subR[el]=(subR[el]-np.mean(result[where_control],0))/sig
#                 else:
#                     IQR=(scoreatpercentile(a=result[where_control], per=75,axis=0)-scoreatpercentile(a=result[where_control], per=25,axis=0))
#                     print np.median(result[where_control],0)[8], IQR[8]
#                     IQR[np.where(IQR==0)]+=0.001
#                     subR[el]=(subR[el]-np.median(result[where_control],0))/IQR
#                 
#         else:
#             for el,j in enumerate(where_compound):
#                 pl=who[j][0]
#                 if ctrl=='neg_plate':
#                     where_control=filter(lambda x:x in where_compound, np.where(who[:,0]==pl)[0])
#                 elif ctrl=="neg_all":
#                     where_control=where_compound
#                 elif ctrl=="plate":
#                     where_control=np.where(who[:,0]==pl)[0]
#                 else:
#                     raise AttributeError
#                 
#                 if normal_z_score:
#                     sig=np.std(result[where_control],0)
#                     sig[np.where(sig==0)]+=0.001
#                     subR[el]=(subR[el]-np.mean(result[where_control],0))/sig
#                 else:
#                     IQR=(scoreatpercentile(a=result[where_control], per=75,axis=0)-scoreatpercentile(a=result[where_control], per=25,axis=0))
#                     print np.median(result[where_control],0)[8], IQR[8]
#                     IQR[np.where(IQR==0)]+=0.001
#                     subR[el]=(subR[el]-np.median(result[where_control],0))/IQR
#                 
#         subD = np.array(doses[np.where(compounds==compound)], dtype=int) 
#         
#         if compound in xbL:
#             dose_range=range(np.max(subD)-1) if (compound=='BPA' and np.max(subD)==10) else range(np.max(subD))
#             for dose in dose_range:
#                 subW=who[np.where((compounds==compound) & (doses==dose+1))]
#                 try:
#                     axes[i,dose].matshow(subR[np.where(subD==dose+1)[0]][:,np.newaxis], cmap=cmap, norm=mpl.colors.Normalize(-norm,norm))
#                 except ValueError:
#                     axes[i,dose].matshow(subR[np.where(subD==dose+1)[0]], cmap=cmap, norm=mpl.colors.Normalize(-norm,norm))
#                 for j in range(subR[np.where(subD==dose+1)].shape[0]):
#                     axes[i,dose].text(0,j,subW[j],fontsize=7)
#                 axes[i,dose].set_ylim(0,6)
#                 axes[i,dose].set_xlim(0,16)
#                 axes[i,dose].text(0,-5,"Dose {}".format(dose+1))
#                 axes[i,dose].set_xticklabels([])
#             axes[i,0].set_title(compound)
#         else:
#             subW =  np.array(who[np.where(compounds==compound)])
#             axes[k].set_ylim(0,27)
#             try:
#                 axes[k].pcolor(subR[:,np.newaxis], cmap=cmap, norm=mpl.colors.Normalize(-norm,norm))
#             except ValueError:
#                 axes[k].pcolor(subR, cmap=cmap, norm=mpl.colors.Normalize(-norm,norm))
#             for j in range(subR.shape[0]):
#                 axes[k].text(0,j,subW[j])
#             axes[k].set_title(compound); k+=1
#         
#     if sh:
#         p.show()
#     else:
#         p.savefig(os.path.join(outputFolder, filename))
#     if L!=compL:
#         return plotDistances_XB(result, compounds, who,doses, combination, cmap, norm, compL=filter(lambda x: x not in xbL, compL), sh=sh,\
#                                 outputFolder=outputFolder, filename='controls{}'.format(filename),normal_z_score=normal_z_score, ctrl=ctrl)
#     
# def plotDistances_TWEvol(result, who,compounds, doses, combination=(lambda x: np.max(x,0)), compL=None, sh=False,
#                     outputFolder = '../../../projects/Xb_screen/dry_lab_results/track_predictions__settings2/features_on_films',
#                     filename='stat_comb_max_time_windows_{}',
#                     legend=True,
#                     CTRL_ONLY=False
#                      ):
#     if compL is None:
#         compL=compoundL
#         
#     r_comb={el: combination(result[el]) for el in result}
#     
#     for compound in filter(lambda x: x not in ['DMSO', 'Nonane'], compL):
#         print compound        
#         f,axes=p.subplots(1,5, sharex=True, sharey=True,figsize=(24,16))
#         control=CONTROLS[compound]
#         for k in range(len(time_windows)):
#             if not CTRL_ONLY:
#                 ind=np.where(compounds[k]==compound)[0]
#                 subR=r_comb[k][ind[0]:ind[-1]]
#                 for i,el in enumerate(subR):
#                     if not np.isnan(el):
#         #This because there are wells for which some features gave nan pvalues : 'movement_type' for well 61 pl 201214
#         #It is due to util.listFileManagement.correct_Nan 
#                         axes[k].scatter(i, el, color=couleurs[plates.index(who[k][ind[0]+i,0])], marker=markers[plates.index(who[k][ind[0]+i,0])], label=who[k][ind[0]+i,0])                                                                                      
#                         axes[k].set_title(time_windows[k])
#                         axes[k].text(i,el,int(doses[k][ind[0]+i]))
#                         axes[k].set_ylim(0,max(1000, np.max(subR)))
#             
#                 indC=np.where(compounds[k]==control)[0]
#                 subR=r_comb[k][indC[0]:indC[-1]]
#                 for i,el in enumerate(subR):
#                     axes[k].scatter(ind[-1]-ind[0]+4, el, color='black', marker=markers[plates.index(who[k][indC[0]+i,0])])
#             else:
#                 indC=np.where(compounds[k]==control)[0]
#                 subR=r_comb[k][indC[0]:indC[-1]]
#                 for i,el in enumerate(subR):
#                     axes[k].scatter(i, el, color=couleurs[plates.index(who[k][indC[0]+i,0])], marker=markers[plates.index(who[k][indC[0]+i,0])],label=who[k][indC[0]+i,0])
#                     axes[k].text(i,el,int(who[k][indC[0]+i,1]))
#                 axes[k].set_ylim(0,350)
#         if legend:
#             p.legend(fontsize=10)
#         if sh:
#             p.show()
#         else:
#             p.savefig(os.path.join(outputFolder, filename.format(compound)))

def collectingDistances(filename, folder, 
                        key_name = 'distances_length_5CtrlC',
                        qc_filename='../data/mapping_2014/qc_export.txt',mapping_filename='../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt', 
                        testCtrl =False,
                        iter_range=range(5),
                        redo=False, siRNAFilterList=None,long_version=False, usable_file='../resultData/features_on_films/usable_experiments_whole_mitocheck.pkl'):
    '''
    Function to collect p-values from applying CellExtractor to the Mitocheck data, hardly generalizing to the xenobiotic screen data (sorry Will) 
    '''
    
    
    global parameters
    if filename not in os.listdir(folder) or redo:
        if not testCtrl:
            files = filter(lambda x: key_name in x and 'CTRL' not in x and 'all' not in x, os.listdir(folder))
            yqualDict=expSi(qc_filename, sens=0)
            dictSiEntrez=siEntrez(mapping_filename, yqualDict.keys())
#             if long_version:
#                 f=open(usable_file)
#                 usable=pickle.load(f); f.close()
        else:
            files = filter(lambda x: key_name in x and 'CTRL' in x  and 'all' not in x, os.listdir(folder))
            
        print len(files)
        files.sort()
        result=[None for k in iter_range]
        who=defaultdict(list)
        siRNAList=[]
        genes=[]
        
        for file_ in files:
            try:
                f=open(os.path.join(folder, file_))
                d=pickle.load(f)
                f.close()
            except OSError, IOError:
                print "Issue with ", file_
                continue
            if not testCtrl:
                siRNA = file_.split('_')[-1][:-4]
                shape = list(set([d[el].shape[0] for el in d]))
                try:
                    assert(len(shape)==1)
                except:
                    if len(d)==6:#meaning we have one time iteratiom 0 with less wells
                        print [d[el].shape[0] for el in d]
                        key_ = [el for el in d if d[el].shape[0]==np.min(shape)][0]
                        del d[key_]
                        f=open(os.path.join(folder, file_), 'w')
                        pickle.dump(d,f); f.close()
                        
                shape = np.max(shape)
                siRNAList.extend([siRNA for k in range(shape)])
                gene = dictSiEntrez[siRNA]
                genes.extend([gene for k in range(shape)])
                
            else:
                siRNA=file_[-13:-4]

            for param_set in d:
                currParams = dict(param_set)
                
                l= Counter(np.where(~np.isnan(d[param_set]))[0]).keys()
                if l==[]:
                    continue
                else:
#                    if len(l)==1 and not testCtrl:
#                        #meaning it is an orphan siRNA with only one experiment. But do we really want to take them out? The experiment is real after all
#                        continue 
                    who[currParams['iter']].extend(currParams['wells'])
                    result[currParams['iter']] = d[param_set] if result[currParams['iter']]==None else np.vstack((result[currParams['iter']], d[param_set]))
        
#         if not testCtrl:
#             result = np.
        
        f=open(os.path.join(folder, filename), 'w')
        pickle.dump((result, who, siRNAList, genes), f); f.close()

    else:
        f=open(os.path.join(folder, filename))
        result=pickle.load(f); f.close()
    
    return result

def empiricalPvalues(dist_controls, dist_exp, folder, name, sup=False, also_pval=False):
    empirical_pval = np.zeros(shape = (dist_exp.shape[0], dist_exp.shape[1]), dtype=float)
    empirical_qval = np.zeros(shape = (dist_exp.shape[0], dist_exp.shape[1]), dtype=float)

    for k in range(dist_exp.shape[1]):
        #for all features
        print '---------------',k
        for j in range(dist_exp.shape[0]):
            #for all experiments
            
            if sup:
#if we are dealing with distances to calculate empirical distributions,
#we're looking for distances that are big, so a small p-value should be for big distances.
#Hence our empirical p-values are 100 - the percentage of distances that are lower
                empirical_pval[j,k] = (100 - percentileofscore(dist_controls[:,k], dist_exp[j,k]))/100.0
            else:
                empirical_pval[j,k] = percentileofscore(dist_controls[:,k], dist_exp[j,k])/100.0
                
#The minimum precision is 1/#controls. It is not possible that there is p-value 0 here        
        empirical_pval+=1/float(dist_controls.shape[0])
#Calling R function p-adjust, applying Benjamini-Hochberg method
        empirical_qval = stats.p_adjust(FloatVector(list(empirical_pval[:,k])), method = 'BH')
    #plotting empirical p-value distribution
    f=p.figure()
    ax=f.add_subplot(111)
    nb,bins,patches=ax.hist(empirical_pval, bins=50, alpha=0.2)
    nb,bins,patches=ax.hist(empirical_qval, bins=50, alpha=0.2, color='red')
    ax.set_xlim(0,1)
    p.savefig(os.path.join(folder, 'pvalParamKS_{}.png'.format(name)))
    p.close('all')
    if also_pval:
        return np.array(empirical_pval), np.array(empirical_qval)
    
    return np.array(empirical_qval)

#def empiricalDistributions(dist_controls, dist_exp, folder,iteration=1, sup=False, union=False, redo=False, renorm_first_statistic=False, renorm_second_statistic=True,
#                           without_mean_persistence=True):
#    '''
#    This function computes empirical p-value distributions
#    dist_controls et dist_exp are dictionaries with keys=parameter sets
#    dist_controls[param] array of size nb experiments x nb of features
#    
#    Parameter without_mean_persistence says if we are going to take into account the mean persistence feature
#    
#    '''
##here we define what the features are
#    features = list(featuresNumeriques)
#    features.append('mean persistence')
#    features.append('mean straight')
#    
#    param=parameters[0]
#    print "First, univariate p-values"
#    if renorm_first_statistic:
#        if 'empirical_p_qval.pkl' not in os.listdir(folder) or redo:
#            ctrl_pval, ctrl_qval = empiricalPvalues(dist_controls, dist_controls, folder,name='ctrlPval{}'.format(iteration), sup=sup)
#            empirical_pval, empirical_qval = empiricalPvalues(dist_controls, dist_exp, folder,name='expPval{}'.format(iteration), sup=sup)
#            
#            
#            f = open(os.path.join(folder, 'empirical_p_qval{}.pkl'.format(iteration)), 'w')
#            pickle.dump((ctrl_pval, ctrl_qval, empirical_pval, empirical_qval),f); f.close()
#        else:
#            f=open(os.path.join(folder, 'empirical_p_qval{}.pkl'.format(iteration)))
#            ctrl_pval, ctrl_qval, empirical_pval, empirical_qval = pickle.load(f); f.close()
#    else:
#        ctrl_pval = dist_controls; empirical_pval = dist_exp
#    
#    print 'Second, combining p-values'
#    # statistical test according to Fisher's method http://en.wikipedia.org/wiki/Fisher%27s_method
#    ctrl_combined_pval = {param : np.zeros(shape = (ctrl_pval[param].shape[0],1), dtype=float) for param in parameters}
#    combined_pval = {param : np.zeros(shape = (empirical_pval[param].shape[0],1), dtype=float) for param in parameters}
#
#    ctrl_combined_qval = {param : np.zeros(shape = (ctrl_pval[param].shape[0],), dtype=float) for param in parameters}
#    combined_qval = {param : np.zeros(shape = (empirical_pval[param].shape[0],), dtype=float) for param in parameters}
#    for param in parameters:
#        if without_mean_persistence:
#            ctrl_pval[param]=np.delete(ctrl_pval[param], features.index('mean persistence'),1)
#            empirical_pval[param] = np.delete(empirical_pval[param], features.index('mean persistence'),1)
#            
#        pdb.set_trace()
#            
#        stat = -2*np.sum(np.log(ctrl_pval[param]),1) 
#        ctrl_combined_pval[param] = stat[:,np.newaxis]
#        
#        stat = -2*np.sum(np.log(empirical_pval[param]),1)
#        combined_pval[param] = stat[:,np.newaxis]
#    
#    
#    if renorm_second_statistic:
#        if 'comb_empirical_p_qval{}.pkl'.format(iteration) not in os.listdir(folder) or redo:
#            ctrl_combined_pval2, ctrl_combined_qval = empiricalPvalues(ctrl_combined_pval, ctrl_combined_pval, folder,name='ctrlCombinedStat{}'.format(iteration), sup=True)
#            combined_pval, combined_qval = empiricalPvalues(ctrl_combined_pval, combined_pval, folder,name='expCombinedStat{}'.format(iteration), sup=True)
#            
#            f = open(os.path.join(folder, 'comb_empirical_p_qval{}.pkl'.format(iteration)), 'w')
#            pickle.dump((ctrl_combined_pval2, ctrl_combined_qval, combined_pval, combined_qval),f); f.close()
#        else:
#            f = open(os.path.join(folder, 'comb_empirical_p_qval{}.pkl'.format(iteration)), 'r')
#            ctrl_combined_pval2, ctrl_combined_qval, combined_pval, combined_qval=pickle.load(f); f.close()
#
#    else:
#        for param in parameters:
#            ctrl_combined_pval[param] = chi2.sf(ctrl_combined_pval[param], 2*ctrl_pval[param].shape[1])
#            ctrl_combined_qval[param]= ctrl_combined_pval[param]*ctrl_combined_pval[param].shape[0]/(1+np.argsort(ctrl_combined_pval[param]))
#            
#            combined_pval[param] = chi2.sf(combined_pval[param], 2*empirical_pval[param].shape[1])
#            combined_qval[param]= combined_pval[param]*combined_pval[param].shape[0]/(1+np.argsort(combined_pval[param]))
#    print 'End'
#    return ctrl_combined_pval2, ctrl_combined_qval, combined_pval, combined_qval


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
                 interestFeatures = ['length'],
                 testCtrl = False,
                 div_name='total_variation', 
                 time_window=None,
                 iter_=0,
                 verbose=0):
        
        '''
MITOCHECK data: if doing ctrl-ctrl p-values, the plate that is under study is in testCtrl. In this case, self.siRNA takes values CTRL_[plate] 
Otherwise testCtrl is 0 and self.siRNA is the siRNA under study

XB SCREEN data: if doing ctrl-ctrl p-values, the plate that is under study is in testCtrl. In this case self.siRNA takes values CTRL_[plate]_solvent
Otherwise testCtrl is 0 and self.siRNA is xenobiotic_dose. If we are going for plate normalization, as opposed to negative control normalization,
self.siRNA takes value CTRL_[plate]_plate
    '''

        self.siRNA = siRNA
        self.settings = settings.Settings(settings_file, globals())
        self.inputFile = self.settings.inputFile if time_window== None else self.settings.filename_twindow.format(time_window=time_window)
        
        self.time_window = time_window
        self.verbose=verbose
        #default value: we compare experiments to negative controls on the same plate
        self.ctrl='neg_ctrl'
        
        if not testCtrl:
            self.plate = None
            if self.settings.xb_screen:
                if self.settings.norm=='neg_ctrl':
                    self.ctrl = CONTROLS[siRNA.split('_')[0]]
                elif self.settings.norm=='plate':
                    self.ctrl='plate'
                else:
                    raise ValueError
                print "Control ", self.ctrl
                
        else:
            if self.verbose:
                print "Testing controls for plate {}".format(testCtrl)
            if self.settings.xb_screen==False:
                assert (testCtrl in os.listdir(self.settings.data_folder))
                self.plate = testCtrl
                self.siRNA = 'CTRL_{}'.format(self.plate[:9])
            else:
                self.ctrl=testCtrl.split('_')[1]
                self.plate=testCtrl.split('_')[0]
                assert (self.plate in os.listdir(self.settings.data_folder))
                self.siRNA = 'CTRL_{}'.format(testCtrl)
            
        self.div_name =div_name
#        self.cost_type=cost_type
#        self.bins_type = bin_type
#        self.bin_size = bin_size
#        self.lambda_=lambda_   
        self.iter=iter_     
        if interestFeatures is None:
            self.currInterestFeatures = list(featuresNumeriques)
            self.currInterestFeatures.extend(['mean straight'])
            self.featuresSaved=featuresSaved
        else:
            self.currInterestFeatures=interestFeatures
            self.featuresSaved=interestFeatures
            
        if self.settings.histDataAsWell:
            raise AttributeError
    
    def _findExperiment(self):
        '''
        Purpose : find experiments that are related to this siRNA
        '''
        if self.settings.xb_screen==False:
            yqualDict=expSi(self.settings.quality_control_file, sens=0)
            l= strToTuple(yqualDict[self.siRNA], os.listdir(self.settings.data_folder))
        else:
            xb=self.siRNA.split('_')[0]; dose=self.siRNA.split('_')[1]
            try:
                f=open(self.settings.ok_wells_asDICT)
                info_dict=pickle.load(f); f.close()
                l=[(pl, "{:>05}".format(w)) for pl, w in info_dict[xb][int(dose)]]
            except:
                info_dict, _=fromXBToWells([xb],confDir=self.settings.plate_setups_folder, dose_filter=dose, verbose=self.verbose)
                l=info_dict[xb][dose]
            
        if self.verbose>0:
            print l
        return l
    
    def _concatenation(self, expList):
        '''
        This function returns the matrix of trajectory features of a given experiment list
        '''
        if self.settings.xb_screen==False:
            return histConcatenation(self.settings.data_folder, expList, self.settings.mitocheck_file,
                                    self.settings.quality_control_file,
                                    filename=self.inputFile,
                                    features=self.currInterestFeatures,
                                    min_size=self.settings.min_size,
                                    verbose=self.verbose, perMovie=True, hist=self.settings.histDataAsWell)
        else:
            return xbConcatenation(os.path.abspath(os.path.join(self.settings.data_folder, os.pardir)), expList, xb_list='processedDictResult_P{}.pkl', 
                                   visual_qc=self.settings.visual_qc,
                                   flou_qc=self.settings.flou_qc, 
                                   filename=self.inputFile,
                                   verbose=self.verbose, perMovie = True, 
                                   track_folder="track_predictions__settings2")
            
    def _usableControl(self, indication=None):
        '''
        This function returns the list of control experiment for a given plate, that have passed the quality control.
        For the xb screen data, it will probably be here that we are going to add the simulated control experiments
        '''
        if self.settings.xb_screen==False:
            ctrl = np.array(appendingControl([self.plate]))
            return ctrl[np.where(usable_MITO(self.settings.data_folder, 
                                ctrl, 
                                filename=self.inputFile,
                                features=self.currInterestFeatures,
                                min_size=self.settings.min_size,
                                qc=self.settings.quality_control_file, 
                                mitocheck=self.settings.mitocheck_file))]
        else:
            #we have put in self.siRNA CTRL_plate_control because there is DMSO and Nonane for the different xbs
            if self.ctrl!='plate':
                return usable_XBSC(self.ctrl, 0, self.plate)
            else:
                if indication==None:
                    try:
                        f=open(self.settings.ok_wells_asLIST)
                        a=pickle.load(f); f.close()
                    except:
                        a=np.array(computingToRedo()[0])
                    return [(pl, "{:>05}".format(w)) for pl, w in a[np.where(a[:,0]==self.plate)]]
                elif indication=='all_ctrl':
                    a=[]
                    try:
                        f=open(self.settings.ok_wells_asDICT)
                        d=pickle.load(f); f.close()
                    except:
                        d=None
                        
                    for control in Counter(CONTROLS.values()).keys():
                        try:
                            arr=np.array(d[control][0])
                            a.extend([(pl, "{:>05}".format(w)) for pl, w in arr[np.where(arr[:,0]==self.plate)]])
                        except:
                            a.extend(usable_XBSC(control, 0, self.plate))
                    return a
                else:
                    raise ValueError
        
    
    def getData(self, histDataAsWell):
        '''
        Ici sur toutes les experiences liees au siRNA self.siRNA on construit l'histogramme des features numeriques
        R contains all the information for all experiments, and length is the list of the length of all experiments,
        meaning the number of trajectories in each
        '''
        histDict = defaultdict(list)
        
        r, self.expList,_, length, _, _, _ = self._concatenation(self.expList)
        
        for i in range(len(length)):
            for k,feature in enumerate(self.currInterestFeatures):
                histDict[feature].append(r[np.sum(length[:i]):np.sum(length[:i+1]),self.featuresSaved.index(feature)])
        
#        f=open(os.path.join(self.settings.result_folder, 'distExp_ctrl_{}_{}.pkl'.format(self.bins_type, self.bin_size)))
#        bins = pickle.load(f); f.close()
        
        return histDict, None
    
    def plateNormCtrlHistograms(self, bins, ctrlExpList):
        histDict = defaultdict(list)
        length=[]

        if self.verbose:
            print 'Controls plate-wise', len(ctrlExpList)
        try:
            curr_r, _,_, curr_length, _, _, _ = self._concatenation(np.array(ctrlExpList))
        except:
            print "Problem with controls from plate {}".format(self.plate)
            length.append(np.nan)
            for k, feature in enumerate(self.currInterestFeatures):
                histDict[feature].append(np.nan)
        else:    
            length.append(np.sum(curr_length))
            assert(np.sum(curr_length)==curr_r.shape[0])                        
            for k, feature in enumerate(self.currInterestFeatures):
                histDict[feature].append(curr_r[:,self.featuresSaved.index(feature)])
                
        assert(len(length)==1)
        assert(len(histDict[feature])==1)
        return histDict
    
    def ctrlHistograms(self, bins,usable_ctrl=None, toDel=[]):
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
                if self.verbose:
                    print 'before cleaning', len(usable_ctrl)
                
                true_ctrl = filter(lambda x: x not in toDel, range(len(usable_ctrl)))
                ctrlExpList = list(np.array(usable_ctrl)[true_ctrl])
                if self.verbose:
                    print 'after cleaning', len(ctrlExpList)
                    
            #saving the controls that are usable, and five permutations of range(len(usable_ctrl)) for the experiments
                permutations = [np.random.permutation(len(usable_ctrl))[:min(len(usable_ctrl)-1,5)] for k in range(5)]
                if not self.settings.xb_screen:
                    f=open(os.path.join(self.settings.result_folder, self.settings.ctrl_exp_filename.format(plate)), 'w')
                    pickle.dump([np.array(usable_ctrl), permutations], f); f.close()
                elif self.settings.ctrl_exp_filename.format(plate, self.ctrl) not in os.listdir(self.settings.result_folder):
            #this ensure that we use the same order of permutations for all control pvalue list across the different time windows.
            #Hence we can directly compare the pvalue lists across time-windows
                    print "Defining usable controls and permutations"
                    f=open(os.path.join(self.settings.result_folder, self.settings.ctrl_exp_filename.format(plate, self.ctrl)), 'w')
                    pickle.dump([np.array(usable_ctrl), permutations], f); f.close()
                
            else:
            #loading the controls that were used to compute control control p-values
                try:
                    if not self.settings.xb_screen:
                        f=open(os.path.join(self.settings.result_folder, self.settings.ctrl_exp_filename.format(plate)))
                        usable_ctrl, permutations = pickle.load(f); f.close()
                    else:
                        f=open(os.path.join(self.settings.result_folder, self.settings.ctrl_exp_filename.format(plate, self.ctrl)))
                        usable_ctrl, permutations = pickle.load(f); f.close()
                    ctrlExpList=usable_ctrl[permutations[self.iter]]
                except IOError:
                    sys.stderr.write('No file registering used ctrls for plate {}'.format(plate))
                    ctrlExpList=None
    #DONC LES CONTROLES SONT LOADED SANS LES nan POUR TOUTES LES FEATURES POUR 5Ctrl1 je relance 5Ctrl2 AVEC les nan pour les features ou ils ne sont pas nan
            try:
                curr_r, _,_, curr_length, _, _, _ = self._concatenation(np.array(ctrlExpList))
            except:
                print "Problem with controls from plate {}".format(plate)
                length.append(np.nan)
                for k, feature in enumerate(self.currInterestFeatures):
                    histDict[feature].append(np.nan)
            else:    
                length.append(np.sum(curr_length))
                assert(np.sum(curr_length)==curr_r.shape[0])                        
                for k, feature in enumerate(self.currInterestFeatures):
                    histDict[feature].append(curr_r[:,self.featuresSaved.index(feature)])
        assert(len(length)==len(plates))
        assert(len(histDict[feature])==len(plates))
        return plates, histDict
    
    def parameters(self):
        r={
           'div_name':self.div_name,
           'iter':self.iter,
           'wells':tuple(self.expList)
           }
        if self.time_window is not None:
            r['time_window']=self.time_window
        if self.verbose:
            print "Parameters: iter", self.iter, "time window ", self.time_window 
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
                h1 = histogrammes[feature][i]
                h2 = ctrl_histogrammes[feature][corresponding_ctrl]
                if self.div_name =='KS':
                    if not (type(h2)==float and np.isnan(h2)):
                #taking the p-value because the distance does not take into account the sample sizes explicitly
                        distances[i,k]=ks_2samp(h1[~np.isnan(h1)], h2[~np.isnan(h2)])[1]
                    else:
                        pdb.set_trace()
                        distances[i,k]=np.nan
                elif self.div_name =='MW':
                    if not (type(h2)==float and np.isnan(h2)):
                #taking the p-value because the distance does not take into account the sample sizes explicitly
                        distances[i,k]=mannwhitneyu(h1[~np.isnan(h1)], h2[~np.isnan(h2)])[1]*2
                    else:
                        distances[i,k]=np.nan
                else:
                    raise AttributeError

        return distances
    
    def alreadyDone(self):
        if self.settings.outputFile.format(self.siRNA) in os.listdir(self.settings.result_folder):
            f=open(os.path.join(self.settings.result_folder, self.settings.outputFile.format(self.siRNA)))
            d=pickle.load(f); f.close()
            if self.parameters() in d:
                return True
            else:
                return False
        else:
            return False
    
    def saveResults(self, distances):
        if self.plate is None:
            if self.settings.outputFile.format(self.siRNA) in os.listdir(self.settings.result_folder):
                f=open(os.path.join(self.settings.result_folder, self.settings.outputFile.format(self.siRNA)))
                d=pickle.load(f); f.close()
                d.update({self.parameters():distances})
            else:
                d={self.parameters():distances}
        else:
            if self.settings.outputFile.format(self.siRNA) in os.listdir(self.settings.result_folder):
                f=open(os.path.join(self.settings.result_folder, self.settings.outputFile.format(self.siRNA)))
                d=pickle.load(f)
                f.close()
                if self.parameters() not in d:
                    d[self.parameters()]=None
            else:
                d={self.parameters():None}
            d[self.parameters()]=np.vstack((d[self.parameters()], distances)) if d[self.parameters()] is not None else distances

        f=open(os.path.join(self.settings.result_folder, self.settings.outputFile.format(self.siRNA)), 'w')
        pickle.dump(d, f)
        f.close()

        return
    
    def __call__(self):
        if not os.path.isdir(self.settings.result_folder):
            os.mkdir(self.settings.result_folder)
    #i. getting experiments corresponding to siRNA if we're not looking for control experiments only
        if self.plate is None:
            self.expList=self._findExperiment()
            #before anything, testing for existence if not redo anyway
            if not self.settings.redo and self.alreadyDone():
                print 'Already done'
                return
            if self.expList==[]:
                sys.stderr.write("No experiments for siRNA {}".format(self.siRNA))
                return
            #ii.calculate the histograms and binnings of experiments, using pre-computed binnings, ON THE EXPERIMENTS ONLY
            histogrammes, bins = self.getData(self.settings.histDataAsWell)
        
            #iii. calculate the histograms of corresponding controls. Should be per plate then, a dictionary
            plates, ctrl_histogrammes = self.ctrlHistograms(bins)
            
            #iv. calculate the distance from each experiment to its control
            distances = self.calculateDistances(plates, histogrammes, ctrl_histogrammes)
            if self.verbose:
                print distances
            #v. save results
            self.saveResults(distances)
            
            
        else:
            usable_ctrl = self._usableControl() 
            print len(usable_ctrl)
            if self.ctrl !='plate':
                if len(usable_ctrl)<2:
                    return
                elif len(usable_ctrl)==2:
                    #only one comparison possible if only two control movies
                    different_controls=[[0]]
                    
                elif len(usable_ctrl)<=6:
                    #randomly selecting one well of the plate that will be used to be compared to the others, and do it once
                    different_controls=[[i] for i in range(len(usable_ctrl))]
                    
                else:
                    #randomly selecting two wells of the plate that will be used to be compared to the others, and do it twice
                    different_controls=np.array([[a,b] for a,b in combinations(range(len(usable_ctrl)), 2)])[np.random.permutation(len(usable_ctrl)*(len(usable_ctrl)-1)/2)[:21]]
                    
                    
                if self.verbose:
                    print different_controls
    
                for false_experiments in different_controls:
                    self.expList = [usable_ctrl[i] for i in false_experiments]
                    if not self.settings.redo and self.alreadyDone():
                        print 'Already done'
                        return
    
                    
                    if self.verbose:
                        print 'false experiments', self.expList
    
                    histogrammes, bins = self.getData(self.settings.histDataAsWell)
                    
                    if len(usable_ctrl)==8:
                        left=[j for j in range(len(usable_ctrl)) if j not in false_experiments][np.random.permutation(len(usable_ctrl)-len(false_experiments))[0]]
                        false_experiments=np.hstack((false_experiments,left))
                        
                    if self.verbose:
                        print 'to del on this plate ', false_experiments
                        
                    plates, ctrl_histogrammes = self.ctrlHistograms(bins,usable_ctrl, toDel=false_experiments)
                        
                    #ii.calculate the histograms and binnings of experiments, using pre-computed binnings, ON all control experiments 
                    #for the plate
            
                    #iv. calculate the distance from each experiment to its control
                    distances = self.calculateDistances(plates, histogrammes, ctrl_histogrammes)
                    if self.verbose:
                        print distances
                    #v. save results
                    self.saveResults(distances)
            else:
            #in this case in usable_ctrl is the list of usable experiments on the plate
                different_controls = np.array([np.random.permutation(len(usable_ctrl))[:self.settings.n_n] for k in range(5)])
                print "Defining plate-wise controls and permutations"
                f=open(os.path.join(self.settings.result_folder, self.settings.ctrl_exp_filename.format(self.plate, self.ctrl)), 'w')
                pickle.dump([np.array(usable_ctrl), different_controls], f); f.close()

                self.expList=self._usableControl('all_ctrl')
                histogrammes, bins = self.getData(self.settings.histDataAsWell)
                
                if self.verbose:
                    print self.expList
                
                distances=None
                for currNormalization in different_controls:
                    if not self.settings.redo:
                        if self.alreadyDone():
                            print 'Already done'
                            return
                    #ii.calculate the histograms and binnings of experiments, using pre-computed binnings, ON all control experiments 
                    #for the plate                        
                    ctrl_histogrammes = self.plateNormCtrlHistograms(bins,np.array(usable_ctrl)[currNormalization])
                     
                    #iv. calculate the distance from each experiment to its control
                    currDistances = self.calculateDistances([self.plate], histogrammes, ctrl_histogrammes) 
                    distances = currDistances if distances ==None else np.vstack((distances, currDistances))
                    if self.verbose:
                        print currDistances
                #v. save results
                self.saveResults(currDistances)

        
        return
    
if __name__ == '__main__':
    '''
    Idea here is to calculate the distance of an experiment histogram to its control, for each numerical feature
    nb_exp should be the minimum number of experiments to get a stable distribution of the feature among
    considered experiments, so that the calculation is parallelizable
    
    -siRNA: the experimental condition, real siRNA for Mitocheck data, or xb_dose for xb screen data
    - plate: indicate the plate for sutyding ctrl-ctrl relations on this plate, else leave at None
    - ctrl: if dealing with xb screen data, indicate the control solvent between DMSO and Nonane
    
    
    '''
    parser = OptionParser(usage="usage: %prog [options]") 
    parser.add_option('--settings_file', type=str, default='tracking/settings/settings_cycle_lengthpval_extraction.py')   
    parser.add_option('--action', type=str, default=None)
    parser.add_option('--siRNA', type=str, dest='siRNA', default=None)
    parser.add_option('--testCtrl', type=str, dest='testCtrl', default=0)
    parser.add_option('--solvent', type=str, dest='solvent', default=None)
    parser.add_option('--div_name', type=str, dest='div_name', default='KS')
#    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
#    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
#    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)
    parser.add_option('--iter', type=int, dest="iter", default=0)
    parser.add_option('--time_window', type=int, dest="time_window", default=None)
#    parser.add_option("-l",type=int, dest="lambda_", default=10)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    (options, args) = parser.parse_args()
    
    print options.time_window
    if options.action=='collectDistances':
        collectingDistances('all_distances_whole_5Ctrl2.pkl', '../resultData/features_on_films/', 
                            '../data/qc_export.txt', '../data/mitocheck_siRNAs_target_genes_Ens75.txt', testCtrl=False, redo=True,
                            long_version=True)
        
    elif options.action=='collectSelectedDistances':
        f=open('../resultData/features_on_films/siRNAhighconf_5Ctrl2.pkl', 'r')
        siRNAFilterList =pickle.load(f)
        f.close()
        collectingDistances('all_distances_whole_5Ctrl_highconfSubset.pkl', '../resultData/features_on_films/', 
                            '../data/qc_export.txt', '../data/mitocheck_siRNAs_target_genes_Ens75.txt', testCtrl=False, redo=True,
                            long_version=True,
                            siRNAFilterList=siRNAFilterList)
        
    else:
        if options.solvent is not None and options.testCtrl:
#in case we're dealing with xb screen data we need to specify which solvent we want to study
            testCtrl='{}_{}'.format(options.testCtrl, options.solvent)
        else:
            testCtrl=options.testCtrl
        extractor=cellExtractor(options.siRNA, options.settings_file,testCtrl=testCtrl,div_name= options.div_name,
                                 time_window=options.time_window,
                                 iter_=options.iter,
                                 verbose=options.verbose)
        extractor()
        
#redone experiments after hdf5 catastrophe
#redone_experiments =np.array([['123438', 'LT0061_03--233'],
#       ['123438', 'LT0061_05--233'],
#       ['148427', 'LT0061_03--284'],
#       ['148427', 'LT0061_05--284'],
#       ['122325', 'LT0061_03--295'],
#       ['122325', 'LT0061_05--295'],
#       ['28902', 'LT0079_05--131'],
#       ['28902', 'LT0079_13--131']], 
#      dtype='|S14')
##experiments for which hdf5 are fine but the controls on the same plate are not
#tofilter_experiments=['LT0007_26--001',
# 'LT0003_15--194',
# 'LT0003_15--270',
# 'LT0003_15--271',
# 'LT0006_44--234',
# 'LT0006_44--236',
# 'LT0006_44--252',
# 'LT0007_26--152',
# 'LT0007_26--174',
# 'LT0010_27--080',
# 'LT0010_27--096',
# 'LT0010_27--097',
# 'LT0010_27--098',
# 'LT0010_27--099',
# 'LT0015_22--039']


    #        if filtering_bad_hdf5:
    #            f=open(os.path.join('../data/wrong_trajectories.pkl'))
    #            wrong_hdf5_files=pickle.load(f); f.close()
    #            
    #            wrong_hdf5_files=['{}--{}'.format(el[0][:9], el[1][2:5]) for el in wrong_hdf5_files]
    #            wrong_hdf5_files.extend(tofilter_experiments)
    #            #filtering bad hdf5 files
    #            toDel=[np.where(expL==wrong_file)[0][0] for wrong_file in wrong_hdf5_files if wrong_file in expL]
    #
    #            siRNAL=np.delete(siRNAL, toDel,0)            
    #            expL=np.delete(expL, toDel, 0)
    #            geneL=np.delete(geneL, toDel, 0)
    #            global_result=np.delete(global_result, toDel,0)
    #            
    #            l=Counter(siRNAL)
    #            if filtering_lonely_siRNAs:
    #                toDel=[np.where(np.array(siRNAL)==sirna)[0][0] for sirna in filter(lambda x: l[x]==1,l)]
    #                siRNAL=np.delete(siRNAL, toDel,0)
    #                expL=np.delete(expL, toDel, 0)
    #                geneL=np.delete(geneL, toDel, 0)
    #                global_result=np.delete(global_result, toDel,0)
    #
    #    
    #            #changing values for redone experiments
    #            f=open('../resultData/features_on_films/all_distances_whole_5CtrlC.pkl')
    #            redone=pickle.load(f)[parameters[0]];f.close()
    #            new_expL = np.array(redone[1])
    #            new_arr=-2*np.sum(np.log(redone[-1]),1)
    #            for redone_exp in redone_experiments[:,1]:
    #                for k in range(global_result.shape[1]):
    #                    global_result[np.where(expL==redone_exp),k]=new_arr[np.where(new_expL==redone_exp)]
