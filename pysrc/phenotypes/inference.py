import os, pdb
from util.listFileManagement import expSi, siEntrez, multipleGeneListsToFile,\
    EnsemblEntrezTrad
from collections import Counter
import numpy as np
import cPickle as pickle
from analyzer.rank_product import computeRPpvalues

from drug_screen_utils import CLASSES, plotInternalConsistency, KNOWN_TARGETS
from scipy.spatial.distance import cdist
from scipy.stats.mstats_basic import scoreatpercentile
from operator import itemgetter
from scipy.stats.stats import pearsonr
from itertools import product
from _collections import defaultdict
from util import hierarchical_clustering
import matplotlib as mpl
from phenotypes import drug_screen_utils

#indices for the qc wells to delete as well as plate 5
indices_to_delete=[117, 148, 157, 185, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200,
       201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213,
       214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226,
       227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
       240, 241, 242, 243, 244, 245, 246, 247]

def evaluate_inference(result, distance_name=None, hitlist_file='/media/lalil0u/New/workspace2/Xb_screen/data/mitocheck_exp_hitlist_perPheno.pkl', 
                       print_=False,folder='/media/lalil0u/New/projects/drug_screen/results/',
                       threshold=0.0001):
    f=open(hitlist_file, 'r')
    hitlist=pickle.load(f); f.close()
    yqualdict=expSi('../data/mapping_2014/qc_export.txt')
    dictSiEntrez=siEntrez('../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt')
    trad=EnsemblEntrezTrad('../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt')
    
    to_plot_info=defaultdict(dict)
    to_GO_info=defaultdict(list)
    
    for cond in sorted(result):
        if print_:
            print '-----------', cond
        l=np.array([(dictSiEntrez[el[0]], el[-1]) for el in result[cond]])
        to_plot_info[cond]=defaultdict(list)
        for el in sorted(KNOWN_TARGETS[cond.split('--')[0]]):
            if np.where(l[:,0]==el)[0].shape[0]>0:
                rank_=np.where(l[:,0]==el)[0][0]
                to_plot_info[cond]['genes'].append(rank_)
                to_plot_info[cond]['gene_list'].append(el)
                if print_:
                    print el, l[np.where(l[:,0]==el)], rank_
        
        if type(threshold)==int:
            gene_lim=l[:threshold][:,0]
        else:
            lim=np.where(np.array(l[:,1], dtype=float)==threshold)[0]
            gene_lim=l[lim][:,0]
        to_GO_info[cond]=[trad[el] for el in gene_lim]
        res=[]
        for gene in gene_lim:
            res.extend(hitlist[gene])
    #STEP2 here : we add writing gene list files
            
        to_plot_info[cond]['type']= Counter(res)
    
    multipleGeneListsToFile([to_GO_info[el] for el in to_GO_info], [el for el in to_GO_info], name=os.path.join(folder, 'GO_{}.txt'.format(distance_name)))
    
    return to_plot_info

def global_evaluate_inference(distance_name_list,folder='/media/lalil0u/New/projects/drug_screen/results/',
                              cmap=mpl.cm.OrRd,num_replicates=3,
                              threshold=0.0001):
    total_plot_info={}
    for distance in distance_name_list:
        print distance
        f=open(os.path.join(folder, 'inference_{}replicates/inference_{}.pkl'.format(num_replicates,distance)))
        result=pickle.load(f); f.close()
        
        total_plot_info[distance]=evaluate_inference(result, threshold=threshold, distance_name=distance,
                                            folder=os.path.join(folder, 'inference_{}replicates'.format(num_replicates)))
    
    r = drug_screen_utils.plotInferenceResult(distance_name_list,total_plot_info, 
                                                cmap=cmap)
    
    return r

def condition_clustering(distance_name, folder='/media/lalil0u/New/projects/drug_screen/results/', color_gradient='YlOrRd',
                         hit_only=False,
                          level=0.4):
    '''
    We can do drug clustering based on their distances to eachother but we can also do drug clustering based on their distances to Mitocheck
    THINK ABOUT THAT
    '''
    distances, who_, exposure_, _=_return_right_distance(distance_name, folder, check_internal= True)
    
    plates=[int(el.split('--')[0].split('_')[1]) for el in who_]
    exposure_wPL=np.array(['{}{:>10}'.format(exposure_[i], plates[i]) for i in range(len(exposure_))])

    f=open(os.path.join(folder, 'DS_hits_1.5IQR.pkl'))
    _, who_hits, drug_hits, dose_hits=pickle.load(f)
    f.close()
    if hit_only:
        wh_=np.hstack((np.where(who_==el)[0] for el in who_hits))
        distances=distances[wh_]; distances=distances[:,wh_]
        print wh_
        who_=who_[wh_]
    
    print distances.shape
    clusters=hierarchical_clustering.heatmap(distances, row_header=exposure_wPL, column_header=exposure_wPL, 
                                    row_method='ward', column_method='ward', 
                                    row_metric='euclidean', column_metric='euclidean', 
                                    color_gradient=color_gradient, filename='clustering_{}'.format(distance_name), 
                                    other_data=None, level=level,
                                    colorbar_ticklabels=['', '', ''],
                                    colorbar_title='Distance',
                                    range_normalization=(scoreatpercentile(distances.flatten(),10), scoreatpercentile(distances.flatten(), per=90)))

    exposure_hits=np.array(['{}--{}'.format(drug_hits[i], dose_hits[i]) for i in range(len(who_hits))])
    global_=np.bincount(clusters)
    hit_clusters = clusters[np.where(np.array([el in who_hits for el in who_]))]
    
    print global_, np.bincount(hit_clusters)
    who_cluster_hits={k: Counter(exposure_hits[np.where(hit_clusters==k)]) for k in range(1,np.max(clusters)+1)}
    
    return hit_clusters, who_cluster_hits



def check_distance_consistency(distance_name,
                               folder='/media/lalil0u/New/projects/drug_screen/results/', 
                               check_internal=True):
    '''
   Idea: check if from one plate to another, the experiments for the same exposure are close to eachother for this particular distance
   - M: distance matrix of size (hits, hits)
   - who: list of length hits
   - exposure: list of length hits with drug--dose (discrete and not actual concentration) 
'''
    distances, who_hits, exposure_hits, mito_who=_return_right_distance(distance_name, folder, check_internal)
    plates=[int(el.split('--')[0].split('_')[1]) for el in who_hits]
    exposure_hits_wPL=np.array(['{}{:>10}'.format(exposure_hits[i], plates[i]) for i in range(len(exposure_hits))])
    ind=np.argsort(plates)
    
    if check_internal:
        distances=distances[ind]; distances=distances[:,ind]
        who_hits=who_hits[ind]
        exposure_hits=exposure_hits[ind]
        exposure_hits_wPL=exposure_hits_wPL[ind]
        print who_hits
        
    exposure_hits2=[(el.split('--')[0], int(el.split('--')[1])) for el in exposure_hits]
    d=Counter(exposure_hits2)
    distinct_exposure=filter(lambda x:d[x]>=2, d)
    ind=np.hstack((np.where(np.array(exposure_hits)=='{}--{}'.format(*cond))[0] for cond in sorted(distinct_exposure, key=itemgetter(0,1))))
    if check_internal:
        M=distances[ind]; M=M[:,ind]
        print exposure_hits_wPL[ind]
        plotInternalConsistency(M, exposure_hits_wPL[ind])
        
    else:
        r=[]
        for cond in sorted(distinct_exposure, key=itemgetter(0,1)):
            print cond,
            currM=distances[np.where(np.array(exposure_hits)=='{}--{}'.format(*cond))]
            
            corr=np.mean(np.reshape([pearsonr(currM[i], currM[j])[0] for (i,j) in product(range(currM.shape[0]),repeat=2)],
                                    newshape=(currM.shape[0], currM.shape[0]))[np.triu_indices(currM.shape[0], 1)])
            print corr
            r.append(corr)
        return sorted(distinct_exposure, key=itemgetter(0,1)), np.array(r)
    
def measure_condition_separation(distance_name_list,
                                 folder='/media/lalil0u/New/projects/drug_screen/results/'):
    '''
    This function compares the different distances in distance_name_list, and more precisely their ability to distinguish
    between the different conditions of the drug screen, but not between different replicates of the conditions. Hence it computes
    the sum of distances(replicates)/distances(other conditions)
'''
    result={}
    
    for distance_name in distance_name_list:
        print distance_name
        distances, _, exposure_, _=_return_right_distance(distance_name, folder, check_internal=True)
        
        result[distance_name]= _compute_replicate_dist_vs_else(distances, exposure_)
    return result

def _compute_replicate_dist_vs_else(distances, exposure_list):
    distinct_exposure=sorted(list(set(exposure_list)))
    result=[]
    for exposure in distinct_exposure:
        replicates=np.where(exposure_list==exposure)[0]
        
        curr_distances=distances[replicates]
        
        others=filter(lambda x: x not in replicates, range(len(exposure_list)))
        
        numerateur = np.sum(curr_distances[:,replicates][np.triu_indices(replicates.shape[0], k=1)])
        denominateur = curr_distances[:,others]/float(replicates.shape[0])
        assert(denominateur.shape[0]==replicates.shape[0])
        assert(denominateur.shape[1]==len(others))
        if numerateur>0:
            result.append(numerateur/float(np.sum(denominateur)))
        
    return np.array(result)

def functional_inference(M, who_hits, exposure_hits, who_Mitocheck, num_permutations, threshold, random_result, taking_siRNAs=False, num_replicates=3):
    '''
    - M: distance matrix of size (hits, mitocheck)
    - taking_siRNAs: indicates if you want to consider different values of the same siRNAs independently (False) or as replicates of the same condition
'''
    r={}
    yqualdict=expSi('../data/mapping_2014/qc_export.txt')
    yqualdict.update(expSi('../data/mapping_2014/qc_validation_exp.txt', primary_screen=False))
    dictSiEntrez=siEntrez('../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt')
    
    siRNAs=[yqualdict[e] for e in who_Mitocheck]
    genes=[dictSiEntrez[e] for e in siRNAs]
    
    count=Counter(exposure_hits)
    #This way I look at exposures that are hits at least 50% of the times/plates
    for el in filter(lambda x: count[x]>=num_replicates, count):
        print el
        where_=np.where(exposure_hits==el)[0]
        
        batch_names = who_hits[where_]
        
        curr_dist=np.hstack((M[j] for j in where_))
        curr_who=[(batch_names[0], '') for k in range(M.shape[1])]
        
        curr_conditions=list(who_Mitocheck) if not taking_siRNAs else list(siRNAs)
        for name in batch_names[1:]:
            curr_who.extend([(name, '') for k in range(M.shape[1])])
            if not taking_siRNAs:
                curr_conditions.extend(who_Mitocheck)
            else:
                curr_conditions.extend(siRNAs)
            
        r[el], random_result=computeRPpvalues(curr_dist, np.array(curr_who), conditions=np.array(curr_conditions), technical_replicates_key=np.median, 
                     num_permutations=num_permutations, reverse=False, 
                     batch_names=batch_names, random_result=random_result,
                     signed=False)
        r[el]=sorted(r[el], key=itemgetter(2))
        pval=np.array(np.array(r[el])[:,-1], dtype=float)
        if not taking_siRNAs:
            currG=[dictSiEntrez[yqualdict[e]] for e in np.array(r[el])[np.where(pval<=threshold)][:,0]]
        else:
            currG=[dictSiEntrez[e] for e in np.array(r[el])[np.where(pval<=threshold)][:,0]]
        print sorted(Counter(currG).keys())
        
        
    return r

def _return_right_distance(distance_name, folder, check_internal):
    '''
    Possible distances:
    - transport: transport distance global on time
    - 'U_pheno_score': not normalized pheno scores
    - 'N_pheno_score': normalized pheno scores
    - 'ttransport_MAX': maximal transport distance on timepoints
    - 'ttransport_INT': sum of all transport distance on timepoints
    
    Check_internal:
    - True: then we want distances from DS exp to DS exp, including hits but not only => 904 datapoints
    - False: we're looking at correlation of distances to Mitocheck hits of DS hits, so only 184 experiments
    
    Little reminder: we have 904 movies in the drug screen that passed the QC and 806 experiments, 98 controls
    
    '''    
    mito_who=None
    if not check_internal:
        f=open(os.path.join(folder, 'DS_hits_1.5IQR.pkl'))
        res_, who_, drugs_, doses_=pickle.load(f)
        f.close()
        
        print 'Dealing with ', len(who_), 'experimental hits'
    else:
        f=open(os.path.join(folder, 'DS_pheno_scores.pkl'))
        res_, who_, _, drugs_, doses_=pickle.load(f)
        f.close()
        
        print '{} experiments from the drug screen'.format(len(who_))
    
    exposure_hits=['{}--{}'.format(drugs_[i], doses_[i]) for i in range(len(who_))]
    
    if distance_name == 'transport':
        #time aggregated transport distance. Here the drug screen is at the end
        f=open(os.path.join(folder, 'all_Mitocheck_DS_agg_transport_10.pkl'))
        distance, who=pickle.load(f); f.close()

        mito_positions=(806,distance.shape[0])

        if not check_internal:
            distances=np.vstack((distance[np.where(np.array(who)==el), mito_positions[0]:mito_positions[1]] for el in who_))[:,0]
            mito_who=who[mito_positions[0]:mito_positions[1]]#_mito_who_transform(who, *mito_positions)
        else:
            distances=np.vstack((distance[np.where(np.array(who)==el)] for el in who_))
            distances=np.hstack((distances[:, np.where(np.array(who)==el)] for el in who_))[:,:,0]
        
    elif 'ttransport' in distance_name:
        #non time aggregated transport distance
        if not check_internal:
            if distance_name=='ttransport_MAX':
                f=open(os.path.join(folder, 'all_Mitocheck_DS_UNagg_transport_10_MAX.pkl'))
                distance, who=pickle.load(f); f.close()
            elif distance_name=='ttransport_INT':
                f=open(os.path.join(folder, 'all_Mitocheck_DS_UNagg_transport_10_INT.pkl'))
                distance, who=pickle.load(f); f.close()
            else:
                raise ValueError(distance_name)
            mito_positions=(184,distance.shape[0])
            
            assert (len(who)==6904) 
            assert (len(who_)==184)

            distances=np.vstack((distance[np.where(np.array(who)==el)] for el in who_))
            distances=distances[:,mito_positions[0]:mito_positions[1]]
            mito_who=mito_who=who[mito_positions[0]:mito_positions[1]]#_mito_who_transform(who, *mito_positions)
        else:
            if distance_name=='ttransport_MAX':
                f=open(os.path.join(folder, 'DS_UNagg_transport_10_MAX.pkl'))
                distance=pickle.load(f); f.close()
            elif distance_name=='ttransport_INT':
                f=open(os.path.join(folder, 'DS_UNagg_transport_10_INT.pkl'))
                distance=pickle.load(f); f.close()
            else:
                raise ValueError(distance_name)
        #Here we have the controls as well
            f=open(os.path.join(folder, 'pheno_count_ALL_DS_time.pkl'))
            _,who=pickle.load(f); f.close()
            
            assert (len(who)==904) 
            assert (len(who_)==806)

            distances=np.vstack((distance[np.where(np.array(who)==el)] for el in who_))
            distances=np.hstack((distances[:, np.where(np.array(who)==el)] for el in who_))[:,:,0]
            
    elif 'pheno_score' in distance_name:
        f=open(os.path.join(folder, 'MITO_pheno_scores.pkl'))
        mito_scores,mito_who=pickle.load(f); f.close()
        
        scores=np.vstack((res_, mito_scores))
        scores=np.delete(scores, CLASSES.index('Anaphase'),1)
        scores=np.delete(scores, CLASSES.index('Interphase'),1)
#We need to normalize in all cases to be able to do meaningful distances. Otherwise some phenotypes have more importance than others

        if distance_name=='N_pheno_score':  
            scores=(scores-np.mean(scores,0))/np.std(scores,0)
            
        if not check_internal:
            distances = cdist(scores[:res_.shape[0]], scores[res_.shape[0]:], 'euclidean')
        else:
            distances = cdist(scores[:res_.shape[0]], scores[:res_.shape[0]], 'euclidean')
            
    elif distance_name=='nature':
        f=open(os.path.join(folder, 'all_Mitocheck_DS_nature_distance.pkl'))
        distance,who=pickle.load(f); f.close()
        mito_positions=(806,distance.shape[0])

        distances=np.vstack((distance[np.where(np.array(who)==el)] for el in who_))
        
        if check_internal:
            distances=np.hstack((distances[:, np.where(np.array(who)==el)] for el in who_))[:,:,0]
        else:
            distances=distances[:,mito_positions[0]:mito_positions[1]]
            mito_who=mito_who=who[mito_positions[0]:mito_positions[1]]#_mito_who_transform(who, *mito_positions)
    else:
        raise ValueError('Wrong distance name')
    print distances.shape
    return distances, np.array(who_), np.array(exposure_hits), np.array(mito_who)

def inference(distance_name, folder='/media/lalil0u/New/projects/drug_screen/results/', num_permutations=10000,
              taking_siRNAs=True, num_replicates=3,
               threshold=0.1, random_result=None):
    
    distances, who_hits, exposure_hits, mito_who=_return_right_distance(distance_name, folder, check_internal=False)
    
    r= functional_inference(distances, np.array(who_hits), np.array(exposure_hits), np.array(mito_who), 
                                num_permutations, threshold, random_result, taking_siRNAs, num_replicates=num_replicates)
        
    f=open(os.path.join(folder, 'inference_{}replicates/inference_{}.pkl'.format(num_replicates, distance_name)), 'w')
    pickle.dump(r, f); f.close()
    
    return r  
    

