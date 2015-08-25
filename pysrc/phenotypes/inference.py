import os, pdb
from util.listFileManagement import expSi, siEntrez
from collections import Counter
import numpy as np
import cPickle as pickle
from analyzer.rank_product import computeRPpvalues

from drug_screen_utils import CLASSES, plotInternalConsistency
from scipy.spatial.distance import cdist
from scipy.stats.mstats_basic import scoreatpercentile
from phenotypes.drug_screen_utils import lim_Mito
from operator import itemgetter
from scipy.stats.stats import pearsonr
from itertools import product

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
    plates=[int(el.split('--')[0].split('_')[1]) for el in who_hits]; ind=np.argsort(plates)
    
    if check_internal:
        distances=distances[ind]; distances=distances[:,ind]
        who_hits=who_hits[ind]
        exposure_hits=exposure_hits[ind]
        print who_hits
        
    exposure_hits2=[(el.split('--')[0], int(el.split('--')[1])) for el in exposure_hits]
    d=Counter(exposure_hits2)
    distinct_exposure=filter(lambda x:d[x]>2, d)
    ind=np.hstack((np.where(np.array(exposure_hits)=='{}--{}'.format(*cond))[0] for cond in sorted(distinct_exposure, key=itemgetter(0,1))))
    if check_internal:
        M=distances[ind]; M=M[:,ind]
        print exposure_hits[ind], who_hits[ind]
        plotInternalConsistency(M, exposure_hits[ind])
        
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
        

def functional_inference(M, who_hits, exposure_hits, who_Mitocheck, num_permutations, threshold, random_result, taking_siRNAs=False):
    '''
    - M: distance matrix of size (hits, mitocheck)
    - taking_siRNAs: indicates if you want to consider different values of the same siRNAs independently (False) or as replicates of the same condition
'''
    r={}
    yqualdict=expSi('../data/mapping_2014/qc_export.txt')
    dictSiEntrez=siEntrez('../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt')
    
    siRNAs=[yqualdict[e] for e in who_Mitocheck]
    genes=[dictSiEntrez[e] for e in siRNAs]
    
    count=Counter(exposure_hits)
    #This way I look at exposures that are hits at least 50% of the times/plates
    for el in filter(lambda x: count[x]>2, count):
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
        print Counter(currG)
        
    return r

def _return_right_distance(distance_name, folder, check_internal):
    f=open(os.path.join(folder, 'DS_hits_1.5IQR.pkl'))
    res_hits, who_hits, drugs_hits, doses_hits=pickle.load(f)
    f.close()
    print 'Dealing with ', len(who_hits), 'experimental hits'
    exposure_hits=['{}--{}'.format(drugs_hits[i], doses_hits[i]) for i in range(len(who_hits))]
    
    if distance_name == 'transport':
        #time aggregated transport distance. Here the drug screen is at the end
        f=open(os.path.join(folder, 'all_Mitocheck_DS_agg_transport_10.pkl'))
        distance=pickle.load(f); f.close()
        f=open(os.path.join(folder, 'all_Mitocheck_DS_phenohit.pkl'))
        _,who=pickle.load(f); f.close()
        if not check_internal:
            distances=np.vstack((distance[np.where(np.array(who)==el), :lim_Mito] for el in who_hits))[:,0]
        else:
            distances=np.vstack((distance[np.where(np.array(who)==el)] for el in who_hits))
            
            distances=np.hstack((distances[:, np.where(np.array(who)==el)] for el in who_hits))[:,:,0]
        mito_who=['{}--{:>03}'.format(el.split('--')[0], int(el.split('--')[1])) for el in who[:lim_Mito]]
        
    elif distance_name=='ttransport':
        #non time aggregated transport distance
        OLD_SIZE=248
        f=open(os.path.join(folder, 'all_Mitocheck_DS_UNagg_transport_10_MAX.pkl'))
        distance=pickle.load(f); f.close()
        print distance.shape
        f=open(os.path.join(folder, 'all_Mitocheck_DS_phenohit_perFrame.pkl'))
        _,who=pickle.load(f); f.close()
        
        distances=np.vstack((distance[np.where(np.array(who)==el)] for el in who_hits))
        
        if not check_internal:
            distances=distances[:,OLD_SIZE:]
        else:
            distances=np.hstack((distances[:, np.where(np.array(who)==el)] for el in who_hits))[:,:,0]
        print distances.shape
            
        mito_who=['{}--{:>03}'.format(el.split('--')[0], int(el.split('--')[1])) for el in who[OLD_SIZE:]]
    
    elif distance_name=='pheno_score':
        f=open(os.path.join(folder, 'MITO_pheno_scores.pkl'))
        mito_scores,mito_who=pickle.load(f); f.close()
        
        scores=np.vstack((res_hits, mito_scores))
        scores=np.delete(scores, CLASSES.index('Anaphase'),1)
        scores=np.delete(scores, CLASSES.index('Interphase'),1)
#We need to normalize in all cases to be able to do meaningful distances. Otherwise some phenotypes have more importance than others
        scores=(scores-np.mean(scores,0))/np.std(scores,0)
        if not check_internal:
            distances = cdist(scores[:res_hits.shape[0]], scores[res_hits.shape[0]:], 'euclidean')
        else:
            distances = cdist(scores[:res_hits.shape[0]], scores[:res_hits.shape[0]], 'euclidean')
        
    return distances, np.array(who_hits), np.array(exposure_hits), np.array(mito_who)

def inference(distance_name, folder='/media/lalil0u/New/projects/drug_screen/results/', num_permutations=1000,
              taking_siRNAs=False,
               threshold=0.1, random_result=None):
    
    distances, who_hits, exposure_hits, mito_who=_return_right_distance(distance_name, folder, check_internal=False)
    
    return functional_inference(distances, np.array(who_hits), np.array(exposure_hits), np.array(mito_who), 
                                num_permutations, threshold, random_result, taking_siRNAs)
        
        
    

