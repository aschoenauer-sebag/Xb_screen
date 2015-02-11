import pdb
import numpy as np

from collections import Counter


from analyzer import plates
from _collections import defaultdict

def rank_product(data, who, conditions, technical_replicates_key):
    '''
    This function should compute the rank product of all conditions.
    Given that on some plates the same condition may have been plated twice or thrice,
    technical_replicates_key says how to deal with that.
    
    Mode: reverse=True ie small ranks will correspond to big values in data.
    This is correct if data contains statistics or distances. 

    Technical_replicates_key: should indicate which rank to give in this case: median,
    max, min. Min means that we take the smallest possible rank for this condition on the plate, very much too optimistic.
    Max means we take the biggest possible rank for this condition on the plate, conservative.
    '''
    all_conditions = Counter(conditions).keys()
    num_technical_replicates=defaultdict(list)
    ranks = np.zeros(shape=(len(all_conditions), len(plates)), dtype=float)
    
    for j,plate in enumerate(plates):
        print plate
    
        local_cond=conditions[np.where(who[:,0]==plate)][np.argsort((-data)[np.where(who[:,0]==plate)])]
        
        for i,condition in enumerate(all_conditions):
            ranks[i,j]=technical_replicates_key(np.where(local_cond==condition)[0]+1) if condition in local_cond else -1

        num_technical_replicates[plate]=[len(np.where((conditions==condition)&(who[:,0]==plate))[0]) for condition in all_conditions]
       
        ranks[np.where(ranks[:,j]<0),j]=len(np.where(who[:,0]==plate)[0])
        ranks[:,j]/=len(np.where(who[:,0]==plate)[0])
        
    return all_conditions,np.prod(ranks, axis=1, dtype=float),num_technical_replicates

def computeRPpvalues(data, who, conditions, technical_replicates_key, num_permutations):
    all_conditions,real_result, num_technical_replicates = rank_product(data, who, conditions, technical_replicates_key)
    
    rrp=randomRankProduct(num_permutations)
    random_result = rrp(num_technical_replicates, technical_replicates_key)
    
    pvals=np.zeros(shape=len(all_conditions,), dtype=float)
    for i in range(real_result.shape[0]):
        pvals[i]=max(1, len(np.where(random_result[:,i]<real_result[i])[0]))/float(num_permutations)
        
    return zip(all_conditions,real_result, pvals)

class randomRankProduct(object):
    def __init__(self, num_permutations=1000):
        self.N = num_permutations
        
    def __call__(self,num_technical_replicates, technical_replicates_key):
        num_conditions = len(num_technical_replicates[plates[0]])
        result = np.zeros(shape=(self.N,num_conditions), dtype=float)
        
        for k in range(self.N):
            print k,
            result[k]=self.randomRankProduct(num_conditions,num_technical_replicates, technical_replicates_key)
            
        return result
            

    def randomRankProduct(self,num_conditions,num_technical_replicates, technical_replicates_key):
        '''
        This function simulates rank products for random orders, as close as possible to our experimental setting.
        '''
        ranks = np.zeros(shape=(num_conditions, len(plates)), dtype=float)
        
        for j,plate in enumerate(plates):
            curr_num=np.sum(num_technical_replicates[plate])
            
            local_cond=np.repeat(range(num_conditions), repeats=num_technical_replicates[plate])
            local_cond=np.random.permutation(local_cond)
            
            for condition in range(num_conditions):
                ranks[condition,j]=technical_replicates_key(np.where(local_cond==condition)[0]+1) if condition in local_cond else -1
            
            ranks[np.where(ranks[:,j]<0),j]=curr_num
            ranks[:,j]/=curr_num
            
        return np.prod(ranks, axis=1, dtype=float)
        
        
        
        