from math import pow
outFile = 'result_{}.pkl'

methods = ('mini-batchk-means',
           'kernelk-means',
           #'c-means',
           'gaussianmixture',
           'isomapping',
           'spectralclustering')

extended_params = {'mini-batchk-means':['n_clusters', 'max_iter','batch_size', 'init_size', 'init', 'n_init'],
          'kernelk-means':['n_clusters', 'max_iter', 'kernel','gamma', 'degree', 'coef0'],
          'c-means':['n_clusters', 'max_iter', 'n_iter', 'fuzzyfier'],
          'gaussianmixture':['n_components', 'n_iter', 'covariance_type', 'n_init'],
        #WARNING this is not a clustering algorithm but a non-linear dimension reductionality
          'isomapping':['n_components', 'max_iter', 'n_neighbors'],
          'spectralclustering':['n_clusters', 'n_init', 'assign_labels','affinity','gamma','n_neighbors']
          }
 

           #possible values for init : k-means++, random or ndarray of initialization centers
params = {'mini-batchk-means':['n_clusters', 'max_iter','batch_size', 'init_size', 'init', 'n_init'],
          
          #possible values for kernel ['rbf', 'sigmoid', 'polynomial', 'poly', 'linear', 'cosine']
          'kernelk-means':['n_clusters', 'max_iter', {'kernel' : ['gamma', 'degree', 'coef0']}],
          
          'c-means':['n_clusters', 'max_iter', 'n_iter', 'fuzzyfier'],
          
          #possible values for covariance_type : full, spherical, diag, tied
          'gaussianmixture':['n_components', 'n_iter', 'covariance_type', 'n_init'],
          
        #WARNING this is not a clustering algorithm but a non-linear dimension reductionality
          'isomapping':['n_components', 'max_iter', 'n_neighbors', 'k_means_n_clusters'],
          
          #possible values for affinity : 'nearest_neighbors', 'precomputed' or any kernel, or a callable for homemade kernel
          #possible values for assign_labels : kmeans, discretize. n_init is then the number of initialization for kmeans
          'spectralclustering':['n_clusters', 'n_init', 'assign_labels',{'affinity': ['gamma','n_neighbors']}]
          #for home made spectral clustering
          #['n_clusters', 'n_neighbors', 'sigma']
          }

param_sets = {'mini-batchk-means': [[5,10,17,20,25,30], 1000, 1000, 3000, 'k-means++', 5],
              'gaussianmixture':[[5,10,17,20,25,30], 1000, ['diag', 'full', 'tied', 'spherical'], 5],
              'isomapping':[range(5, 50, 5), 1000, [2,5, 10, 15, 20, 25,50], [5,10,17,20,25,30]],
              'kernelk-means': [[5,10,17,20,25,30], 1000, 
                                 {'rbf':[[pow(10, k) for k in range(-5,1)], 0,0], 
                                  'linear':[0,0,0], 
               #                   'polynomial':[[pow(10, k) for k in range(-5,6)], range(1,6), [pow(10, k) for k in range(-1,2)]],
               #                   'sigmoid':[[pow(10, k) for k in range(-5,6)], 0, [pow(10, k) for k in range(-1,2)]], 
                                  'cosine':[0,0,0]}],
              'spectralclustering':[[5,10,17,20,25,30],5,['kmeans', 'discretize'], 
                                    {'rbf':[[pow(10, k) for k in range(-5,1)],0], 'nearest_neighbors':[0,[5, 10, 15, 20, 25,50]]}]     
              }

pca_param = [False, 10, 25, 50]
whiten_param = [True, False]

indices={'mini-batch k-means':['cohesion', 'silhouette score'],
         
         #il faut refaire le silhouette score dans le RKHS
         'kernel k-means':['cohesion', 'silhouette score'],
         
         'c-means':[]         
         
         }