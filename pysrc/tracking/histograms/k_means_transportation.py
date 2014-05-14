import warnings, pdb, time, sys, os

#sys.path.insert(0, '/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc/') #CHAUD
import cPickle as pickle
import numpy as np
from operator import itemgetter
import scipy.sparse as sp
from optparse import OptionParser
from sklearn.base import BaseEstimator, ClusterMixin, TransformerMixin
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster.k_means_ import _mini_batch_convergence, _tolerance
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.utils.sparsefuncs import mean_variance_axis0
from sklearn.utils import check_arrays
from sklearn.utils import check_random_state
from sklearn.utils import atleast2d_or_csr
from sklearn.utils import as_float_array
from sklearn.externals.joblib import Parallel
from sklearn.externals.joblib import delayed

from transportation import multSinkhorn,multEMD1d, costMatrix, ddcostMatrix, computingBins, computingddBins, findWassersteinBarycenter, findSinkhornBarycenter
from tracking.histograms import *
from util.sandbox import histLogTrsforming


def calculDistances(data, center, div_name, weight, lambda_=10, M=None, mat_hist_sizes=None, nb_feat_num=12, simulated=False):
    print 'divergence ', div_name

    print 'weight {}'.format(weight)
    
    distMat=_distances(data,center[np.newaxis,:], div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights=weight)
#_distances(X, centers,div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights)
    distances=np.sum(distMat**2)
    
    return distances


def k_means(X, n_clusters,lambda_, mat_hist_sizes, div_name = 'transportation',init='k-means++',
            n_jobs=1, M=None, nb_feat_num=12, dist_weights=None, previous_centers=None, n_init=10, 
            max_iter=300, verbose=False, tol=1e-4, random_state=None, copy_x=True):
    """K-means clustering algorithm.

    Parameters
    ----------
    X : array-like or sparse matrix, shape (n_samples, n_features)
        The observations to cluster.

    n_clusters : int
        The number of clusters to form as well as the number of
        centroids to generate.
        
        
    lambda_: constant to use for Sinkhorn version of transportation distance
    
    M: cost matrix for transportation distance
    
    mat_hist_sizes: matrix of shape (nb_hist,nb_features_hist) containing the number
        of bins in each case (e.g. nb of bins for the histogram of acceleration values
        after divisions)
        
    div_name: name of chosen divergence on histograms
        
    nb_feat_num: number of numerical features ie features that are just numbers as opposed
        to histograms.
        
    dist_weights: default, vector of ones. But one can choose to overweight some
        distances, eg between persistence histograms of normal displacements

    max_iter : int, optional, default 300
        Maximum number of iterations of the k-means algorithm to run.

    n_init : int, optional, default: 10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    init : {'k-means++', 'random' or 'incremental'}, optional
        Method for initialization, default to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': generate k centroids from a Gaussian with mean and
        variance estimated from the data.

        If an ndarray is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

        If a callable is passed, it should take arguments X, k and
        and a random state and return an initialization.

    tol : float, optional
        The relative increment in the results before declaring convergence.

    verbose : boolean, optional
        Verbosity mode.

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    copy_x : boolean, optional
        When pre-computing distances it is more numerically accurate to center
        the data first.  If copy_x is True, then the original data is not
        modified.  If False, the original data is modified, and put back before
        the function returns, but small numerical differences may be introduced
        by subtracting and then adding the data mean.

    n_jobs : int
        The number of jobs to use for the computation. This works by breaking
        down the pairwise matrix into n_jobs even slices and computing them in
        parallel.

        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debuging. For n_jobs below -1,
        (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one
        are used.

    Returns
    -------
    centroid : float ndarray with shape (k, n_features)
        Centroids found at the last iteration of k-means.

    label : integer ndarray with shape (n_samples,)
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    inertia : float
        The final value of the inertia criterion (sum of squared distances to
        the closest centroid for all observations in the training set).

    """
    random_state = check_random_state(random_state)
    assert div_name in DIVERGENCES, \
           'Given divergence {} is not implemented yet. Choose something in {}'.format(div_name, DIVERGENCES)
    best_inertia = np.infty
    if 'transportation' in div_name:
        M=_costMatrix(mat_hist_sizes, M)

    if dist_weights==None:
        dist_weights =np.ones(shape=(mat_hist_sizes.shape[1]+1))#np.reshape( np.ones(shape=(1+mat_hist_sizes.shape[0]*mat_hist_sizes.shape[1],1)),1+mat_hist_sizes.shape[0]*mat_hist_sizes.shape[1])
    print "CURRENT WEIGHTS", dist_weights
    
#    X = as_float_array(X, copy=copy_x)
    tol = _tolerance(X, tol)

    best_labels, best_inertia, best_centers = None, None, None
    if n_jobs == 1:
        # For a single thread, less memory is needed if we just store one set
        # of the best results (as opposed to one set per run per thread).
        for it in range(n_init):
            # run a k-means once
            labels, inertia, centers = _kmeans_single(
                X, n_clusters, init,
            #NOUVEAUX PARAMS
                div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,
                previous_centers,
                max_iter=max_iter, verbose=verbose,tol=tol, random_state=random_state)
            # determine if these results are the best so far
            if best_inertia is None or inertia < best_inertia:
                best_labels = labels.copy()
                best_centers = centers.copy()
                best_inertia = inertia
    else:
#BON
        # parallelisation of k-means runs
        seeds = random_state.randint(np.iinfo(np.int32).max, size=n_init)
        results = Parallel(n_jobs=n_jobs, verbose=0)(
            delayed(_kmeans_single)(X, n_clusters, init,
                                    div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,
                                    max_iter=max_iter,
                                    verbose=verbose, tol=tol,
                                    #x_squared_norms=x_squared_norms,
                                    # Change seed to ensure variety
                                    random_state=seed)
            for seed in seeds)
        # Get results with the lowest inertia
        labels, inertia, centers = zip(*results)
        best = np.argmin(inertia)
        best_labels = labels[best]
        best_inertia = inertia[best]
        best_centers = centers[best]

    return best_centers, best_labels, best_inertia

###############################################################################
# Initialization heuristic

def _incremental_init(X, n_clusters, previous_centers, div_name, lambda_, M, mat_hist_sizes, nb_feat_num, 
            dist_weights,previous_labels=None, n_local_trials=None, random_state=None):
    """Init n_clusters seeds according to k-means++

    Parameters
    -----------
    X: array or sparse matrix, shape (n_samples, n_features)
        The data to pick seeds for. To avoid memory copy, the input data
        should be double precision (dtype=np.float64).

    n_clusters: integer
        The number of seeds to choose
    previous_centers: the centers from previous initialization

    n_local_trials: integer, optional
        The number of seeding trials for each center (except the first),
        of which the one reducing inertia the most is greedily chosen.
        Set to None to make the number of trials depend logarithmically
        on the number of seeds (2+log(k)); this is the default.

    random_state: integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    Notes
    -----
    Selects initial cluster centers for k-mean clustering in a smart way
    to speed up convergence. see: Arthur, D. and Vassilvitskii, S.
    "k-means++: the advantages of careful seeding". ACM-SIAM symposium
    on Discrete algorithms. 2007

    Version ported from http://www.stanford.edu/~darthur/kMeansppTest.zip,
    which is the implementation used in the aforementioned paper.
    """
    
    if n_clusters==2:
        print 'k=2 so it is just a k-means++ init'
        return _k_init(X, n_clusters, div_name, lambda_, 
                       M, mat_hist_sizes, nb_feat_num, dist_weights, n_local_trials, random_state)
        
    n_previous_clusters = previous_centers.shape[0]
    n_clusters_to_find=n_clusters-n_previous_clusters+1

    random_state = check_random_state(random_state)
    
    if previous_labels==None:
    #We compute the label assignement to the previous centers
        previous_labels, _, _ = \
            _labels_inertia_precompute_dense(X, div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,previous_centers)

    biggest = np.argmax(np.bincount(previous_labels))
    new_centers = _k_init(X[np.where(previous_labels==biggest)], n_clusters_to_find, div_name, lambda_, 
            M, mat_hist_sizes, nb_feat_num, dist_weights, n_local_trials, random_state)
    
    centers = np.vstack((np.delete(previous_centers, biggest, 0), new_centers))
    
    return centers

def _k_init(X, n_clusters,div_name, lambda_, M, mat_hist_sizes, nb_feat_num, 
            dist_weights, n_local_trials=None, random_state=None):
    """Init n_clusters seeds according to k-means++

    Parameters
    -----------
    X: array or sparse matrix, shape (n_samples, n_features)
        The data to pick seeds for. To avoid memory copy, the input data
        should be double precision (dtype=np.float64).

    n_clusters: integer
        The number of seeds to choose

    n_local_trials: integer, optional
        The number of seeding trials for each center (except the first),
        of which the one reducing inertia the most is greedily chosen.
        Set to None to make the number of trials depend logarithmically
        on the number of seeds (2+log(k)); this is the default.

    random_state: integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    Notes
    -----
    Selects initial cluster centers for k-mean clustering in a smart way
    to speed up convergence. see: Arthur, D. and Vassilvitskii, S.
    "k-means++: the advantages of careful seeding". ACM-SIAM symposium
    on Discrete algorithms. 2007

    Version ported from http://www.stanford.edu/~darthur/kMeansppTest.zip,
    which is the implementation used in the aforementioned paper.
    """
    n_samples, n_features = X.shape
    random_state = check_random_state(random_state)

    centers = np.empty((n_clusters, n_features))

    # Set the number of local seeding trials if none is given
    if n_local_trials is None:
        # This is what Arthur/Vassilvitskii tried, but did not report
        # specific results for other than mentioning in the conclusion
        # that it helped.
        n_local_trials = 2 + int(np.log(n_clusters))

    # Pick first center randomly
    center_id = random_state.randint(n_samples)
    if sp.issparse(X):
        centers[0] = X[center_id].toarray()
    else:
        centers[0] = X[center_id]

    # Initialize list of closest distances and calculate current potential
#    if x_squared_norms is None:
#        x_squared_norms = _squared_norms(X)
    
    closest_dist_sq = _distances(X, centers[0][np.newaxis,:], div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights)**2
#ATTENTION THIS FUNCTION DOES NOT RETURN SQUARED DISTANCES
    current_pot = np.sum(closest_dist_sq)

    # Pick the remaining n_clusters-1 points
    for c in xrange(1, n_clusters):
        # Choose center candidates by sampling with probability proportional
        # to the squared distance to the closest existing center
        rand_vals = random_state.random_sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(closest_dist_sq.cumsum(), rand_vals)

        # Compute distances to center candidates
        distance_to_candidates = _distances(X,X[candidate_ids], div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights)

        # Decide which candidate is the best
        best_candidate = None
        best_pot = None
        best_dist_sq = None
        for trial in xrange(n_local_trials):
            # Compute potential when including center candidate
            new_dist_sq = np.minimum(closest_dist_sq,
                                     distance_to_candidates[trial])
            new_pot = new_dist_sq.sum()

            # Store result if it is the best local trial so far
            if (best_candidate is None) or (new_pot < best_pot):
                best_candidate = candidate_ids[trial]
                best_pot = new_pot
                best_dist_sq = new_dist_sq

        # Permanently add best center candidate found in local tries
        if sp.issparse(X):
            centers[c] = X[best_candidate].toarray()
        else:
            centers[c] = X[best_candidate]
        current_pot = best_pot
        closest_dist_sq = best_dist_sq

    return centers

def _distances(X, centers,div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights):
    #print 'X and centers shapes', X.shape, centers.shape
#get the centers matrix in the right shape
    centers=centers.T

    result=np.zeros(shape=(centers.shape[1], X.shape[0]))
    if 'transportation' in div_name and M is None:
        M={}
        for m in range(mat_hist_sizes.shape[1]):
            for n in range(mat_hist_sizes.shape[0]):
                M[(n,m)]=costMatrix(mat_hist_sizes[n,m], 1)

    #So here we take care of euclidean distance in numerical features space, except
    #if the corresponding weight is set to 0
    if dist_weights[0]!=0:
        for l in range(X.shape[0]):
            result[:,l]=dist_weights[0]*np.sqrt(np.sum((centers[:nb_feat_num]-X[l,:nb_feat_num, np.newaxis])**2, 0))
    #Here taking care of distances between histograms
    #mat_hist_sizes is typically array([[50, 50, 50, 50, 50]]) so the second shape is the number of histogram features
    up_to=nb_feat_num
    
    for m in range(mat_hist_sizes.shape[1]):         
        if dist_weights[1+m]==0:
            up_to+=mat_hist_sizes[0,m]
            continue
        for n in range(mat_hist_sizes.shape[0]):
            cur_hist_size=mat_hist_sizes[n,m]
            #we skip this histogram feature if the weight is 0
            if div_name=='transportation':#in fact here it's the Sinkhorn divergence. The exact transportation distance is coded with "etransportation"
                for l in range(X.shape[0]):
                    if np.sum(X[l, up_to:up_to+cur_hist_size])==0:
                        raise AttributeError
                    else:
                        result[:,l]+=dist_weights[1+m]*multSinkhorn(M[(n,m)], lambda_,X[l, up_to:up_to+cur_hist_size], centers[up_to:up_to+cur_hist_size])
            elif div_name == 'etransportation':
                for l in range(X.shape[0]):
                    if np.sum(X[l, up_to:up_to+cur_hist_size])==0:
                        raise AttributeError
                    else:
                        result[:,l]+=dist_weights[1+m]*multEMD1d(M[(n,m)], X[l, up_to:up_to+cur_hist_size], centers[up_to:up_to+cur_hist_size])
            elif div_name=='hellinger':
                for p in range(centers.shape[1]):
                    result[p]+=dist_weights[1+m]*np.sqrt(0.5*np.sum((np.sqrt(X[:,up_to:up_to+cur_hist_size])-np.sqrt(centers[up_to:up_to+cur_hist_size,p]))**2,1))
            elif div_name=='total_variation':
                for p in range(centers.shape[1]):
                    result[p]+=dist_weights[1+m]*np.sum(np.abs(X[:,up_to:up_to+cur_hist_size]-centers[up_to:up_to+cur_hist_size,p]),1)
                
        up_to+=cur_hist_size
    #print "DUREE", time.clock()-debut
    return result

def _kmeans_single(X, n_clusters,init, div_name,
                   lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,
                   previous_centers,
                   max_iter=300,
                   verbose=False, random_state=None,
                   tol=1e-4):
    """A single run of k-means, assumes preparation completed prior.

    Parameters
    ----------
    X: array-like of floats, shape (n_samples, n_features)
        The observations to cluster.

    k: int
        The number of clusters to form as well as the number of
        centroids to generate.
        
    lambda_: constant to use for Sinkhorn version of transportation distance
    
    M: cost matrix for transportation distance
    
    mat_hist_sizes: matrix of shape (nb_hist,nb_features_hist) containing the number
        of bins in each case (e.g. nb of bins for the histogram of acceleration values
        after divisions)
        
    nb_feat_num: number of numerical features ie features that are just numbers as opposed
        to histograms.
        
    dist_weights: default, vector of ones. But one can choose to overweight some
        distances, eg between persistence histograms of normal displacements

    max_iter: int, optional, default 300
        Maximum number of iterations of the k-means algorithm to run.

    init: {'k-means++', 'random', or ndarray, or a callable}, optional
        Method for initialization, default to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': generate k centroids from a Gaussian with mean and
        variance estimated from the data.

        If an ndarray is passed, it should be of shape (k, p) and gives
        the initial centers.

        If a callable is passed, it should take arguments X, k and
        and a random state and return an initialization.
        
        ONLY RANDOM INIT POUR L'INSTANT

    tol: float, optional
        The relative increment in the results before declaring convergence.

    verbose: boolean, optional
        Verbosity mode

    x_squared_norms: array, optional
        Precomputed x_squared_norms. Calculated if not given.

    random_state: integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    Returns
    -------
    centroid: float ndarray with shape (k, n_features)
        Centroids found at the last iteration of k-means.

    label: integer ndarray with shape (n_samples,)
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    inertia: float
        The final value of the inertia criterion (sum of squared distances to
        the closest centroid for all observations in the training set).
    """
    
    random_state = check_random_state(random_state)
    n_samples = X.shape[0]
    best_labels, best_inertia, best_centers = None, None, None
    # init

    centers = _init_centroids(X, n_clusters, init,
                div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights, 
                previous_centers=previous_centers,
                random_state=random_state)
    if verbose:
        print 'Initialization complete, type ', init

    # Allocate memory to store the distances for each sample to its
    # closer center for reallocation in case of ties PAS COMPRIS
    distances = np.zeros(shape=(X.shape[0],), dtype=np.float64)

    # iterations
    for i in range(max_iter):
        print "ITERATION ", i
        centers_old = centers.copy()
        # labels assignement is also called the E-step of EM
        labels, inertia, distances = \
            _labels_inertia_precompute_dense(X, div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,centers)

        # computation of the means is also called the M-step of EM
        centers = _centers_dense(X, labels, n_clusters, distances)

        if verbose:
            print 'Iteration %i, inertia %s' % (i, inertia)

        if best_inertia is None or inertia < best_inertia:
            best_labels = labels.copy()
            best_centers = centers.copy()
            best_inertia = inertia
        distances_old_new_centers = _distances(centers_old, centers, 
                    div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights)
        
        if np.sum(np.diagonal(distances_old_new_centers)** 2) < tol:
            if verbose:
                print 'Converged to similar centers at iteration', i
            break
    return best_labels, best_inertia, best_centers

def _init_centroids(X, k, init, div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,
                    random_state=None, init_size=None, previous_centers=None):
    """Compute the initial centroids

    Parameters
    ----------

    X: array, shape (n_samples, n_features)

    k: int
        number of centroids

    init: {'k-means++', 'random' or ndarray or callable} optional
        Method for initialization

    random_state: integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    init_size : int, optional
        Number of samples to randomly sample for speeding up the
        initialization (sometimes at the expense of accurracy): the
        only algorithm is initialized by running a batch KMeans on a
        random subset of the data. This needs to be larger than k.

    Returns
    -------
    centers: array, shape(k, n_features)
    """
    random_state = check_random_state(random_state)
    n_samples = X.shape[0]

    if init_size is not None and init_size < n_samples:
        if init_size < k:
            warnings.warn(
                "init_size=%d should be larger than k=%d. "
                "Setting it to 3*k" % (init_size, k),
                RuntimeWarning, stacklevel=2)
            init_size = 3 * k
        init_indices = random_state.permutation(n_samples)[:init_size]#random_integers(0, n_samples - 1, init_size)
        X = X[init_indices]

        n_samples = X.shape[0]
    elif n_samples < k:
            raise ValueError(
                "n_samples=%d should be larger than k=%d" % (n_samples, k))

    if init == 'k-means++':
        centers = _k_init(X, k, div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,
                        random_state=random_state)
    elif init == 'random':
        seeds = random_state.permutation(n_samples)[:k]
        centers = X[seeds]
    elif init=='incremental':
        if previous_centers==None and k>2:
            raise ValueError("Previous centers array should not be empty")
        else:
            centers = _incremental_init(X, k, previous_centers, div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,
                        random_state=random_state)
    elif hasattr(init, '__array__'):
        centers = init
    elif callable(init):
        centers = init(X, k, random_state=random_state)
    else:
        raise ValueError("the init parameter for the k-means should "
                         "be 'k-means++' or 'random' or an ndarray, "
                         "'%s' (type '%s') was passed." % (init, type(init)))

    if sp.issparse(centers):
        centers = centers.toarray()

    if len(centers) != k:
        raise ValueError('The shape of the inital centers (%s) '
                         'does not match the number of clusters %i'
                         % (centers.shape, k))

    return centers

#def _labels_inertia(X, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights,
#                    centers, distances=None):
#    """E step of the K-means EM algorithm
#
#    Compute the labels and the inertia of the given samples and centers
#
#    Parameters
#    ----------
#    X: float64 array-like or CSR sparse matrix, shape (n_samples, n_features)
#        The input samples to assign to the labels.
#
#    x_squared_norms: array, shape (n_samples,)
#        Precomputed squared euclidean norm of each data point, to speed up
#        computations.
#
#    centers: float64 array, shape ( n_clusters,n_features)
#        The cluster centers.
#
#    distances: float64 array, shape (n_samples,)
#        Distances for each sample to its closest center.
#
#    Returns
#    -------
#    labels: int array of shape(n)
#        The resulting assignment
#
#    inertia: float
#        The value of the inertia criterion with the assignment
#    """
##    n_samples = X.shape[0]
#    # set the default value of centers to -1 to be able to detect any anomaly easily
##    labels = - np.ones(n_samples, np.int32)
##    if distances is None:
##        distances = np.zeros(shape=(0,), dtype=np.float64)
##    if sp.issparse(X):
##        inertia = _k_means._assign_labels_csr(
##            centers, labels, distances=distances)
##    else:
#    
#    return _labels_inertia_precompute_dense(X, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights, centers)
##    inertia = _k_means._assign_labels_array(X, centers, labels, distances=distances)
##    
##    return labels, inertia


def _labels_inertia_precompute_dense(X, div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights, centers):
    n_samples = X.shape[0]
    k = centers.shape[0]
#OK SO HERE WE COMPUTE THE DISTANCES AS WE WANT
    distances = _distances(X, centers,div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights)

    labels = np.empty(n_samples, dtype=np.int32)
    labels.fill(-1)
    mindist = np.empty(n_samples)
    mindist.fill(np.infty)
    for center_id in range(k):
        dist = distances[center_id]
        labels[dist < mindist] = center_id
        mindist = np.minimum(dist, mindist)
    inertia = (mindist**2).sum()
    #print inertia
    return labels, inertia, mindist

def _centers_dense(X, labels, n_clusters, distances):
    """M step of the K-means EM algorithm

    Computation of cluster centers / means.
    
    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
    
    labels: array of integers, shape (n_samples)
    Current label assignment
    
    n_clusters: int
    Number of desired clusters
    
    distances: array-like, shape (n_samples)
    Distance to closest cluster for each sample.
    
    Returns
    -------
    centers: array, shape (n_clusters, n_features)
    The resulting centers
    """
    ## TODO: add support for CSR input
    n_samples = X.shape[0]
    n_features = X.shape[1]
    centers = np.zeros(shape=(n_clusters, n_features))
    n_samples_in_cluster = np.bincount(labels, minlength=n_clusters)
    empty_clusters = np.where(n_samples_in_cluster == 0)[0]
    # maybe also relocate small clusters?

    if len(empty_clusters):
        # find points to reassign empty clusters to
        far_from_centers = distances.argsort()[::-1]

    for i, cluster_id in enumerate(empty_clusters):
        # XXX two relocated clusters could be close to each other
        new_center = X[far_from_centers[i]]
        centers[cluster_id] = new_center
        n_samples_in_cluster[cluster_id] = 1

    for i in range(n_samples):
#        for j in range(n_features):
        centers[labels[i]] += X[i]

    centers /= n_samples_in_cluster[:, np.newaxis]

    return centers

def _costMatrix(mat_hist_sizes,*args):   
    '''
    Computes the ground metric or cost matrix for transportation distance. If the
    input args is of length 0 then default behaviour is cost matrix with bin number
    differences, power 1. If args[0] is an int then cost matrix with bin number
    differences power int. And if args[0] is an array then cost matrix with absolute bin values
    differences at power 1.
    
    ''' 
    if len(args)==0 or args[0]==None:
        M={}
        for m in range(mat_hist_sizes.shape[1]):
            for n in range(mat_hist_sizes.shape[0]):
                M[(n,m)]=costMatrix(mat_hist_sizes[n,m], power=1)
        return M
    elif type(args[0])==int:
        M={}
        for m in range(mat_hist_sizes.shape[1]):
            for n in range(mat_hist_sizes.shape[0]):
                M[(n,m)]=costMatrix(mat_hist_sizes[n,m], power=args[0])
        return M
    elif type(args[0])==np.ndarray:
        print "YOUPI"
        bins=list(args[0])
        M={}
        for m in range(mat_hist_sizes.shape[1]):
            for n in range(mat_hist_sizes.shape[0]):
                M[(n,m)]=costMatrix(mat_hist_sizes[n,m], bins=bins[m])
        return M
    elif type(args[0])==dict:
        return args[0]

    

def _mini_batch_step(X, div_name, lambda_, M, mat_hist_sizes,
                nb_feat_num, dist_weights, centers, counts,
                     old_center_buffer, compute_squared_diff,
                     distances=None, random_reassign=False,
                     random_state=None, reassignment_ratio=.01,
                     verbose=False):
    """Incremental update of the centers for the Minibatch K-Means algorithm

    Parameters
    ----------

    X: array, shape (n_samples, n_features)
        The original data array.

    x_squared_norms: array, shape (n_samples,)
        Squared euclidean norm of each data point.

    centers: array, shape (k, n_features)
        The cluster centers. This array is MODIFIED IN PLACE

    counts: array, shape (k,)
         The vector in which we keep track of the numbers of elements in a
         cluster. This array is MODIFIED IN PLACE

    distances: array, dtype float64, shape (n_samples), optional
        If not None, should be a pre-allocated array that will be used to store
        the distances of each sample to its closest center.

    random_state: integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    random_reassign: boolean, optional
        If True, centers with very low counts are
        randomly-reassigned to observations in dense areas.

    reassignment_ratio: float, optional
        Control the fraction of the maximum number of counts for a
        center to be reassigned. A higher value means that low count
        centers are more easily reassigned, which means that the
        model will take longer to converge, but should converge in a
        better clustering.

    verbose: bool, optional
        Controls the verbosity

    """
    # Perform label assignement to nearest centers
    nearest_center, inertia, _ = _labels_inertia_precompute_dense(X, div_name,
                lambda_, M, mat_hist_sizes,
                nb_feat_num, dist_weights, centers)
    if random_reassign and reassignment_ratio > 0:
#TODO
        random_state = check_random_state(random_state)
        # Reassign clusters that have very low counts
        to_reassign = np.logical_or(
            (counts <= 1), counts <= reassignment_ratio * counts.max())
        number_of_reassignments = to_reassign.sum()
        if number_of_reassignments:
            sys.stderr.write('Attention ici on est donc rentre dans la boucle de random_reassign')
            # Pick new clusters amongst observations with a probability
            # proportional to their closeness to their center
            distance_to_centers = np.asarray(centers[nearest_center] - X)
            distance_to_centers **= 2
            distance_to_centers = distance_to_centers.sum(axis=1)
            # Flip the ordering of the distances
            distance_to_centers -= distance_to_centers.max()
            distance_to_centers *= -1
            rand_vals = random_state.rand(number_of_reassignments)
            rand_vals *= distance_to_centers.sum()
            new_centers = np.searchsorted(distance_to_centers.cumsum(),
                                          rand_vals)
            new_centers = X[new_centers]
            if verbose:
                n_reassigns = to_reassign.sum()
                if n_reassigns:
                    print("[_mini_batch_step] Reassigning %i cluster centers."
                          % n_reassigns)
            if sp.issparse(new_centers) and not sp.issparse(centers):
                new_centers = new_centers.toarray()
            centers[to_reassign] = new_centers

    # implementation for the sparse CSR reprensation completely written in
    # cython
#WE'RE NOT THERE YET
#    if sp.issparse(X):
#        return inertia, _k_means._mini_batch_update_csr(
#            X, x_squared_norms, centers, counts, nearest_center,
#            old_center_buffer, compute_squared_diff)

    # dense variant in mostly numpy (not as memory efficient though)
    k = centers.shape[0]
    squared_diff = 0.0
    for center_idx in range(k):
        # find points from minibatch that are assigned to this center
        center_mask = nearest_center == center_idx
        count = center_mask.sum()

        if count > 0:
        #getting current points and weight of the present center
            points = X[center_mask]
            old_count = counts[center_idx]
            weights = list(np.ones(shape=(count,))); weights.append(old_count)
            
            if compute_squared_diff:
                old_center_buffer[:] = centers[center_idx]

        #i. calculating barycenter with euclidean distance for the numeric features if needed
            if dist_weights[0]!=0:
                centers[center_idx][:nb_feat_num] = _barycenter(centers[center_idx][:nb_feat_num],
                                                               points[:,:nb_feat_num], 
                                                               weights, 
                                                               distance='euclidean')        
        #ii. then for the histogram features THAT WE ARE CURRENTLY USING we need to calculate Wasserstein barycenter
            up_to=nb_feat_num
            for m in range(mat_hist_sizes.shape[1]):
             #we skip this histogram feature if the weight is 0
                if dist_weights[1+m]==0:
                    up_to+=mat_hist_sizes[0,m]
                    continue         
                for n in range(mat_hist_sizes.shape[0]):
                    cur_hist_size=mat_hist_sizes[n,m]
                    centers[center_idx][up_to:up_to+cur_hist_size] = _barycenter(centers[center_idx][up_to:up_to+cur_hist_size],
                                                               points[:,up_to:up_to+cur_hist_size], 
                                                               weights, 
                                                               distance=div_name, lambda_=lambda_, M=M[(n,m)]) 
                up_to+=cur_hist_size
        #FINALLY
        # update the count statistics for this center
            counts[center_idx] += count

            # update the squared diff if necessary NEED TO CHANGE THE DISTANCE IF WE ARE USING ANOTHER ONE
            if compute_squared_diff:
                c=centers[center_idx][np.newaxis,:]
                squared_diff += _distances(c, old_center_buffer[np.newaxis,:], div_name, lambda_, M, mat_hist_sizes, nb_feat_num, dist_weights) 
                #np.sum((centers[center_idx] - old_center_buffer) ** 2)
    return inertia, squared_diff

def _barycenter(center, points, weights, distance, lambda_=None, M=None, max_iter=2000):
    
    old_count = weights[-1]
    new_count = np.sum(weights[:-1])
    if distance=="euclidean":
#Barycenter is then classical mean of the points
        # inplace remove previous count scaling
        center *= old_count

        # inplace sum with new points members of this cluster
        center += np.sum(points, axis=0)

        # inplace rescale to compute mean of all points (old and new)
        center /= (old_count + new_count)
        
        return center
    
    elif distance == 'etransportation':
 #wasserstein barycenter
        p=np.vstack((points, center))
        return findWassersteinBarycenter(p, weights, M, max_iter=max_iter)
    
    elif distance =='transportation':
        p=np.vstack((points, center))
        return findSinkhornBarycenter(p, weights,lambda_, M, max_iter=max_iter)
    
    elif distance == 'total_variation':
        #mediane
        p=np.vstack((points, [center]*old_count)) if old_count>0 else points
        return np.median(p, 0)
    else:
        raise AttributeError("This distance does not exist")

class clusteringHistograms():
    
    def __init__(self,folder, datafile, datasize, bins_type, cost_type, bin_size, algo_type, n_clusters, init,
                 div_name, lambda_=10, M=None, dist_weights=None, batch_size=1000, init_size=2000,
                 n_init=5, verbose=0, ddim = 0, preprocessed=0, random_state=None,max_no_improvement=None):
        
        assert(div_name in DIVERGENCES)
        assert datasize>=init_size, "not enough data for doing mini batch k-means : datasize is smaller than init_size"
#        assert ((init=='incremental' and previous_centers_file is not None and n_clusters>2) or (init!='incremental')),\
#         "You can't ask for incremental initialization and not provide previous_centers"
            
        self.folder=folder
        self.datafile=datafile
        self.datasize=datasize
        self.bins_type=bins_type
        self.mat_hist_sizes=np.array([[bin_size for k in range(5)]]) if ddim==0 else np.array([[bin_size for k in range(2)]])
        self.bin_size=bin_size
        self.algo=algo_type
        self.n_clusters=n_clusters
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
        self.preprocessed=preprocessed
        self.numericFile = 'geigerNumericData.pkl'
        self.histogramFile = 'geigerHistogram_{}_bs{}_Data.pkl'.format(bins_type, bin_size)
    #plus tard : demander a Nell comment faire ce genre d'initialisation efficacement
        
    def _dataPrep(self):
        if not self.preprocessed:
            f=open(os.path.join(self.folder, self.datafile))
            raw_data = pickle.load(f); f.close()
            
            if len(raw_data)>6:
                #then it means the raw data comes from clustering.histConcatenation. Hence the log-trsformation has already been done
                r, histNC = raw_data[:2]    
            else:
                #then it means the raw data comes from trajPack.trajFeatures
                tabFeatures = raw_data[0]; histNC = raw_data[2]
                r = histLogTrsforming(tabFeatures)
                r=r[:,:-1]    
                
            if self.verbose:     
                print r.shape, 'not normalized'
                
            if self.verbose: print 'computing bins, bin type {}, d-dimensional? {}'.format(self.bins_type, bool(self.ddimensional))
            if self.ddimensional:
                histogrammeMatrix, bins = computingddBins(histNC, self.mat_hist_sizes[0], bin_type=self.bins_type)
            else:
                histogrammeMatrix, bins = computingBins(histNC, self.mat_hist_sizes[0], bin_type=self.bins_type)

        else:
            print 'loading processed data'
            f=open(os.path.join(self.folder, self.numericFile), 'r')
            r=pickle.load(f); f.close()
            f=open(os.path.join(self.folder,self.histogramFile), 'r')
            histogrammeMatrix, bins = pickle.load(f); f.close()
        
        r=(r-np.mean(r,0))/np.std(r,0)
        self.nb_feat_num = r.shape[1]  
        r=np.hstack((r, histogrammeMatrix))
        
        return r, bins
    
    def _seedsPrep(self, data, iteration, num_iterations_stability, fraction):
        self.random_state = check_random_state(self.random_state)
        
    #for k-means: if I do incremental initialization, I want to be on the same subset
    #in which case I don't need to recompute the total sum of squares. Let's do that for all init then
   

        filename = 'seeds_{}_{}_s{}_{}_w{}_{}_{}.pkl'.format(self.cost_type[:3],self.bins_type, self.bin_size,self.div_name[:5],self.dist_weights, self.datasize, iteration)
        if self.ddimensional:
            filename = 'dd_'+filename 
        if filename in os.listdir(self.folder):
            if self.verbose: 
                print "loading seeds and total sum of squares from file", os.path.join(self.folder,filename)
            
            f=open(os.path.join(self.folder, filename))
            set_, total_squared_dist= pickle.load(f); f.close()
            data=data[set_]
            if self.verbose: print "working with data set size ", data.shape
        else:
            #here we compute local_data_size, as indicated by the user through -s : on how many points does one want to work                
            local_data_size = min(self.datasize, data.shape[0])
            
            #if we're going to work on stability then we need only a fraction (80%) of the data
            if num_iterations_stability>0:
                stab_data_size=fraction*local_data_size
            else:
                stab_data_size=local_data_size
                
            set_ = self.random_state.permutation(data.shape[0])[:local_data_size]#random_integers(0, data.shape[0] - 1, min(self.datasize, data.shape[0]))
            data=data[set_]
            if self.verbose: print "working with data set size ", data.shape
            
            #computing inertia denominator WITH THE CORRECT DATA BARYCENTER, AND THIS DEPENDS ON THE DIVERGENCE THAT IS USED
            total_squared_dist = self._inertiaCalc(data[:stab_data_size])
            #saving
            print 'saving seeds and total sum of squares in file ',  os.path.join(self.folder,filename)
            f=open(os.path.join(self.folder, filename), 'w')
            pickle.dump([set_, total_squared_dist], f); f.close()
        return data, total_squared_dist

    def _inertiaCalc(self, data):
        #here we calculate the inertia of the total data set.
        data_center = np.ones(shape = (data.shape[1],))
        weights = np.ones(shape = (data.shape[0]+1,)); weights[-1]=0
        dist_weights = WEIGHTS[self.dist_weights] if not self.ddimensional else ddWEIGHTS[self.dist_weights]
#i. calculating barycenter with euclidean distance for the numeric features if needed
        if dist_weights[0]!=0:
            data_center[:self.nb_feat_num] = _barycenter(data_center[:self.nb_feat_num],
                                                               data[:,:self.nb_feat_num], 
                                                               weights, 
                                                               distance='euclidean')        
    #ii. then for the histogram features THAT WE ARE CURRENTLY USING we need to calculate data barycenter with the right distance/divergence
        up_to=self.nb_feat_num
        for m in range(self.mat_hist_sizes.shape[1]):
#we skip this histogram feature if the weight is 0
            if dist_weights[1+m]==0:
                up_to+=self.mat_hist_sizes[0,m]
                continue         
            for n in range(self.mat_hist_sizes.shape[0]):
                cur_hist_size=self.mat_hist_sizes[n,m]
                if self.verbose>5:
                    print 'up to', up_to, 'curr hist size ', cur_hist_size
                if 'transportation' not in self.div_name:
                    data_center[up_to:up_to+cur_hist_size] = _barycenter(data_center[up_to:up_to+cur_hist_size]/float(cur_hist_size),
                                                           data[:,up_to:up_to+cur_hist_size], 
                                                           weights, 
                                                           distance=self.div_name, lambda_=self.lambda_, M=self.cost_matrix[(n,m)])
                else:
                    data_center[up_to:up_to+cur_hist_size] = _barycenter(data_center[up_to:up_to+cur_hist_size]/float(cur_hist_size),
                                                           data[self.random_state.permutation(data.shape[0])[:500],up_to:up_to+cur_hist_size], 
                                                           weights[-501:], 
                                                           distance=self.div_name, lambda_=self.lambda_, M=self.cost_matrix[(n,m)]) 
            up_to+=cur_hist_size

        return calculDistances(data, data_center,self.div_name, dist_weights, self.lambda_, self.cost_matrix,\
                        self.mat_hist_sizes, nb_feat_num=self.nb_feat_num)
    
    def __call__(self, filename, iteration, only_dataprep, num_iterations_stability):
        filename+='_a{}_k{}_d{}_w{}_bt{}_bs{}_ct{}_{}.pkl'
        fraction = 0.8
        #Data preparation: putting histogram data in bins and normalizing
        
        #if self.cost_type==value then returning the bins for further ground metric calculation
        data, bins=self._dataPrep()
        
        if self.ddimensional:
            self.cost_matrix={(0,0):ddcostMatrix(np.product(self.mat_hist_sizes[0]), self.mat_hist_sizes[0])}
            self.mat_hist_sizes=np.array([[np.product(self.mat_hist_sizes)]])
            filename = 'dd_'+filename 
        else:   
            if self.cost_type=='value':
                print 'cost type value'
                self.cost_matrix=bins
            self.cost_matrix=_costMatrix(self.mat_hist_sizes, self.cost_matrix)
        
        #Data subset selection and total squared distances computation
        data, total_squared_dist=self._seedsPrep(data, iteration, num_iterations_stability, fraction)        
        if only_dataprep:
            return 1
        if self.algo==0:
            #K-means
#            if self.verbose:
#                print 'Launching k-means, k={}, divergence {}'.format(self.n_clusters, self.div_name)
#        #In case of incremental init I need to load the centers from previous k-means
#            if self.init=='incremental' and self.n_clusters>2:
#                previous_centers_file=filename.format(self.algo, self.n_clusters-1,self.div_name[:5], self.dist_weights, self.bins_type[:3], self.bin_size, self.cost_type[:3], iteration)
#                f=open(os.path.join(folder, previous_centers_file))
#                previous_centers, _,_=pickle.load(f); f.close()
#            else:
#                previous_centers=None            
#            
#            centers, labels, inertia=k_means(data, self.n_clusters,div_name=self.div_name,init=self.init,\
#                        lambda_=self.lambda_, mat_hist_sizes=self.mat_hist_sizes,\
#                        M=self.cost_matrix, nb_feat_num=self.nb_feat_num,\
#                        dist_weights=WEIGHTS[self.dist_weights], 
#                        previous_centers = previous_centers,n_init=self.n_init,
#                        verbose=self.verbose
#                        )
#            f=open(os.path.join(self.folder,\
#                    filename.format(self.algo, self.n_clusters,self.div_name[:5], self.dist_weights, self.bins_type[:3], self.bin_size,self.cost_type[:3], iteration)), 'w')
#            pickle.dump([centers, labels, inertia/float(total_squared_dist)], f)
#            f.close()
            
            return 1
        
        elif self.algo==1:
            #mini batch k-means
            if self.verbose:
                print 'Launching minibatch k-means, k={}, divergence {}, bin_type {}, bin_size {}, cost_size {}, dataset size {}, init size {}, batch size {}'.format(self.n_clusters, 
                                                    self.div_name, self.bins_type, self.bin_size, self.cost_type, self.datasize,\
                                                    self.init_size, self.batch_size)
    #In case of incremental init I need to load the centers from previous k-means
            if self.init=='incremental' and self.n_clusters>2:
                previous_centers_file=filename.format(self.algo, self.n_clusters-1,self.div_name[:5], self.dist_weights, self.bins_type[:3], self.bin_size, self.cost_type[:3], iteration)
                f=open(os.path.join(folder, previous_centers_file))
                previous_centers, _,_=pickle.load(f); f.close()
            else:
                previous_centers=None     
                
            dist_weights = WEIGHTS[self.dist_weights] if not self.ddimensional else ddWEIGHTS[self.dist_weights]
            if num_iterations_stability==0:        
                model = histogramMiniKMeans(self.n_clusters,self.lambda_, self.mat_hist_sizes,
                         div_name=self.div_name, M=self.cost_matrix, 
                         dist_weights=dist_weights, nb_feat_num = self.nb_feat_num,
                         init=self.init, previous_centers = previous_centers,
                         batch_size=self.batch_size,
                         init_size =self.init_size, verbose=self.verbose,n_init=self.n_init)            
                model.fit(data)
                centers, labels, inertia=model.cluster_centers_, model.labels_, model.inertia_
                
                f=open(os.path.join(self.folder, 
                        filename.format(self.algo, self.n_clusters,self.div_name[:5], self.dist_weights, self.bins_type[:3], self.bin_size,self.cost_type[:3], iteration)), 'w')
                pickle.dump([centers, labels, inertia/float(total_squared_dist)], f)
                f.close()
            elif num_iterations_stability>0:
#IMPLEMENTATION DU TEST POUR LA STABILITE DU CLUSTERING TEL QU'IL FIGURE DANS BEN-HUR ET AL 2002
                stability=[]; inertia=[]
                for it_ in range(num_iterations_stability):
                    if self.verbose>0:
                        print 'stability iteration ', it_
                    
                        print 'ONE'
                    model1 = histogramMiniKMeans(self.n_clusters,self.lambda_, self.mat_hist_sizes,
                         div_name=self.div_name, M=self.cost_matrix, 
                         dist_weights=dist_weights, nb_feat_num = self.nb_feat_num,
                         init=self.init, previous_centers = previous_centers,
                         batch_size=self.batch_size,
                         init_size =self.init_size, verbose=self.verbose,n_init=self.n_init)
                    set1= self.random_state.permutation(data.shape[0])[:int(fraction*data.shape[0])]; set1.sort()
                    local_data1=data[set1]
                    model1.fit(local_data1)
                    labels1=np.array(model1.labels_)
                    inertia.append(model1.inertia_)
                    
                    if self.verbose>0:print 'TWO'
                    
                    model2 = histogramMiniKMeans(self.n_clusters,self.lambda_, self.mat_hist_sizes,
                         div_name=self.div_name, M=self.cost_matrix, 
                         dist_weights=dist_weights, nb_feat_num = self.nb_feat_num,
                         init=self.init, previous_centers = previous_centers,
                         batch_size=self.batch_size,
                         init_size =self.init_size, verbose=self.verbose,n_init=self.n_init)
                    set2= self.random_state.permutation(data.shape[0])[:int(fraction*data.shape[0])]; set2.sort()
                    local_data2=data[set2]
                    model2.fit(local_data2)
                    labels2=np.array(model2.labels_)
                    inertia.append(model2.inertia_)
                    
                    if self.verbose>0:print 'fini'
                    intersection = filter(lambda x: x in set2, set1)
                    indices1=[np.where(set1==i)[0][0] for i in intersection]
                    indices2=[np.where(set2==i)[0][0] for i in intersection]
                    labels1=labels1[indices1]; labels2=labels2[indices2]
                    
                    corr1=np.array([labels1==label for label in labels1])
                    corr2=np.array([labels2==label for label in labels2])
                    for i in range(corr1.shape[0]):
                        corr1[i,i]=0; corr2[i,i]=0
                    
                    stability.append(np.sum(corr1*corr2)/np.sqrt(np.sum(corr1*corr1)*np.sum(corr2*corr2)))
                    if self.verbose>0:print 'stability calculation over'
                    
                    
                f=open(os.path.join(self.folder, 
                        filename.format(self.algo, self.n_clusters,self.div_name[:5], self.dist_weights, self.bins_type[:3], self.bin_size,self.cost_type[:3], "ST")), 'w')
                pickle.dump((stability, np.array(inertia)/float(total_squared_dist)), f)
                f.close()
            return 1
        raise ValueError("Self.algo takes only value 0 or 1")
            


class histogramMiniKMeans(MiniBatchKMeans):
    """
    Three different initializations:
    {random, k-means++ or incremental}
    implemented for three histogram divergences:
    {Hellinger, total variation and Sinkhorn distances}
    """
    
    def __init__(self, n_clusters,lambda_,mat_hist_sizes,
                 div_name='transportation',
                 M=None,dist_weights=None, 
                 nb_feat_num=12, init='k-means++', previous_centers=None, max_iter=300,
                 batch_size=500, verbose=0, compute_labels=True,
                 random_state=None, tol=1e-4, max_no_improvement=10,
                 init_size=1000, n_init=5, reassignment_ratio=0.01):
        
        super(histogramMiniKMeans, self).__init__(
            n_clusters=n_clusters, init=init, max_iter=max_iter, batch_size=batch_size,
            compute_labels=compute_labels, max_no_improvement=max_no_improvement,
            init_size=init_size,reassignment_ratio=reassignment_ratio,
            verbose=verbose, random_state=random_state, tol=tol, n_init=n_init)
        
        assert(div_name in DIVERGENCES)
        assert ((init=='incremental' and previous_centers is not None and n_clusters>2) or (init!='incremental')),\
         "You can't ask for incremental initialization and not provide previous_centers"
        
        if dist_weights==None:
            dist_weights =np.ones(shape=(mat_hist_sizes.shape[1]+1))#np.reshape( np.ones(shape=(1+mat_hist_sizes.shape[0]*mat_hist_sizes.shape[1],1)),1+mat_hist_sizes.shape[0]*mat_hist_sizes.shape[1])
        print "Init mini-batch k-means, CURRENT WEIGHTS", dist_weights
            
        self.previous_centers=previous_centers
        self.div_name=div_name
        self.cost_matrix=M
        self.lambda_=lambda_
        self.mat_hist_sizes=mat_hist_sizes
        self.dist_weights=dist_weights
        self.nb_feat_num=nb_feat_num

    def fit(self, X, y=None):
        """Compute the centroids on X by chunking it into mini-batches.

        Parameters
        ----------
        X: array-like, shape = [n_samples, n_features]
            Coordinates of the data points to cluster
        """
        random_state = check_random_state(self.random_state)
#SO THIS SHOULD NOT HAPPEN
        X = check_arrays(X, sparse_format="csr", copy=False,
                         check_ccontiguous=True, dtype=np.float64)[0]
        n_samples, n_features = X.shape
        if n_samples < self.n_clusters:
            raise ValueError("Number of samples smaller than number "
                             "of clusters.")

        if hasattr(self.init, '__array__'):
            self.init = np.ascontiguousarray(self.init, dtype=np.float64)

        if self.tol > 0.0:
            tol = _tolerance(X, self.tol)

            # using tol-based early stopping needs the allocation of a
            # dedicated before which can be expensive for high dim data:
            # hence we allocate it outside of the main loop
            old_center_buffer = np.zeros(n_features, np.double)
        else:
            tol = 0.0
            # no need for the center buffer if tol-based early stopping is
            # disabled
            old_center_buffer = np.zeros(0, np.double)

        distances = np.zeros(self.batch_size, dtype=np.float64)
        n_batches = int(np.ceil(float(n_samples) / self.batch_size))
        n_iter = int(self.max_iter * n_batches)

        init_size = self.init_size
        if init_size is None:
            init_size = 3 * self.batch_size
        if init_size > n_samples:
            init_size = n_samples
        self.init_size_ = init_size
        if self.verbose:
            print 'init_size', self.init_size_
            print 'batch_size', self.batch_size
            
        validation_indices = random_state.permutation(n_samples)[:init_size]#random_integers(0, n_samples - 1, init_size)
        X_valid = X[validation_indices]

        # perform several inits with random sub-sets
        best_inertia = None
        for init_idx in range(self.n_init):
            if self.verbose:
                print "Init %d/%d with method: %s" % (
                    init_idx + 1, self.n_init, self.init)
            counts = np.zeros(self.n_clusters, dtype=np.int32)

            # TODO: once the `k_means` function works with sparse input we
            # should refactor the following init to use it instead.

            # Initialize the centers using only a fraction of the data as we
            # expect n_samples to be very large when using MiniBatchKMeans
            cluster_centers = _init_centroids(X, self.n_clusters, self.init,
                self.div_name, self.lambda_, self.cost_matrix, self.mat_hist_sizes, 
                self.nb_feat_num, self.dist_weights,
                previous_centers=self.previous_centers,
                random_state=random_state,
                init_size=init_size)

#TO CHECK IS THAT WHAT WE WANT TO DO
            if self.verbose: print "Compute the label assignement on the init dataset"
            
#HERE we basically initialize the vector **counts** and update **cluster_centers** once, on a new init dataset
#necessary ??

            # Sample a minibatch from the full dataset
#            init_indices = random_state.permutation(n_samples)[:self.init_size_]#random_integers(0, n_samples - 1, self.init_size_)
#            #calculate the inertia on the init dataset
#            batch_inertia, centers_squared_diff = _mini_batch_step(
#                X[init_indices], self.div_name,
#                self.lambda_, self.cost_matrix, self.mat_hist_sizes, 
#                self.nb_feat_num, self.dist_weights,
#                cluster_centers, counts, old_center_buffer, compute_squared_diff=False,
#                verbose=self.verbose)#par defaut random_reassign=False
            
            # Keep only the best cluster centers across independent inits on
            # the common validation set

            _, inertia, _ = _labels_inertia_precompute_dense(X_valid,self.div_name, self.lambda_, self.cost_matrix, self.mat_hist_sizes, 
                self.nb_feat_num, self.dist_weights, cluster_centers)
            if self.verbose:
                print "Inertia for init %d/%d: %f" % (
                    init_idx + 1, self.n_init, inertia)
                
#TODO manquent l'initialisation de count et la mise a jour des cluster_centers
#batch_inertia n'est pas reutilise
#centers_squared_diff, old_center_buffer valent zero puirsque compute_squared_diff vaut False

            if best_inertia is None or inertia < best_inertia:
                self.cluster_centers_ = cluster_centers
                self.counts_ = counts
                best_inertia = inertia

        # Empty context to be used inplace by the convergence check routine
        convergence_context = {}

        # Perform the iterative optimization until the final convergence
        # criterion
        for iteration_idx in xrange(n_iter):
            # Sample a minibatch from the full dataset
            minibatch_indices = random_state.permutation(n_samples)[:self.batch_size]

            # Perform the actual update step on the minibatch data
            batch_inertia, centers_squared_diff = _mini_batch_step(
                X[minibatch_indices], self.div_name,
                self.lambda_, self.cost_matrix, self.mat_hist_sizes, 
                self.nb_feat_num, self.dist_weights,
                self.cluster_centers_, self.counts_,
                old_center_buffer, compute_squared_diff=(tol > 0.0), distances=distances,
                # Here we randomly choose whether to perform
                # random reassignment: the choice is done as a function
                # of the iteration index, and the minimum number of
                # counts, in order to force this reassignment to happen
                # every once in a while
                random_reassign=((iteration_idx + 1)
                                 % (10 + self.counts_.min()) == 0),
                random_state=random_state,
                reassignment_ratio=self.reassignment_ratio,
                verbose=self.verbose)

            # Monitor convergence and do early stopping if necessary
            if _mini_batch_convergence(
                    self, iteration_idx, n_iter, tol, n_samples,
                    centers_squared_diff, batch_inertia, convergence_context,
                    verbose=self.verbose):
                break

        if self.compute_labels:
            if self.verbose:
                print 'Computing label assignements and total inertia'
            self.labels_, self.inertia_, self.mindist = _labels_inertia_precompute_dense(
                X,self.div_name, self.lambda_, self.cost_matrix, self.mat_hist_sizes, self.nb_feat_num, self.dist_weights,
                    self.cluster_centers_)

        return self
    
    def find_representatives(self, N):
        '''
        This function returns N trajectories per cluster, the N that are the closest to the cluster center.
        More precisely, it returns a list of indices in the array of trajectories from the experiment at hand,
        of length N*n_clusters.
        '''
        representatives = []
        dist_data = sorted(np.vstack((range(len(self.labels_)), self.labels_, self.mindist)).T, key=itemgetter(2))
        
        count={num_cluster:0 for num_cluster in range(self.n_clusters)}
        for el in dist_data:
            if count[el[1]]<10:
                representatives.append(el[0])
                count[el[1]]+=1
        
        return np.array(sorted(representatives), dtype=int)
            
            

if __name__ == '__main__':
    
    parser = OptionParser(usage="usage: %prog [options]")
    
    parser.add_option('--only_dataprep', type=int, dest='only_dataprep', default=0)
    
    parser.add_option('--stability', type=int, dest='stability', default=0)
    parser.add_option('--preprocessed', type=int, dest='preprocessed', default=0)
    parser.add_option('--ddimensional', type=int, dest='ddim', default=0)
    
    parser.add_option('--sim', type=int, dest='simulated', default=0)
    parser.add_option('--iter', type=str, dest='iter', default=0)
    parser.add_option('-a', type=int, dest='algo', default=1)
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)
    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option("-k", type=int, dest="n_cluster")
    parser.add_option("-w", type=int, dest="weights", default=0)
    parser.add_option("-l",type=int, dest="lambda_", default=10)
    parser.add_option("-s", dest="size", type=int,default=1000)
    parser.add_option("--batch_size", dest="batch_size", type=int,default=500)
    parser.add_option("--init", dest="init", type=str,default='k-means++')
    parser.add_option("--n_init", dest="n_init", type=int,default=5)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    (options, args) = parser.parse_args()
    
    #donc les parametres pas initialises pour minibatch k-means : init_size = 5,000
    #pour les deux : max_iter = 300
    filename='stab1'
    datafile=None 
    if options.simulated:
        datafile = 'hist_tabFeatures_WO.pkl'
        folder = '../resultData/simulated_traj'
    else:
        if not options.preprocessed:
            datafile = 'MAR_histdata_PCAtrsf.pkl'
        folder='../resultData/sinkhornKMeans'
    if not options.preprocessed:
        print 'using datafile ', datafile, 'going for stability testing ', options.stability
#    sys.path.append('/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc/')
    if options.n_cluster<0:
        for k in range(4,12):
            print 'k ---------------- = ', k
            model = clusteringHistograms(folder, datafile, options.size, options.bins_type, options.cost_type, options.bin_size, options.algo, k, options.init,
                         div_name=options.div_name, lambda_=options.lambda_, dist_weights=options.weights, batch_size=options.batch_size,
                         n_init=options.n_init, verbose=options.verbose, ddim =options.ddim, preprocessed = options.preprocessed)
            
            model(filename, options.iter, options.only_dataprep, options.stability)
            if options.only_dataprep:
                break
    else:
        model = clusteringHistograms(folder, datafile, options.size, options.bins_type,options.cost_type, options.bin_size, options.algo, options.n_cluster, options.init,
                     div_name=options.div_name, lambda_=options.lambda_, dist_weights=options.weights, batch_size=options.batch_size,
                     n_init=options.n_init, verbose=options.verbose, ddim =options.ddim,preprocessed = options.preprocessed)
        
        model(filename, options.iter, options.only_dataprep, options.stability)
#    if options.intrand=='0':
#        fichier= '../histogramsNtotfirst.pkl'
#    elif options.intrand=='SIM':
#        fichier = '../resultData/simulated_traj/histogramsNtotSim.pkl'
#    else:
#        fichier = '../rand{}_histogramsNtotfirst.pkl'.format(options.intrand)
#    f=open(fichier)
#    #r, histNtot, histApStot, histAvStot, histAvMtot, who,ctrlStatus, length, genes, sirna=pickle.load(f)
#    X=pickle.load(f)
#    f.close()
#    nb_bins_list=[50,50,50,50,50]; mat_hist_sizes=np.reshape(np.array(nb_bins_list), (1,5))
#    if options.algo==0:
#        print "Doing k-means using transportation distance on file {}, k={}, nb of data points = {}".format(fichier,options.n_cluster, options.size)
#    
#        centers, labels, inertia=k_means(X[:int(options.size)], options.n_cluster,options.lambda_, mat_hist_sizes,n_jobs=1, M=None, nb_feat_num=12,\
#                                         dist_weights=WEIGHTS[options.weights], div_name=options.div_name)
#        
#        f=open('../resultData/first_{}means{}_f{}_w{}_l{}.pkl'.format( options.n_cluster,options.div_name[:5],options.intrand, options.weights, options.lambda_), 'w')
#        pickle.dump([centers, labels, inertia], f)
#        f.close()
#    else:
#        print "Doing mini batch k-means using {} on file {}, k={}, nb of data points = {}".format(options.div_name,fichier,options.n_cluster, options.size)
#        print 'The options are the following: init {}, max_iter {}, batch_size {}, init_size {}, nb of initializations {}'.format('random',
#                    100, options.batch_size, 1000, options.n_init)
#        if options.n_cluster>0:
#            model = histogramMiniKMeans(options.n_cluster,options.lambda_,mat_hist_sizes,
#                     div_name=options.div_name, 
#                     dist_weights=WEIGHTS[options.weights], batch_size=options.batch_size, n_init=options.n_init)
#            model.fit(X[:int(options.size)])
#            centers, labels, inertia=model.cluster_centers_, model.labels_, model.inertia_
#            
#            f=open('../first_batch{}means{}_f{}_w{}_l{}.pkl'.format(options.n_cluster, options.div_name,options.intrand, options.weights, options.lambda_), 'w')
#            pickle.dump([centers, labels, inertia], f)
#            f.close()
#        else:
#            for k in range(2,19):
#                print 'k= {}'.format(k)
#                model = histogramMiniKMeans(k,options.lambda_,mat_hist_sizes,
#                         div_name=options.div_name, 
#                         dist_weights=WEIGHTS[options.weights], batch_size=options.batch_size, n_init=options.n_init)
#                model.fit(X[:int(options.size)])
#                centers, labels, inertia=model.cluster_centers_, model.labels_, model.inertia_
#                
#                f=open('../first_batch{}means{}_f{}_w{}.pkl'.format(k, options.div_name[:5],options.intrand, options.weights), 'w')
#                pickle.dump([centers, labels, inertia], f)
#                f.close()


