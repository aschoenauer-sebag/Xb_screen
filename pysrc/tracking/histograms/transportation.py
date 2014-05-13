import pdb, time

import cplex as c
import cPickle as pickle
import numpy as np

from collections import Counter
from itertools import product
from math import fabs
from matplotlib import pyplot as p
from numpy import diag, dot
from numpy.linalg import inv, norm
from sklearn.utils import check_random_state
from scipy.stats import pearsonr, scoreatpercentile
from optparse import OptionParser
from sklearn.externals.joblib import Parallel
from sklearn.externals.joblib import delayed

from tracking.trajPack import featuresHisto
from tracking.plots import plotCompaDiv

#def testWassersteinBarycenter(points, barycenter, M):
#    distances = multEMD1d(M, barycenter, points.T)
#    print distances
#    return distances

def plotOptWasserstein(iter_val, filename='iter_val', title=None):
    import matplotlib.pyplot as p
    from tracking.trajPack.clustering import couleurs
    if type(iter_val)==list:
        fig = p.figure()
        ax=fig.add_subplot(111)
        for i, el in enumerate(iter_val):
                f=open('{}{}.pkl'.format(filename, el), 'r')
                gaps, bornes, norm_grad, obj= pickle.load(f)
                f.close()
                range_ = len(gaps)
                ax.plot(range(range_), gaps, color = '#'+couleurs[i], label = '{}'.format( el))
        p.grid(True)
        ax.set_yscale('log')
        p.legend()
        p.suptitle(title)
        p.show()
    
    elif type(iter_val)==int:        
        f=open('iter_val{}.pkl'.format(iter_val), 'r')
        gaps, bornes, norm_grad, obj= pickle.load(f)
        f.close()
        range_ = len(gaps)
        print min(gaps)
        f=p.figure()
        ax=f.add_subplot(221)
        ax.plot(range(range_), gaps)
        ax.set_title('Relative gaps between a_k and a(k-1)')
        ax=f.add_subplot(222)
        ax.plot(range(range_), bornes)
        ax.set_title("Borne theorique sur f(a_k) - f(x*)")
        ax=f.add_subplot(223)
        ax.plot(range(range_), norm_grad)
        ax.set_title("Norme L2 de grad f(a_k)")
    #    ax=f.add_subplot(224)
    #    ax.plot(range(range_), obj)
    #    ax.set_title("f(a_k)")
        p.suptitle(title)
        p.show()
    return
    
def findSinkhornBarycenter(points, weights,lamb, M, max_iter = 1000, eps=0.0001, iter_val = 0):
    #implementation beck and teboulle 2003
    
    #points.shape it should be num_points x nb of bins
    assert len(weights)==points.shape[0], "Not giving the weights to calculate barycenter"
    assert(M is not None), 'the cost matrix should be initialized'
    
    nb_points = points.shape[0]
    nb_bins = points[0].shape[0]
    gaps=[]; bornes=[]; norm_grad=[]; obj=[]
    debut = time.clock()

    Minf = np.max(M)
    t_0=np.sqrt(2*np.log(nb_bins))/float(Minf)

    print nb_points
    a=np.ones(shape=nb_bins)/nb_bins; a_old=None

    for k in range(1, max_iter):
        beta=t_0*1.0/np.sqrt(k)
        grad_ak=np.zeros(shape=nb_bins)
        alphas = multSinkhorn(M, lamb, a, points.T, returnAlphas=True)
        for num in range(nb_points):
            grad_ak=grad_ak + weights[num]*alphas[num] if weights[num]!=0 else grad_ak
            
        grad_ak/=np.sum(weights)
        a_old=np.array(a)
        a=a*np.exp(-t_0*beta*grad_ak); a/=np.sum(a)
        
        if k%10==0 and k>0:
            gap = norm(np.abs(a-a_old))/float(norm(a))
            #borne = np.sqrt(2*np.log(nb_bins)/float(k))*np.max(grad_ak)
            #gaps.append(gap); bornes.append(borne); norm_grad.append(norm(grad_ak));# obj.append(curr_value)
    #        print 'iteration', k, 'gap entre a_(k-1) et a_k', gap, 'borne sur f(x_k) - f(x*)', borne, "norme du gradient ", norm(grad_ak)
    
            #here it's useless to use the theoretical bound cause max(grad_ak)=+inf
#            if borne<eps:
#                print "borne theorique ok"
#                break
            if gap<eps:
                print "delta a ok"
                break
    print "DUREE", time.clock()- debut
    
#    f=open('iter_val{}.pkl'.format(iter_val), "w")
#    pickle.dump([gaps, bornes, norm_grad, obj], f)
#    f.close()
    return a

def findWassersteinBarycenter(points, weights, M=None, max_iter = 1000, eps=0.0001, filename='iter_val0'):
    #points.shape it should be num_points x nb of bins
    assert len(weights)==points.shape[0], "Not giving the weights to calculate barycenter"
    
    nb_points = points.shape[0]
    nb_bins = points[0].shape[0]
    gaps=[]; bornes=[]; norm_grad=[]; obj=[]
    debut = time.clock()
    if M==None:
        M=costMatrix(nb_bins, power=1)
    Minf = np.max(M)
    t_0=np.sqrt(2*np.log(nb_bins))/float(Minf)

    print nb_points
    a=np.ones(shape=nb_bins)/nb_bins; a_old=None

    for k in range(1, max_iter):
        beta=t_0*1.0/np.sqrt(k)
        grad_ak=np.zeros(shape=nb_bins)
        for num in range(nb_points):
            grad_ak=grad_ak + weights[num]*findDualOptimum(a, points[num], M) if weights[num]!=0 else grad_ak
            
        grad_ak/=np.sum(weights)
        a_old=np.array(a)
        a=a*np.exp(-t_0*beta*grad_ak); a/=np.sum(a)
        if k%10==0 and k>0:
#            curr_value = np.sum(multEMD1d(M, a, points.T)); curr_value/=float(nb_points)
            gap = norm(np.abs(a-a_old))/float(norm(a))
            borne = np.sqrt(2*np.log(nb_bins)/float(k))*np.max(grad_ak)
            gaps.append(gap); bornes.append(borne); #norm_grad.append(norm(grad_ak)); obj.append(curr_value)
    #        print 'iteration', k, 'gap entre a_(k-1) et a_k', gap, 'borne sur f(x_k) - f(x*)', borne, "norme du gradient ", norm(grad_ak)
            if borne<eps:
                print "borne theorique ok"
                break
            elif gap<eps:
                print "delta a ok"
                break
    print "DUREE", time.clock()- debut
#    f=open('{}.pkl'.format(filename), "w")
#    pickle.dump([gaps, bornes, norm_grad, obj], f)
#    f.close()
    return a

def findWassersteinIsobarycenter(points, M=None, max_iter = 1500, eps=0.001, filename='iter_val0'):
    #points.shape it should be num_points x nb of bins
    nb_points = points.shape[0]
    nb_bins = points[0].shape[0]
    #gaps=[]; bornes=[]; norm_grad=[]; obj=[]
    debut = time.clock()
    if M==None:
        M=costMatrix(nb_bins, power=1)
    Minf = np.max(M)
    t_0=np.sqrt(2*np.log(nb_bins))/float(Minf)

    print M, nb_bins, nb_points
    a=np.ones(shape=nb_bins)/nb_bins; a_old=None
    print "init ", a
    for k in range(1, max_iter):
        beta=t_0*1.0/np.sqrt(k)
        grad_ak=np.zeros(shape=nb_bins)
        for num in range(nb_points):
            grad_ak+=findDualOptimum(a, points[num], M)
            
        grad_ak/=nb_points
        a_old=np.array(a)
        a=a*np.exp(-t_0*beta*grad_ak); a/=np.sum(a)
        if k%10==0 and k>0:
#            curr_value = np.sum(multEMD1d(M, a, points.T)); curr_value/=float(nb_points)
#            gap = norm(np.abs(a-a_old))/float(norm(a))
            borne = np.sqrt(2*np.log(nb_bins)/float(k))*np.max(grad_ak)
            #gaps.append(gap); bornes.append(borne); norm_grad.append(norm(grad_ak)); obj.append(curr_value)
    #        print 'iteration', k, 'gap entre a_(k-1) et a_k', gap, 'borne sur f(x_k) - f(x*)', borne, "norme du gradient ", norm(grad_ak)
            if borne<eps:
                break
    print "DUREE", time.clock()- debut
#    f=open('{}.pkl'.format(filename), "w")
#    pickle.dump([gaps, bornes, norm_grad, obj], f)
#    f.close()
    return a

def findDualOptimum(a,b,M=None):
    if M==None:
        M=costMatrix(a.shape[0], power=1)
        
    objective=np.hstack((a,b))
    alpha=solve_Dual(objective, M)    
    return np.array(alpha)

def solve_Dual(objective, M):
    my_prob = c.Cplex()
    my_prob.set_log_stream(None); my_prob.set_results_stream(None)
    my_prob.objective.set_sense(my_prob.objective.sense.maximize)
    #print my_prob.get_problem_type()# bien LP
   #je passe la fonction objectif
    my_prob.variables.add(obj = objective,
#ATTENTION NE PAS DECLARER LES TYPES CAR SINON ON AURA UN MIP MEME SI TOUTES LES VARIABLES SONT CONTINUES
                          #types=[my_prob.variables.type.continuous for k in range(objective.shape[0])],
                          lb=[-c.infinity for k in range(objective.shape[0])]
                          )

    #je passe les contraintes
    for i,j in product(range(M.shape[0]),repeat=2):
        const=[[i,M.shape[0]+j], [1,1]]
        my_prob.linear_constraints.add(lin_expr=[const],
                                    senses=['L'], 
                                    rhs=[M[(i,j)]]
                                    )

    my_prob.solve()

    x = my_prob.solution.get_values()
#NB: MEME SI ON A QUE DES ENTIERS CA VEUT PAS DIRE QUE C'EST FAUX SI LA MATRICE M NE CONTIENT ELLE MEME QUE DES ENTIERS...
    return x[:M.shape[0]]

def plotCompaDivergences(nameList, nom):
    totalVar=[]; hellinger=[]; kullback=[]
    for name in nameList:
        f=open(name)
        data=pickle.load(f); f.close()
        totalVar.append(data[0])
        hellinger.append(data[1])
        kullback.append(data[2])
    totalVar.append('Total Variation vs Sinkhorn')
    hellinger.append('Hellinger divergence vs Sinhorn')
    kullback.append("Kullback-Leibler divergence vs Sinkhorn")
    print nameList
    plotCompaDiv(totalVar, hellinger, kullback, nom)
    
    return

    

def computingComparisons(X, name):
    ll=[(1,1),(1,2), (5, 1),(10,1)]
    
    for lamb_sink, power in ll:
        print 'lambda {}, power {}'.format(lamb_sink, power)
        compaHistogramDivergences(name, X, lamb_sink, None, 100, power)
        
    return

def compaHistogramDivergences(name, X, lamb_sink, M, nb_centers=100, power=1, eps=0.00001):
    """
    Comparing Sinkhorn distances (param lamb_sink and matrix cost M) with chi-squared divergences
    and total variation
    X: array of shape (n_points, n_bins), n_points not bigger than 1,000 cause we are computing distances
    """
    random_state = check_random_state(None)
    seeds = random_state.permutation(X.shape[0])[:nb_centers]
    centers = X[seeds]
    
    if M is None:
        M=costMatrix(X.shape[1], power)
    print 'Sinkhorn distances'
    sinkDist = np.zeros(shape=( X.shape[0],centers.shape[0]))
    for k in filter(lambda x: x not in seeds, range(X.shape[0])):
        sinkDist[k] = multSinkhorn(M, lamb_sink, X[k], centers.T)
        
    print 'Total variation distances'
    totalVarDiv = np.zeros(shape=( X.shape[0],centers.shape[0]))
    for p in range(nb_centers):
        totalVarDiv[:,p]=np.sum(np.abs(X-centers[p]),1)
    
    print 'Hellinger divergences'
    hellingerDiv = np.zeros(shape=( X.shape[0],centers.shape[0]))
    for p in range(nb_centers):
        hellingerDiv[:,p]=np.sqrt(0.5*np.sum((np.sqrt(X)-np.sqrt(centers[p]))**2,1))
#ET KL DIVERGENCE, ET CHI-SQUARED DIVERGENCE MISSING
    
    kullbackDiv =  np.zeros(shape=( X.shape[0],centers.shape[0]))
    for p in range(nb_centers):
        kullbackDiv[:,p]=-np.sum(centers[p]*np.log((X+eps)/(centers[p]+eps)), 1)

    res1=[pearsonr(sinkDist[:,r], totalVarDiv[:,r])[0] for r in range(nb_centers)]
    res2=[pearsonr(sinkDist[:,r], hellingerDiv[:,r])[0] for r in range(nb_centers)]
    res3=[pearsonr(sinkDist[:,r], kullbackDiv[:,r])[0] for r in range(nb_centers)]
    
    f=open('../corr_histo'+name+'_lamb{}_power{}.pkl'.format(lamb_sink, power), 'w')
    pickle.dump([res1, res2, res3], f); f.close()

    return 1


def _sinkhornNumDistance(X, lamb_sinkhorn, M, mat_hist_sizes, nb_feat_num, centers, dist_weights=None, power=1):
    """K-means clustering algorithm. Avoid using it, prefer the matrix implementation multSinkhornNumDistance
    Parameters
    ----------
    X : array-like or sparse matrix, shape (n_samples, n_features)
        The observations to cluster.
 
    lamb_sinkhorn: constant to use for Sinkhorn version of transportation distance
    
    M: cost matrix for transportation distance
        A VOIR 
    
    mat_hist_sizes: matrix of shape (nb_hist,nb_features_hist) containing the number
        of bins in each case (e.g. nb of bins for the histogram of acceleration values
        after divisions)
        
    nb_feat_num: number of numerical features ie features that are just numbers as opposed
        to histograms.
        
    dist_weights: default, vector of ones. But one can choose to overweight some
        distances, eg between persistence histograms of normal displacements
        
    centers: float64 array, shape (n_clusters, n_features)
        The cluster centers.
    """
    debut = time.clock()
    result=np.empty(shape=(centers.shape[0], X.shape[0]))

    if M is None:
        M={}
        for m in range(mat_hist_sizes.shape[1]):
            for n in range(mat_hist_sizes.shape[0]):
                M[(n,m)]=costMatrix(mat_hist_sizes[n,m], power)
    if dist_weights==None:
        dist_weights=np.ones(shape=(mat_hist_sizes.shape[1]+1))
    

    for l in range(X.shape[0]):
        #print 'traj ',l
        for p in range(centers.shape[0]):
         #   print p
            num=np.linalg.norm(centers[p,:nb_feat_num]-X[l,:nb_feat_num])**2
            up_to=nb_feat_num
            dist=np.empty(shape=mat_hist_sizes.shape)

            for m in range(mat_hist_sizes.shape[1]):
          #      print "feature ",featuresHisto[m]
                for n in range(mat_hist_sizes.shape[0]):
                    #print "histo ",n
                    cur_hist_size=mat_hist_sizes[n,m]
                    #try:
#                    if featuresHisto[m]=="straight":
#                        pdb.set_trace()
                    dist[n,m]=_sinkhorn(M[(n,m)], lamb_sinkhorn, centers[p, up_to:up_to+cur_hist_size], X[l, up_to:up_to+cur_hist_size])
#                    except np.linalg.LinAlgError:
#                        print 'PROBLEM LinAlgError'
                    up_to+=cur_hist_size
            result[p,l]=dist_weights[0]*num+np.dot(dist_weights[1:], np.reshape(dist.flatten(), dist_weights[1:].shape))
    print "DUREE", time.clock()- debut
    return result

def _sinkhorn(M, lamb, r,c, eps=0.01):
    if np.testing.assert_approx_equal(np.sum(r), 1, 0.001) or np.testing.assert_approx_equal(np.sum(c), 1, 0.001):
        print np.sum(r), np.sum(c), 'Given vectors are not in the probability simplex'
        raise AttributeError
    if np.all(r==c):
        return 0
    
    K=np.exp(-lamb*M)
#methode sioux a investiguer plus tard
#    ind=np.where(r>0)
#    r=r[ind]; M=M[ind]
#    K=np.exp(-lamb*M)
    u_old=np.reshape(np.ones(shape=(M.shape[0],1)), M.shape[0])  #bien la bonne initialisation suivant Sinkhorn et Knopp, 1967
    # the algorithm is linearly convergent whenever the original matrix K has total support
    c=np.array(c).reshape(M.shape[0])
    r=np.array(r).reshape(M.shape[0])
    
#    goOn=True
#    while goOn==True:
    for k in range(100000):
        u=dot(inv(diag(dot(K, dot(inv(diag(dot(K.T, u_old))), c)))), r)

        if k>0 and k%100==0:
    #        print k,
            v=dot(inv(diag(dot(K.T,u))),c)
            Plamb = dot(dot(diag(u), K), diag(v))
            zou = norm(r-dot(Plamb, np.ones(shape=(Plamb.shape[0],))), 1)
            zou2=norm(c-dot(Plamb.T,np.ones(shape=(Plamb.shape[0],))), 1)
            #print zou,zou2
            if max(zou,zou2)<eps:
                break
        u_old=u
            
    v=dot(inv(diag(dot(K.T,u))),c)
    Plamb = dot(dot(diag(u), K), diag(v))
    result =np.trace(dot(M.T, Plamb))
    #print 'result ', result 
    return result

def multSinkhornNumDistance(X, lamb_sinkhorn, M, mat_hist_sizes, nb_feat_num, centers, dist_weights=None, power=1):
    """K-means clustering algorithm, more efficient since it computes the distance of a point to all centers in the same computation 
    (matrix form of the sinkhornNumDistance).

    Parameters
    ----------
    X : array-like or sparse matrix, shape (n_samples, n_features)
        The observations to cluster.
 
    lamb_sinkhorn: constant to use for Sinkhorn version of transportation distance
    
    M: cost matrix for transportation distance
        A VOIR 
    
    mat_hist_sizes: matrix of shape (nb_hist,nb_features_hist) containing the number
        of bins in each case (e.g. nb of bins for the histogram of acceleration values
        after divisions)
        
    nb_feat_num: number of numerical features ie features that are just numbers as opposed
        to histograms.
        
    dist_weights: default, vector of ones. But one can choose to overweight some
        distances, eg between persistence histograms of normal displacements
        
    centers: float64 array, shape (n_clusters, n_features)
        The cluster centers.
    """
    print 'Warning, this version is not efficient'
    debut=time.clock()
#get the centers matrix in the right shape
    centers=centers.T

    result=np.empty(shape=(centers.shape[1], X.shape[0]))

    if M is None:
        M={}
        for m in range(mat_hist_sizes.shape[1]):
            for n in range(mat_hist_sizes.shape[0]):
                M[(n,m)]=costMatrix(mat_hist_sizes[n,m], power)
    if dist_weights==None:
        dist_weights=np.ones(shape=(mat_hist_sizes.shape[1]+1))
    num_dist=np.zeros(shape=(centers.shape[1], X.shape[0]))
    for l in range(X.shape[0]):
  #      print 'traj ',l
        num=np.sum((centers[:nb_feat_num]-X[l,:nb_feat_num, np.newaxis])**2, 0)#np.reshape(X[l,:nb_feat_num], (nb_feat_num, 1)))**2, 0)
        num_dist[:,l]=num
        up_to=nb_feat_num
        dist=np.empty(shape=(centers.shape[1], mat_hist_sizes.shape[0], mat_hist_sizes.shape[1]))

        for m in range(mat_hist_sizes.shape[1]):
 #           print "feature ",featuresHisto[m]
            for n in range(mat_hist_sizes.shape[0]):
#                print "histo ",n
                cur_hist_size=mat_hist_sizes[n,m]
                if np.sum(X[l, up_to:up_to+cur_hist_size])==0:
                    raise AttributeError
                else:
                    dist[:,n,m]=multSinkhorn(M[(n,m)], lamb_sinkhorn,X[l, up_to:up_to+cur_hist_size], centers[up_to:up_to+cur_hist_size])

                up_to+=cur_hist_size
        result[:,l]=dist_weights[0]*np.sqrt(num)+np.dot(np.reshape(dist.flatten(), (centers.shape[1], mat_hist_sizes.shape[1])), dist_weights[1:])
    print num_dist
    print "DUREE", time.clock()-debut
    return result


def ddcostMatrix(d, nb_bins_list, power=1, bins=None):
    if bins==None:
        
        M=np.array([[(i/nb_bins_list[0] - j/nb_bins_list[0])**2+(i%nb_bins_list[0] - j%nb_bins_list[0])**2 for i in range(d)] for j in range(d)])

        return np.sqrt(M)
    else:
        print 'error'#'computing cost matrix using bin value differences'
#        M=np.zeros(shape=(d,d))
#        for i in range(d):
#            for j in range(i+1,d):
#                M[i,j]=fabs((bins[i+1]-bins[i])/2 - (bins[j+1]-bins[j])/2)
#                M[j,i]=M[i,j]
        return None

def costMatrix(d, power=1, bins=None):
    if bins==None:
        M=np.array([k**power for k in range(d)])#np.zeros(shape=(d,d))
        for p in range(1,d):
            M=np.vstack((M, np.array([np.abs(k-p)**power for k in range(d)])))
        return M
    else:
        print 'computing cost matrix using bin value differences'
        M=np.zeros(shape=(d,d))
        for i in range(d):
            for j in range(i+1,d):
                M[i,j]=fabs((bins[i+1]-bins[i])/2 - (bins[j+1]-bins[j])/2)
                M[j,i]=M[i,j]
        return M

def multSinkhorn(M, lamb, r,C, eps=0.001, returnAlphas = False):
    """
    ATTENTION THE CENTER MATRIX SHOULD BE OF SHAPE N_FEATURES, N_CENTERS
    Function exactly as published in Cuturi 2013 NIPS, ie it includes saving one Schur product by keeping K_tild in memory. Youpi
    """

    if np.testing.assert_approx_equal(np.sum(r), 1, 0.001) or \
            np.testing.assert_allclose(np.sum(C, 0), np.ones(shape=(C.shape[1],)), 0, 0.001, "Center vectors are not in the probability simplex"):
        print np.sum(r), np.sum(C)
        raise AttributeError

    dist_zeros=np.where(np.sum(np.abs(C.T-r), 1)==0)[0]
    # the algorithm is linearly convergent whenever the original matrix K has total support
    
    r=np.array(r).reshape(M.shape[0])    
    
#methode sioux pour sauver un Schur product on va sauver quelque part la multiplication r/K
    ind=np.where(r>0); ind2=np.where(r<=0)
    r=r[ind]; M=M[ind[0], :]
    K=np.exp(-lamb*M)
    U_old=np.ones(shape=(len(ind[0]), C.shape[1]))  #bien la bonne initialisation suivant Sinkhorn et Knopp, 1967
    K_tild = dot(diag(1./r), K)
    
    for k in range(100000):
        U=1./dot(K_tild, C*1./dot(K.T, U_old) )

        if k>0 and k%100==0:
            #print k,
            V=C*1./dot(K.T, U)
            zou=0; zou2=0
            #NEED: more efficient implementation of the test
            for p in range(C.shape[1]):
                Plamb = dot(dot(diag(U[:,p]), K), diag(V[:,p]))
                zou = max(zou, norm(r-dot(Plamb, np.ones(shape=(Plamb.shape[1],))), 1))
                zou2=max(zou2, norm(C[:,p]-dot(Plamb.T,np.ones(shape=(Plamb.shape[0],))), 1))
                
            if max(zou,zou2)<eps:
                break
        U_old=U
            
#    V=C*1./dot(K.T, U)
    if not returnAlphas:
        result =np.sum(U*dot((K*M), C*1./dot(K.T, U)), 0)
    
        result[dist_zeros]=0
        return result
    else: #it means we're computing the gradient of the barycenter objective function
#U should be of size ind[0] x nb_points then we should insert zeros where r was equal to zero
        if len(ind2[0])>0:  
            U=np.insert(U, ind2, 0, 0) 
        alphas= (-1.0/lamb)*np.log(U)
#Here we deal with the fact that for any dual solution (alpha, beta) and t in R, (alpha -t*np.ones(), beta +t*np.ones()) is also a dual solution
#With or without does not change the value of the barycenter in the end
        #alphas-=alphas[0]        
        return alphas.T

def find_(begin, r):
    try:
        return np.where(r>begin)[0][0]
    except IndexError:
        return -1

def EMD1d(r,c,M):
    if np.all(r==c):
        return 0
    
    r=np.cumsum(r)
    c=np.cumsum(c)

    intervalles=Counter(r).keys()
    intervalles.extend(Counter(c).keys())
    intervalles.append(0)
    intervalles=Counter(intervalles).keys()
    intervalles.sort()
    result=0
    for i in range(len(intervalles)-1):
        begin = intervalles[i]
        end = intervalles[i+1]
        
        k=find_(begin, r)
        l=find_(begin, c)
        result+=(end-begin)*M[(k,l)]
    return result

def multEMD1d(M, r,C):
    """
    ATTENTION THE CENTER MATRIX IS OF SHAPE N_FEATURES, N_CENTERS
    """
    print M
    if np.testing.assert_approx_equal(np.sum(r), 1, 0.001) or \
            np.testing.assert_allclose(np.sum(C, 0), np.ones(shape=(C.shape[1],)), 0, 0.001, "Center vectors are not in the probability simplex"):
        print np.sum(r), np.sum(C)
        raise AttributeError
    
    r=np.array(r).reshape(M.shape[0])    
    result = np.zeros(shape=(C.shape[1]))
    for k in range(C.shape[1]):
        result[k]=EMD1d(r, C[:,k], M)
        
    return result

def computingddBins(histogramme, nb_bins_list, bin_type='minmax'):
   
    #Bin list to return in case needed to compute the ground matrix (if bin_type=='quantile')
    quantiles=None
    featuresHisto=["displacements", "straight"]
    if bin_type=='minmax':
        minMax=np.empty(shape=(len(featuresHisto), 2))
#remembering the trajectories that have infinite values, hence deletion
        toDel=[]

        for i, feature in enumerate(featuresHisto):
            lCour=[]; print feature
            for traj_num, ll in enumerate(histogramme[feature]):
                lCour.extend(np.array(ll)[np.where(np.isinf(np.array(ll))==False)])
                
                if len(np.where(np.isinf(np.array(ll))==False)[0])==0:
                    toDel.append(traj_num)
                    
            lCour=np.array(lCour)
            minMax[i,0]=np.min(lCour)
            minMax[i,1]=np.max(lCour)
        print "to delete ", len(toDel)
        
        nb_traj = len(histogramme[feature])
        matrice_result=np.empty(shape=(nb_traj-len(toDel), np.prod(nb_bins_list)), dtype=np.float64)
        matrice_result.fill(-1)
        
        for traj_num in filter(lambda x: x not in toDel, range(nb_traj)):
    #I put this size: nb of points as the minimum accross all features. But it might diminish if there are infinite values
            nbPoints = min([len(histogramme[feature][traj_num]) for feature in featuresHisto])
            currArray=np.zeros(shape=(nbPoints,len(featuresHisto)))
    
            for i,feature in enumerate(featuresHisto):
                currArray[:,i]=np.array(histogramme[feature][traj_num][:nbPoints])
                
            if np.any(np.isinf(currArray)):
#TO TEST
                pdb.set_trace()
                np.delete(currArray, np.where(np.isinf(currArray))[0][0], 0)
                
            zz=np.histogramdd(currArray, bins = nb_bins_list,range = [(minMax[i,0], minMax[i,1]) for i in range(len(featuresHisto))])[0]
        
            matrice_result[traj_num]=zz.flatten()/float(np.sum(zz))
                
    elif bin_type=='quantile':
        quantiles = np.empty(shape=(len(featuresHisto),), dtype=object)
        for i, feature in enumerate(featuresHisto):
            lCour=[]; print feature
            for ll in histogramme[feature]:
                lCour.extend(np.array(ll)[np.where(np.isinf(np.array(ll))==False)])
            lCour=np.array(lCour)
            quantiles[i] = [scoreatpercentile(lCour, per) for per in [k*100/float(nb_bins_list[i]) for k in range(nb_bins_list[i]+1)]]
        #so we should have nb_bins_list[i]+1 values in quantiles[i] since we also give the rightmost edge
        
        nb_traj = len(histogramme[feature])
        matrice_result=np.empty(shape=(nb_traj, np.sum(nb_bins_list)), dtype=np.float64)
        matrice_result.fill(-1)
        #remembering the trajectories that does not have this feature
        toAverage={feature:[] for feature in featuresHisto}
        for traj_num in range(nb_traj):
            up_to=0
            for i,feature in enumerate(featuresHisto):
                ll=np.array(histogramme[feature][traj_num])[np.where(np.isinf(histogramme[feature][traj_num])==False)]
                if len(ll)==0:
                    toAverage[feature].append(traj_num)
                else:
                    zz=np.histogram(ll, bins=quantiles[i])[0]
                    assert np.sum(zz)==len(ll)
                    matrice_result[traj_num, up_to:up_to+nb_bins_list[i]]=zz/float(len(ll))
                up_to+=nb_bins_list[i]
    else:
        raise ValueError('Bin_type variable cannot take {} as argument'.format(bin_type))

                                                                               
    return matrice_result, quantiles  

def computingBins(histogramme, nb_bins_list, bin_type='minmax', previous_binning=None):
   
    assert len(nb_bins_list)==len(featuresHisto)
    #Bin list to return in case needed to compute the ground matrix (if bin_type=='quantile')
    quantiles=None
    
    if bin_type=='minmax':
        if previous_binning is None:
            minMax=np.empty(shape=(len(featuresHisto), 2))
            for i, feature in enumerate(featuresHisto):
                lCour=[]; print feature
                for ll in histogramme[feature]:
                    lCour.extend(np.array(ll)[np.where(np.isinf(np.array(ll))==False)])
                lCour=np.array(lCour)
                minMax[i,0]=np.min(lCour)
                minMax[i,1]=np.max(lCour)
        else:
            minMax = previous_binning
            
        nb_traj = len(histogramme[feature])
        matrice_result=np.empty(shape=(nb_traj, np.sum(nb_bins_list)), dtype=np.float64)
        matrice_result.fill(-1)
        #remembering the trajectories that does not have this feature
        toAverage={feature:[] for feature in featuresHisto}
        for traj_num in range(nb_traj):
            up_to=0
            for i,feature in enumerate(featuresHisto):
                ll=np.array(histogramme[feature][traj_num])[np.where(np.isinf(histogramme[feature][traj_num])==False)]
                if len(ll)==0:
                    toAverage[feature].append(traj_num)
                else:
                    zz=np.histogram(ll, nb_bins_list[i],(minMax[i,0], minMax[i,1]))[0]
                    assert np.sum(zz)==len(ll)
                    matrice_result[traj_num, up_to:up_to+nb_bins_list[i]]=zz/float(len(ll))
                
                up_to+=nb_bins_list[i]
                
    elif bin_type=='quantile':
        if previous_binning is None:
            quantiles = np.empty(shape=(len(featuresHisto),), dtype=object)
            for i, feature in enumerate(featuresHisto):
                lCour=[]; print feature
                for ll in histogramme[feature]:
                    lCour.extend(np.array(ll)[np.where(np.isinf(np.array(ll))==False)])
                lCour=np.array(lCour)
                quantiles[i] = [scoreatpercentile(lCour, per) for per in [k*100/float(nb_bins_list[i]) for k in range(nb_bins_list[i]+1)]]
                
        else:
            quantiles = previous_binning
        #so we should have nb_bins_list[i]+1 values in quantiles[i] since we also give the rightmost edge
        
        nb_traj = len(histogramme[feature])
        matrice_result=np.empty(shape=(nb_traj, np.sum(nb_bins_list)), dtype=np.float64)
        matrice_result.fill(-1)
        #remembering the trajectories that does not have this feature
        toAverage={feature:[] for feature in featuresHisto}
        for traj_num in range(nb_traj):
            up_to=0
            for i,feature in enumerate(featuresHisto):
                ll=np.array(histogramme[feature][traj_num])[np.where(np.isinf(histogramme[feature][traj_num])==False)]
                if len(ll)==0:
                    toAverage[feature].append(traj_num)
                else:
                    zz=np.histogram(ll, bins=quantiles[i])[0]
                    assert np.sum(zz)==len(ll)
                    matrice_result[traj_num, up_to:up_to+nb_bins_list[i]]=zz/float(len(ll))
                up_to+=nb_bins_list[i]
    else:
        raise ValueError('Bin_type variable cannot take {} as argument'.format(bin_type))
    
    for i,feature in enumerate(featuresHisto):
        if toAverage[feature]!=[]:
            fr=np.sum(nb_bins_list[:i]); to=np.sum(nb_bins_list[:i+1])
            moy=np.mean(matrice_result[np.where(matrice_result[:,fr:to]>-0.5)[0],fr:to], 0)
            matrice_result[toAverage[feature], fr:to]=moy       
                                                                               
    if bin_type=='minmax':
        return matrice_result, minMax
    else:
        return matrice_result, quantiles  

if __name__ == '__main__':
    #MAINTENANT ON VOUDRAIT CALCULER LES DISTANCES ENTRE LES POINTS, ON VA REPARTIR CA SUR PLUSIEURS NOEUDS
    parser = OptionParser(usage="usage: %prog [options]")
    
    parser.add_option('--sim', type=int, dest='simulated')
    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option('-c', type=int, dest='choice', default=0)
    parser.add_option('-f', type=int, dest='intrand', default=0)
    
    parser.add_option("-d", type=int, dest="debut")
    parser.add_option("-l",type=int, dest="lamb")
    parser.add_option("-s", dest="size", type=int,default=1000)
    (options, args) = parser.parse_args()
#    if options.simulated:
#        f=open('../resultData/simulated_traj/histogramsNtotSim.pkl')
#    else:
#        print 'pas simulated'
#        f=open('..//histogramsNtotfirst.pkl')
#    data=pickle.load(f); f.close()
#    center=np.mean(data,0)
#    #
##Parallel(n_jobs=3, verbose=5)(delayed(foldingBig)(mesSolu, trainT, testT, n_big_fold, small=False) for n_big_fold in range(nb_big_folds))
#    for k in range(len(WEIGHTS)):
#        calculDistances(data, center, options.div_name, k, options.simulated) 
#    if options.intrand==0:
#        fichier= '../histogramsNtotfirst.pkl'
#    else:
#        fichier = '../rand{}_histogramsNtotfirst.pkl'.format(options.intrand)
#    
#    if options.choice==0:
#        print "We are going to compute squared distances of all points to points from {} to {} on file {}.".format(options.debut, options.debut+100, fichier)
#        f=open(fichier)
#        #r, histNtot, histApStot, histAvStot, histAvMtot, who,ctrlStatus, length, genes, sirna=pickle.load(f)
#        X=pickle.load(f)
#        f.close()
#        nb_bins_list=[50,50,50,50,50]; mat_hist_sizes=np.reshape(np.array(nb_bins_list), (1,5))
#    
#        distances = multSinkhornNumDistance(X[:int(options.size)], options.lamb, None, mat_hist_sizes, 12, X[options.debut:options.debut+100])
#        distances=np.sum(distances**2)
#    
#        f=open('TSS_debut{}_rand{}.pkl'.format(options.debut, options.intrand), 'w')
#        pickle.dump(distances, f); f.close()
#    else:
#        print "Computing correlations between different divergences on histograms on file {}".format(fichier)
#    #POUR CALCULER LES CORRELATIONS ENTRE LES DIFFERENTES DISTANCES SUR LES HISTOGRAMMES
#        f=open(fichier)
#        data=pickle.load(f); f.close()
#        
#        up_to=62
#        for name in featuresHisto[1:]:
#            print '--------{}'.format(name)
#            computingComparisons(data[:1000,up_to:up_to+50], name)
#            up_to+=50
