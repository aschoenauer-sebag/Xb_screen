import os, pdb, sys, time

import numpy as np
import cPickle as pickle

from matplotlib import pyplot as p
from collections import defaultdict
from itertools import product
from math import sqrt
from numpy.linalg import inv
from optparse import OptionParser
from warnings import warn

from peach import FuzzyCMeans
from sklearn.metrics import silhouette_score
from sklearn.cluster import MiniBatchKMeans, KMeans, SpectralClustering
from sklearn.mixture import GMM, DPGMM
from sklearn.manifold import Isomap
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from scipy.sparse.linalg import eigsh
from scipy.sparse import csc_matrix
from scipy.stats import pearsonr
from sklearn.utils.graph import graph_laplacian
from PyIBP import PyIBP as IBP
from util.kkmeans import KernelKMeans

from tracking.trajPack import featuresHisto, featuresNumeriques
from tracking.plots import plotClustInd, makeColorRamp, plotMovies, plotKMeansPerFilm, markers
from util.sandbox import cleaningLength, logTrsforming, subsampling, dist, histLogTrsforming, homeMadeGraphLaplacian
from util.listFileManagement import gettingSiRNA, expSi, siEntrez, typeD, typeD2, is_ctrl,\
    strToTuple
from util.plots import basic_colors, couleurs

from tracking.histograms import *
from util.kkmeans import KernelKMeans

#from joblib import Parallel, delayed, Memory

def correct_from_Nan(arr, perMovie):
#Si movement type a un coeff de correlation inf a 0.7 pour >=2 regressions, on supprime la trajectoire
    arr[np.where(arr[:,14]>1),11]=None
    toDel = []
    
#Si on regarde les trajectoires a l'echelle du movie on garde les trajectoires qui ont quelques NaN que l'on supprimera en regardant les distributions
    if perMovie:
        test=np.any(arr[:,-1]>=5)
        
#Sinon on supprime les trajectoires concernees pcq on ne pourra pas faire le clustering 
    else:
        test=(np.any(arr[:,-1]>=5) or np.any(np.isnan(arr)))
        
    if test:
        toDel = np.where(arr[:,-1]>=5)[0]
        if not perMovie:
            toDel=np.hstack((toDel, np.where(np.isnan(arr))[0]))
        arr=np.delete(arr, toDel, 0)
    arr=np.hstack((arr[:,:len(featuresNumeriques)+len(featuresHisto)], arr[:,-1, np.newaxis]))

    return arr, toDel

def histConcatenation(folder, exp_list, mitocheck, qc, filename = 'hist_tabFeatures_{}.pkl', verbose=0, hist=True, perMovie = False):
    who=[]; length=[]; r=[]; X=[]; Y=[]; ctrlStatus = []; sirna=[]; genes=[]
    time_length=[]; pbl_well=[]
    histNtot={nom:[] for nom in featuresHisto}

    yqualDict=expSi(qc)
    dictSiEntrez=siEntrez(mitocheck)

    for i, exp in enumerate(exp_list):
        print i,
        pl,w=exp
        if pl[:9]+'--'+w[2:5] not in yqualDict:
#i. checking if quality control passed
            sys.stderr.write("Quality control not passed {} {} \n".format(pl[:9], w[2:5]))
            continue   
        elif not is_ctrl((pl,w)) and yqualDict[pl[:9]+'--'+w[2:5]] not in dictSiEntrez:
#ii.checking if siRNA corresponds to a single target in the current state of knowledge
            sys.stderr.write( "SiRNA having no target or multiple target {} {}\n".format(pl[:9], w[2:5]))  
            continue
        try:
#iii.loading data
            f=open(os.path.join(folder, pl, filename.format(w)), 'r')
            arr, coord, histN= pickle.load(f)
            f.close()
        except IOError:
            sys.stderr.write("No file {}\n".format(os.path.join(pl, filename.format(w))))
        except EOFError:
            sys.stderr.write("EOFError with file {}\n".format(os.path.join(pl, filename.format(w))))
            pbl_well.append((pl, w))
        else:
            if arr==None:
                sys.stderr.write( "Array {} is None\n".format(os.path.join(pl, filename.format(w))))
                pbl_well.append((pl, w))
                continue
            elif len(arr.shape)==1 or arr.shape[0]<20:
                sys.stderr.write("Array {} has less than 20 trajectories. One needs to investigate why. \n".format(os.path.join(pl, filename.format(w))))
                pbl_well.append((pl, w))
                continue
            else:
                try:
                    arr, toDel = correct_from_Nan(arr, perMovie)

                    r= arr if r==[] else np.vstack((r, arr))
                    if hist:
                        for nom in featuresHisto:
                            histN[nom]=np.delete(np.array(histN[nom]), toDel)
                            histNtot[nom].extend(list(histN[nom]))
                    ll=arr.shape[0]
    
                except (TypeError, ValueError, AttributeError):
                    sys.stderr.write("Probleme avec le fichier {}".format(os.path.join(pl, filename.format(w))))
                else:   
                    time_length.extend([len(coord[k][0]) for k in filter(lambda x: x not in toDel, range(len(coord)))])
                    siCourant = yqualDict[pl[:9]+'--'+w[2:5]]
                    sirna.append(siCourant)
                    who.append((pl, w))
                    length.append(ll)
                    
                    if is_ctrl((pl,w)):
                        ctrlStatus.append(0)
                        genes.append('ctrl')
                    else:
                        ctrlStatus.append(1)
                        try:
                            genes.append(dictSiEntrez[siCourant])
                        except KeyError:
                            print "Key Error for experiment {} {}".format(pl, w)
                            continue
    if r ==[]:
        return None

#log trsforming data
    r2 = histLogTrsforming(r, verbose=verbose)        

    warn('The data was not normalized. Please check that it will be done before applying any algorithm.')
    
    return pbl_well, r2[:,:-1], histNtot,  who,ctrlStatus, length, genes, sirna, time_length

def concatenation(folder, exp_list, mitocheck, qc):
    who=[]; length=[]; r=[]; X=[]; Y=[]; ctrlStatus = []; genes=[]; sirna=[]

    yqualDict=expSi(qc)
    dictSiEntrez=siEntrez(mitocheck)
    
    whoCtrl=['00074_01', '00315_01']
    if exp_list ==None:
        exp_list=[]
        plates = os.listdir(folder)
        print len(plates)
        for pl in plates:
            listfiles = filter(lambda x: 'd_tabFeatures' in x, os.listdir(os.path.join(folder, pl)))
            
            for fichier in listfiles:
                exp_list.append((pl, fichier.split('_')[2]+'_01'))

        return concatenation(folder, exp_list)
    else:
        print "loading from experiments list"
        i=0
        for pl, w in exp_list:
            print i,
            i+=1
            try:
                f=open(os.path.join(folder, pl, 'd_tabFeatures_{}.pkl'.format(w)), 'r')
                arr= pickle.load(f)
                f.close()
            except IOError:
                print "Pas de fichier {}".format(os.path.join(pl, 'd_tabFeatures_{}.pkl'.format(w)))
            except EOFError:
                print "Probleme EOFError d'ouverture du fichier {}".format(os.path.join(pl, 'd_tabFeatures_{}.pkl'.format(w)))
                pdb.set_trace()
            else:
                if arr==None:
                    print "Probleme avec le fichier {}".format(os.path.join(pl, 'd_tabFeatures_{}.pkl'.format(w)))
                    pdb.set_trace()
                elif pl[:9]+'--'+w[2:5] not in yqualDict:
                    print "Quality control not passed", pl[:9], w[2:5]   
                elif w not in whoCtrl and yqualDict[pl[:9]+'--'+w[2:5]] not in dictSiEntrez:
                    print "SiRNA having no target or multiple target"           
                else:
                    try:
                        arr=arr[np.where(arr[:, -1]<5)]
                        arr=np.delete(arr, np.where(np.isnan(arr))[0], 0)
                    except (TypeError, EOFError, ValueError, AttributeError):
                        print "Probleme avec le fichier {}".format(os.path.join(pl, 'd_tabFeatures_{}.pkl'.format(w)))
                        pdb.set_trace()
                    else:   
                        r= arr if r==[] else np.vstack((r, arr))
      
                        who.append((pl, w))
                        length.append(arr.shape[0])
                        siCourant = yqualDict[pl[:9]+'--'+w[2:5]]
                        sirna.append(siCourant)

                        if w in whoCtrl:
                            ctrlStatus.append(0)
                            genes.append('ctrl')
                        else:
                            ctrlStatus.append(1)
                            try:
                                genes.append(dictSiEntrez[siCourant])
                            except KeyError:
                                pdb.set_trace()
        print r.shape
        
#log trsforming data
        r2 = logTrsforming(r)        
        
        print r2.shape, 'not normalized, cleaned from qc'
        return r2[:,:20], who,ctrlStatus, length, genes, sirna


#def cleaningCoord(traj_indices, arr, x,y):
#    allIndices = []
#    xc=np.array(x); yc=np.array(y)
#    for line in traj_indices[::-1]:
#        allIndices.extend(range(int(np.sum(arr[:line, 35])), int(np.sum(arr[:line+1, 35]))))
#    xc = np.delete(xc, allIndices, 0)
#    yc = np.delete(yc, allIndices, 0)
#    return xc, yc
#                arr=np.delete(arr, np.where(np.isnan(arr))[0], 0)


def silhouette(data, labels, N, metric='euclidean'):
    #N est le nb de fois que l'on veut tirer un echantillon des donnees pour calculer le score silhouette
    #Bien sur si la taille des donnes le permet on calcule le silhouette score sur tout le monde
    if data.shape[0]>12000:
        print 'We are going to compute the average of {} silhouette coefficient scores sampled in the data and labels produced'.format(N)
        result =[]
        for l in range(N):
            print "Time {}".format(l+1)    
            sData, sLabels = subsampling(data, labels) 
            result.append(silhouette_score(sData, sLabels, metric))
        return np.mean(result), np.std(result)
    else:
        print "We are going to compute the silhouette coefficient score."
        return silhouette_score(data, labels, metric)

def IBP_LinGaussianModel(data, num_samp):
    Klist=[]; likelist=[]; sigma_xlist=[]; sigma_alist=[]; alphalist=[]
    # IBP parameter (gamma hyperparameters)
    (alpha, alpha_a, alpha_b) = (1., 1., 1.)
# Observed data Gaussian noise (Gamma hyperparameters)
    (sigma_x, sx_a, sx_b) = (1., 1., 1.)
# Latent feature weight Gaussian noise (Gamma hyperparameters)
    (sigma_a, sa_a, sa_b) = (1., 1., 1.)
    print "data",data.shape
    buffet=IBP(data, (alpha, alpha_a, alpha_b), (sigma_x, sx_a, sx_b),(sigma_a, sa_a, sa_b))
    for s in range(num_samp):
        print "iteration ", s
        Klist.append(buffet.K)
        likelist.append(buffet.logLike())
        alphalist.append(buffet.alpha)
        sigma_xlist.append(buffet.sigma_x)
        sigma_alist.append(buffet.sigma_a)
        
        buffet.fullSample()
        
    return Klist, likelist, sigma_xlist, sigma_alist, alphalist
    

def K_Kmeans(data, mini, maxi, N=5, kernel='rbf', gamma=0.1,max_iter=1000, return_labels = False):
    #mini and maxi are respectively the min and max numbers of clusters that we want to try
    silhouette_r = []
    cohesion_r=[]
    for k in range(mini, maxi+1):
        print "--------------------------------------Cluster number : ", k
        model = KernelKMeans(n_clusters=k, max_iter=1000, tol=1e-3, random_state=None,
                 kernel=kernel, gamma=gamma, degree=3, coef0=1,
                 kernel_params=None, verbose=0)

        labels = model.fit_predict(data)
        final_dist = model.final_dist
#        print 'Computing silhouette score'
#        silhouette_r.append(silhouette(data, labels, N))
        print 'Computing cohesion'
#        centers = zou.cluster_centers_
        total_sum_squares = data.shape[0]*np.trace(model.gram_matrix)-np.sum(model.gram_matrix)
        cohesion_r.append(np.sum((final_dist)**2)/float(total_sum_squares))
    if return_labels:
        return silhouette_r, cohesion_r, labels
    else:
        return silhouette_r, cohesion_r           

def BatchKmeans(data, mini, maxi, N=5, metric='euclidean'):
    #mini and maxi are respectively the min and max numbers of clusters that we want to try
    silhouette_r = []
    cohesion_r=[]
    for k in range(mini, maxi+1):
        print "--------------------------------------Cluster number : ", k
        model = MiniBatchKMeans(n_clusters=k, batch_size = 1000, init='k-means++',n_init=1000,max_iter=1000, max_no_improvement=100, compute_labels = True)#100, 300 puis 500,500
        zou = model.fit(data)
        print 'Computing silhouette score'
        if k>1:
            silhouette_r.append(silhouette(data, zou.labels_, N))
        else: silhouette_r.append("undefined")
        
        print 'Computing cohesion'
        centers = zou.cluster_centers_
        total_sum_squares = np.sum((data - np.mean(data, 0))**2)
        cohesion=0
        for p in range(k):
            cohesion+=np.sum((data[np.where(zou.labels_ == p)]-centers[p])**2)
        cohesion_r.append(cohesion/float(total_sum_squares))
        
    return silhouette_r, cohesion_r    
    
#it is slower but can apparently be better
def Kmeans(data, mini, maxi, N=5, metric='euclidean', return_labels = False):
    #mini and maxi are respectively the min and max numbers of clusters that we want to try
    silhouette_r = []
    cohesion_r=[]
    for k in range(mini, maxi+1):
        print "--------------------------------------Cluster number : ", k
        model = KMeans(n_clusters=k, init='k-means++',n_init=100,max_iter=1000)#100, 300 puis 500,500
        #n_clusters=8, init='k-means++', n_init=10, max_iter=300, tol=0.0001, precompute_distances=True, verbose=0, random_state=None, copy_x=True, n_jobs=1, k=None)
        zou = model.fit(data)
        print 'Computing silhouette score'
        silhouette_r.append(silhouette(data, zou.labels_, N))
        
        print 'Computing cohesion'
        centers = zou.cluster_centers_
        total_sum_squares = np.sum((data - np.mean(data, 0))**2)

        cohesion=0
        for p in range(k):
            cohesion+=np.sum((data[np.where(zou.labels_ == p)]-centers[p])**2)
        cohesion_r.append(cohesion/float(total_sum_squares))
    if return_labels:
        return silhouette_r, cohesion_r, zou.labels_
    else:
        return silhouette_r, cohesion_r     
    
def CMeans(data, mini, maxi, m=2, eps=0.0001, n_iter=100, max_iter=1000, init='random'):
    partition_entropy_r = []#partition-entropy coefficient, see Validity_survey p. 31
    FS_r = []#Fukuyama-Sugeno index, see same article p.33    
    XB_r = []
    for k in range(mini, maxi+1):
        print 'on fit le modele sur les donnees, k= {}'.format(k)
        #we are going to run the algo n_iter times with different random init and then keep the better in terms of Xie-Beni index
        XB_courant= []; FS_courant = []; partition_entropy_courant=[]
        for iteration in range(n_iter):
            print iteration,
            mu=np.zeros(shape=(data.shape[0], k))
            #initialization: matrix Mu must be a stochastic matrix with no class with everyone or no one
            if init=='random':
                mu[:,0:int(k/2)]=0.5*np.random.rand(data.shape[0], int(k/2))
                if k%2==0: mu[:,k/2:]=0.5-mu[:,0:k/2]
                else:
                    mu[:,int(k/2):-1]=0.65*(0.5-mu[:,0:int(k/2)])
                    mu[:,-1]=1-np.sum(mu[:,:-1], 1)
                np.random.shuffle(mu); np.random.shuffle(mu.transpose())
            else:
                print "This initialization is not implemented yet", init
                raise AttributeError
#            if np.any(np.sum(mu, 1)!=np.ones(shape=(mu.shape[0],))):
#                raise
            model = FuzzyCMeans(data, mu, m)
            centers = model(emax=eps, imax=max_iter)
            tree=KDTree(centers)
            min_center_distance = min([max(tree.query(centers[l], 2)[0]) for l in range(k)])
            prob_labels = model.mu+10**(-15)
            FS=0; XB = 0
            for p in range(k):
                center = centers[p]
                FS+=np.sum((prob_labels[:, p]**m)*(np.sum((data-center)**2, 1)- np.linalg.norm(center)**2))
                XB+=np.sum((prob_labels[:, p]**2)*np.sum((data-center)**2, 1))
            XB/=(data.shape[0]*min_center_distance)
            XB_courant.append(XB); FS_courant.append(FS) 
            partition_entropy_courant.append(-1/float(data.shape[0]*np.log2(k))* np.sum(prob_labels * np.log2(prob_labels)))
        #sortie de la boucle d'iterations et choix des meilleurs indicateurs
        FS_r.append(min(FS_courant)); XB_r.append(min(XB_courant))
        #Partition entropy coefficient normalized by its max (not done for real and RR data)
        partition_entropy_r.append(min(partition_entropy_courant))
        
    return  FS_r, XB_r, partition_entropy_r 

def GaussianMixture(data, mini, maxi, covar='full'):
    BIC_r = []
    partition_entropy_r = []#partition-entropy coefficient, see Validity_survey p. 31
    score_r = []
    FS_r = []#Fukuyama-Sugeno index, see same article p.33
    print covar
    for k in range(mini, maxi+1):
        model = GMM(n_components=k, covariance_type=covar, random_state=None, thresh=0.01, min_covar=0.001, 
                        n_iter=10, n_init=10, params='wmc', init_params='wmc')
        print 'on fit le modele sur les donnees, k= {}'.format(k)
        model.fit(data)
        print 'on predit'
        prob_labels = model.predict_proba(data)+10**(-15)
        labels = model.predict(data)
        score = model.score(data); score_r.append((np.sum(score), np.mean(score)))
        BIC_r.append(model.bic(data)) 
        FS=0; sys.stderr.write("calculating FS index")
        for p in range(k):
            zz = data[np.where(labels==p)]
            if len(zz)==0:
                continue
            center = np.mean(zz, 0)
            FS+=np.sum(prob_labels[:, p]*(np.sum((data-center)**2, 1)- np.linalg.norm(center)**2))
        
        FS_r.append(FS);    
        #Partition entropy coefficient normalized by its max (not done for real and RR data)
        if k>1:
            partition_entropy_r.append(-1/float(data.shape[0]*np.log2(k))* np.sum(prob_labels * np.log2(prob_labels)))
        else: partition_entropy_r.append("undefined")
        
    return  BIC_r, FS_r, score_r, partition_entropy_r

def isomapping(data, mini, maxi, neighbours):
#    neighbours = [5, 10, 15, 20]
    error_r = dict(zip(neighbours, [[] for pp in range(4)]))
    corr_r = dict(zip(neighbours, [[] for pp in range(4)]))
    for n in neighbours:
        print 'neighbours :', n
        error_C = error_r[n]
        corr_C = corr_r[n]
        for k in range(mini, maxi+1):
            print 'components :', k
            m = Isomap(n_neighbors=n, n_components=k, eigen_solver='auto', tol=0, max_iter=None, path_method='auto', neighbors_algorithm='kd_tree')
            x = m.fit_transform(data)
            error_C.append(m.reconstruction_error()) 
            geod_d = m.dist_matrix_.flatten()
            new_euclid_d = cdist(x, x, metric='euclidean').flatten()
            corr_C.append(1- pearsonr(geod_d, new_euclid_d)[0]**2)
    #Parallel(n_jobs=3, verbose=5)(delayed(wholeProcess)(2.5, k, -1, args.output, loadingFolder) for k in range(nb_big_folds))
    return corr_r, error_r

def outputBin(data, ctrlSize,nbPheno, lPheno, binSize, sigma):#, nbDim=2, nbNeighbours=20):
    #m = Isomap(n_neighbors=nbNeighbours, n_components=nbDim, eigen_solver='auto', tol=0, max_iter=None, path_method='auto', neighbors_algorithm='kd_tree')
    #D = m.fit_transform(data)
    ctrl = data[:ctrlSize]
    ctrlTree = KDTree(ctrl, leafsize=10)
    length=ctrlSize
    
    mini = np.amin(data, 0); maxi=np.amax(data, 0); 
    nbPointsX = int((maxi[0]-mini[0])/float(binSize))+1
    nbPointsY = int((maxi[1]-mini[1])/float(binSize))+1
    
    result = np.zeros(shape=(nbPheno, nbPointsX, nbPointsY))
    denomCtrl = np.zeros(shape=(nbPointsX, nbPointsY)); denomCtrl.fill(0.001)
    print 'Parti'
    for pointX, pointY in product(range(nbPointsX), range(nbPointsY)):
        x=mini[0]+(pointX+0.5)*binSize; y=mini[1]+(pointY+0.5)*binSize
        ctrldou, ctrli = ctrlTree.query((x, y), ctrlSize, distance_upper_bound=binSize/sqrt(2))
        if min(ctrldou)<100:
            ctrlPoint = filter(lambda t: t[1]<ctrl.shape[0] and np.all(np.abs(ctrl[t[1]]-(x, y))<(binSize/2.0, binSize/2.0)), zip(ctrldou, ctrli))        
            for distance, cPoint in ctrlPoint:
                denomCtrl[pointX, pointY]+=dist((x,y), ctrl[cPoint], sigma)
                
    for ifilm in range(nbPheno):
        print 'film ', ifilm
        pheno = data[length:length+lPheno[ifilm]]
        phenoTree = KDTree(pheno, leafsize=10)
        
        for pointX, pointY in product(range(nbPointsX), range(nbPointsY)):
            x=mini[0]+(pointX+0.5)*binSize; y=mini[1]+(pointY+0.5)*binSize
            denom=denomCtrl[pointX, pointY]
            phenodou, phenoi=phenoTree.query((x, y), data.shape[0]-ctrlSize, distance_upper_bound=binSize/sqrt(2))
            if min(phenodou)<100:
                phenoPoint =filter(lambda t: t[1]<pheno.shape[0] and np.all(np.abs(pheno[t[1]]-(x, y))<(binSize/2.0, binSize/2.0)), zip(phenodou, phenoi))
                for distance, pPoint in phenoPoint:
                    local = dist((x,y), pheno[pPoint], sigma)
                    result[ifilm, pointX, pointY]+=local; denom+=local
        length+=lPheno[ifilm]        
        #if denom>0:result[ifilm, pointX, pointY]/=denom
    plotMovies('/media/lalil0u/New/workspace2/Xb_screen/resultData/features_on_films', result, 'pattern_b{}_s{}'.format(binSize, sigma))
    return result

def homeMadeSpectralClust(data, cluster_nb_min, cluster_nb_max, neighbours_nb, sig, show=False, affinity =None):
    '''
    Affinity if we provide an affinity matrix not using affinity(x,y)=exp(-sigma*||x-y||^2)
    '''
    #here we expect an n_samples x n_features array

    eigengap_r = []; sil_r=[]; coh_r=[]; variance_r = []
    BIC_r = []; FS_r = []; PEC_r=[];
    if affinity == None:
        print "Computing KDTree on data"
        tree = KDTree(data)
        W = np.zeros(shape=(data.shape[0], data.shape[0]), dtype=np.float)
        ## Put the values in the matrix according to k-NN approach + rbf similarity measure
        #ms si j'ai bien compris on peut aussi faire juste l'un ou juste l'autre, il faut regarder
        print 'Computing affinity matrix'
        for n in range(data.shape[0]):
            for k in range(n+1):
                #now using complete affinity matrix instead of kNN fucking bug i notice n months after writing this code. It was affinity between k and k...
                W[n, k]=dist(data[k], data[n], sig)
                W[k,n]=W[n,k]
#            dou, i=tree.query(data[n], k=neighbours_nb)
#            for couple in filter(lambda x: x[1]<np.inf, zip(i,dou)):
#                W[n,couple[0]] = dist(data[n], data[couple[0]], sig)
#                #we use kNN approach and not mutual kNN approach
#                W[couple[0],n] = W[n,couple[0]]
        D = np.diag(np.sum(W,0))
        Lrw = homeMadeGraphLaplacian(W, normed=True)
    else:
        Lrw = homeMadeGraphLaplacian(affinity, normed=True)
        
    if np.any(np.diagonal(Lrw)==0):
        outliersGraves = np.where(np.diagonal(Lrw)==0)
        dataBis = np.delete(data, outliersGraves, 0)
        Wbis = np.delete(np.delete(W, outliersGraves, 0), outliersGraves, 1)
        f=open("outliersGraves{}.pkl".format(time.clock()), 'w'); pickle.dump(data[outliersGraves],f); f.close()
        print "Important outliers that are only connected to themselves, I take them off"
        return homeMadeSpectralClust(dataBis, cluster_nb_min, cluster_nb_max, neighbours_nb, sig, show=show, affinity = Wbis)
    #nL, dd = graph_laplacian(W, normed=True, return_diag=True) #here we are computing Lsym using sklearn.utils. dd is the diagonal of the laplacian
    #we will need it to normalize the obtained matrix of eigenvectors (cf tutorial)
    #uL = D-W        #chaque ligne de u est donc l'embedding d'une ligne de data. Maintenant on peut faire KMeans dessus. On peut aussi plotter sur les deux premieres composantes
        #pour voir ce que ca donne

    print 'Diagonalizing normalized Laplacian'; Lrw=csc_matrix(Lrw)
    #la on est en shift invert mode j'ai copie sur sklearn mais je ne comprends pas
#    vals, vecs = eigsh(-L, k=cluster_nb+15, which = 'LM', tol=0.01, sigma=1.0, return_eigenvectors=True)
    #meme idee mais je comprends comment ca marche
    #nvals, nvecs = eigsh(nL, k=cluster_nb_max+10, which = 'LM', tol=0.01, sigma=0, return_eigenvectors=True)
    vals_rw, vecs_rw = eigsh(Lrw, k=cluster_nb_max+10, which = 'LM', tol=0.01, sigma=0, return_eigenvectors=True)
    
#par contre je ne comprends pas pourquoi la ca ne donne pas les memes valeurs propres
        #chaque ligne de u est donc l'embedding d'une ligne de data. Maintenant on peut faire KMeans dessus. On peut aussi plotter sur les deux premieres composantes
        #pour voir ce que ca donne

#    uvals, uvecs = eigsh(uL, k=cluster_nb+15, M=D, which = 'LM', tol=0.01, sigma=0, return_eigenvectors=True)
    #U = vecs.T[cluster_nb::-1] * dd#U is normalized here
    #chaque ligne de u est donc l'embedding d'une ligne de data. Maintenant on peut faire KMeans dessus. On peut aussi plotter sur les deux premieres composantes
    #pour voir ce que ca donne
    U = vecs_rw[:,:cluster_nb_max]#(nvecs /np.sqrt(np.sum(nvecs*nvecs, 0)))[:,:cluster_nb_max]#forme nb_points*cluster_nb
    for cluster_nb in range(cluster_nb_min, cluster_nb_max+1):
        #first indicator of wise cluster_nb choice : 
        eigengap_r.append(vals_rw[cluster_nb] - vals_rw[cluster_nb-1])
        
        #let's look at the Gaussian mixture result
        BIC, FS, score, partition_entropy = GaussianMixture(U, cluster_nb, cluster_nb)
        BIC_r.extend(BIC); FS_r.extend(FS); PEC_r.extend(partition_entropy)
        
        #let's look at the kmeans result
        sil, coh= Kmeans(U, cluster_nb, cluster_nb, return_labels=False); sil_r.extend(sil); coh_r.extend(coh); 
        
#        #computing correlation between affinity matrix and ideal affinity matrix containing clustering labels
#        ideal_W = np.zeros(shape=(data.shape[0], data.shape[0]))
#        ordre_cluster = []
#        for p in range(cluster_nb):
#            points = np.where(labels==p)[0]
#            ordre_cluster.extend(points)
#            for i, j in combinations(points, 2):  
#                ideal_W[i, i]=1
#                ideal_W[j,j]=1      
#                ideal_W[i,j]=1
#                ideal_W[j,i]=1
#        W = W[ordre_cluster,]; W=W[:,ordre_cluster]
#        #plotting that
#        method = None; metric = None; color = 'red_black_sky'
#        heatmap(W, range(0, W.shape[0]), range(0, W.shape[0]), method, method, metric, metric, color, "SpecClust_{}_{}_{}".format("Reality", neighbours_nb, cluster_nb), show)
#        ideal_W = ideal_W[ordre_cluster,]; ideal_W=ideal_W[:,ordre_cluster]
#        heatmap(ideal_W, range(0, W.shape[0]), range(0, W.shape[0]), method, method, metric, metric, color, "SpecClust_{}_{}_{}".format("Ideal", neighbours_nb, cluster_nb), show)
#        
#        variance_r.append(1-pearsonr(W.flatten(), ideal_W.flatten())[0]**2)
        
        
    return eigengap_r, BIC_r, FS_r, PEC_r, sil_r, coh_r, variance_r
        
        
def printFigures(folder, data, X, Y, length,indexWell):
    until = 0
    print indexWell
    if indexWell>0:
        until = np.sum(length[:indexWell])
    end = np.sum(length[:indexWell+1])
    num = int(end-until)
    f=p.figure(); ax=f.add_subplot(111);
#    mini = (min(X[until:end]), min(Y[until:end])), 
#    maxi=(max(X[until:end]), max(Y[until:end]))
    last=0
    for p in range(min(20, num)):
        indexTraj = p+until
        printTraj(ax, data, X, Y,length, indexTraj,inc=p,last=last)
        last += np.max(Y[np.sum(data[:indexTraj,35]):np.sum(data[:indexTraj+1, 35])])+max(0, -np.min(Y[np.sum(data[:indexTraj+1,35]):np.sum(data[:indexTraj+2, 35])]))+5
#printFigures(folder, data, X, Y,length, indexTraj = p+until, indexWell=None, inc=p, mini = (min(X[until:end]), min(Y[until:end])), 
#maxi=(max(X[until:end]), max(Y[until:end])))
#    f.savefig(os.path.join(folder, 'trajec_{}.png'.format(indexWell)))
    p.show()
    return 

def printTraj(ax, data, X, Y, length, indexTraj, inc=0, last=0):
    colors = makeColorRamp(20, basic_colors)
#    print len(colors)
#    print indexTraj
    until = np.sum(data[:indexTraj,35])
    end = np.sum(data[:indexTraj+1, 35])
    x = X[until:end]; y = Y[until:end]
    ax.plot(x, y+last, color=colors[inc])
    
    return
    
def importIBP(folder, baseName, iterations, sh=False):
    if type(baseName)==list:
        k, l, sx, sa, a=[],[],[],[],[]
        for name in baseName:
            print name
            kl, log, sigmax =importIBP(folder, name, iterations)
            k.append(kl); l.append(log); sx.append(sigmax); #sa.append(sigmaa); a.append(alphas)

        k.append("Feature number"); l.append("Log-likelihood"); sx.append('Variance of the design matrix')
        plotClustInd(folder, range(iterations),baseName, None,sh, k, l, sx)
    else:
        k, l, sx, sa, a=[],[],[],[],[]
        for p in range(50):
            try:
                f=open(os.path.join(folder, baseName+'{}_{}.pkl'.format(iterations, p)), 'r')
                kl, log, sigmax, sigmaa, alphas = pickle.load(f)
                f.close()
                print p,
            except:
                pass
            else:
                k.append(kl); l.append(log); sx.append(sigmax); sa.append(sigmaa); a.append(alphas)
    
        return k, l, sx

def importIso(folder, baseName, neighbours, sh=False):
    if type(baseName)==list:
        sil=[]; coh=[]
        for name in baseName:
            print name
            s,c=importIso(folder, name, neighbours)
            sil.append(s); coh.append(c)
        pdb.set_trace()
        sil.append("Residual variance"); coh.append("Reconstruction error")
        plotClustInd(folder, range(2,16),baseName, None,sh, sil, coh)
    else:
        sil=[]; coh=[]#pour corr_r, error_r
        for k in range(50):
            try:
                f=open(os.path.join(folder, baseName+'{}_{}.pkl'.format(neighbours,k)), 'r')
                s, c = pickle.load(f)
                f.close()
                print k,
            except:
                pass
            else:
                sil.append(s[neighbours]); coh.append(c[neighbours])
    
        return sil, coh
    
def compTimes(div_name,filenames, algo=1, iteration=0, simulated=False, sh=False,weight=0, folderList=['/media/lalil0u/New/workspace2/Tracking/resultData/sinkhornKMeans/results_PCA1',
            '/media/lalil0u/New/workspace2/Tracking/resultData/sinkhornKMeans/results_second']):
#    filename+='_a{}_k{}_d{}_w{}_bt{}_bs{}_{}.pkl'
    filename='{}_a{}_k{}_d{}_w{}_bt{}_bs{}_ct{}_{}.pkl'
#    elif file_type==-1:
#        #ici tout ce qui est en quantile a le cost_type value
#        filename+='_a{}_k{}_d{}_w{}_bt{}_bs{}_ik-means++_{}.pkl'
#    elif file_type==-2:
#        filename+='_a{}_k{}_d{}_w{}_bt{}_bs{}_{}.pkl'
    
    algo_names = ['K-means', "Mini batch K-means"]
    range_=range(2,20)

    if not simulated:
        data_name = 'real data'
        output_file = '{}_a{}_notSimPCA1.png'.format(div_name, algo)
    else:
        folder = '/media/lalil0u/New/workspace2/Tracking/resultData/simulated_traj/results_second'
        data_name = 'simulated data'
        output_file = '{}_a{}_Sim2.png'.format(div_name, algo)
    
    result={}
    min_ = 1
    for bt in ['quantile', 'minmax']:
        result[bt]={}
        if bt=='quantile': cost_choices=['num', 'val']
        else: cost_choices=["num"]
        for cost_type in cost_choices:
            result[bt][cost_type]={}
            for bs in [10,50]:
                result[bt][cost_type][bs]={k:np.zeros(shape=(len(range_),)) for k in range(len(folderList))}
                for k,folder_ in enumerate(folderList):
#                    print folder_
                    for l in range_:
                        try:
                            f=open(os.path.join(folder_, filename.format(filenames[k], algo, l,div_name[:5], weight, bt[:3], bs, cost_type, iteration)))
                            _,_,inertia=pickle.load(f); f.close()
                        except IOError:
                            print 'pas {},{}'.format(l,k)
                        else:
                            result[bt][cost_type][bs][k][l-2]=inertia
                            if inertia<min_: min_=inertia
                        
    fig=p.figure(figsize=(24,13)); i=0
    for bt in ['quantile', 'minmax']:
        if bt=='quantile': cost_choices=['num', 'val']
        else: cost_choices=["num"]
        for cost_type in cost_choices:
            for bs in [10,50]:
                ax=fig.add_subplot(2,3,i)
                for k,folder_ in enumerate(folderList):                                                                                                             
                    ax.scatter(range_, result[bt][cost_type][bs][k], color='#'+couleurs[k],marker=markers[k], label=str(folder_.split("/")[-1]))
                    ax.grid(linestyle='-')
                    ax.set_ylim(bottom=min_, top=1)
                    ax.set_title('Bin type {}, bin size {}, cost type {}'.format(bt, bs, cost_type))
                i+=1
    fig.suptitle("Divergence {} on {}, algo {}".format(div_name, data_name, algo_names[algo]))
    p.legend()
    if sh:
        p.show()
    else: 
        print 'image saved as {}'.format(os.path.join(folder, output_file))
        p.savefig(os.path.join(folder, output_file))
        
def calcStability(filename, div_name = 'transportation_exact', iteration = 5, sh=False, output_file=None,
                  ddimensional = False,
                 folder = "/media/lalil0u/New/workspace2/Xb_screen/resultData/sinkhornKMeans/results_PCA_MAR1"):
    
    datasize=29600
    weights = WEIGHTS if not ddimensional else ddWEIGHTS
    range_=range(2,15)
    filename=filename+'_a1_k{}_d{}_w{}_bt{}_bs{}_ct{}_{}.pkl' if not ddimensional else 'dd_'+filename+'_a1_k{}_d{}_w{}_bt{}_bs{}_ct{}_{}.pkl'
    pairs = [(0,1), (2,3), (4,0)]
    if output_file == None:
        result = {}
        
        for bt in ['quantile', 'minmax']:
            result[bt]={}
            if bt=='quantile' and "transportation" in div_name: cost_choices=['num', 'val']
            else: cost_choices=["num"]
            for cost_type in cost_choices:
                result[bt][cost_type]={}
                for bs in [10]:
                    result[bt][cost_type][bs]={k:np.zeros(shape=(len(range_)), dtype=list) for k in range(len(WEIGHTS))}
                    for k in range(len(weights)):
                        result[bt][cost_type][bs][k].fill([0])
                        for pair in pairs:
                            try:
                                seedname = 'seeds' if not ddimensional else 'dd_seeds'
                                seeds = '{}_{}_{}_s{}_{}_w{}_{}_{}.pkl'.format(seedname, cost_type[:3],bt, bs,div_name[:5],k,datasize, pair[0])
                                f=open(os.path.join(folder, seeds))
                                set1, _ = pickle.load(f)
                                f.close()
                                seeds = '{}_{}_{}_s{}_{}_w{}_{}_{}.pkl'.format(seedname, cost_type[:3],bt, bs,div_name[:5],k,datasize, pair[1])
                                f=open(os.path.join(folder, seeds))
                                set2, _ = pickle.load(f)
                                f.close()
                                
                            except IOError:
                                continue
                            
                            for l in range_:                            
                                try:
                                    f=open(os.path.join(folder, filename.format(l,div_name[:5], k, bt[:3], bs, cost_type, pair[0])))
                                    _,labels1,_=pickle.load(f); f.close()
                                    f=open(os.path.join(folder, filename.format(l,div_name[:5], k, bt[:3], bs, cost_type, pair[1])))
                                    _,labels2,_=pickle.load(f); f.close()
                                    
                                except IOError:
                                    pass
                                    #print 'pas {},{}, iter {}'.format(l,k, pair[0])
                                else:
                                    print 'working on {},{}, iter {}, {}'.format(l,k, pair[0], pair[1])
                                    intersection = filter(lambda x: x in set2, set1)
                                    indices1=[np.where(set1==i)[0][0] for i in intersection]
                                    indices2=[np.where(set2==i)[0][0] for i in intersection]
                                    labels1=labels1[indices1]; labels2=labels2[indices2]
                                    
                                    corr1=np.array([labels1==label for label in labels1])
                                    corr2=np.array([labels2==label for label in labels2])
                                    for i in range(corr1.shape[0]):
                                        corr1[i,i]=0; corr2[i,i]=0
                                    if result[bt][cost_type][bs][k][l-2] == [0]:
                                        result[bt][cost_type][bs][k][l-2]=[np.sum(corr1*corr2)/np.sqrt(np.sum(corr1*corr1)*np.sum(corr2*corr2))]
                                    else:
                                        result[bt][cost_type][bs][k][l-2].append(np.sum(corr1*corr2)/np.sqrt(np.sum(corr1*corr1)*np.sum(corr2*corr2)))
        f=open(os.path.join(folder, 'stability_calc.pkl'), 'w')
        pickle.dump(result, f); f.close()
    else:
        f=open(os.path.join(folder, output_file), 'r'); result = pickle.load(f); f.close()
    pdb.set_trace()   
    fig=p.figure(figsize=(24,13)); i=0
    for bt in ['quantile', 'minmax']:
        if bt=='quantile'  and "transportation" in div_name: cost_choices=['num', 'val']
        else: cost_choices=["num"]
        for cost_type in cost_choices:
            for bs in [10]:
                ax=fig.add_subplot(1,3,i)
                for k in range(len(weights)): 
                    if np.any(np.array([np.mean(result[bt][cost_type][bs][k][el-2]) for el in range_])==0):
                        ax.errorbar(range_, [np.mean(result[bt][cost_type][bs][k][el-2]) for el in range_], [np.std(result[bt][cost_type][bs][k][el-2]) for el in range_], 
                                color='#'+couleurs[k],fmt=markers[k], label=str(WEIGHTNAMES[k]))
                    else:
                        ax.errorbar(range_, [np.mean(result[bt][cost_type][bs][k][el-2]) for el in range_], [np.std(result[bt][cost_type][bs][k][el-2]) for el in range_], 
                               color='#'+couleurs[k],marker=markers[k], label=str(WEIGHTNAMES[k]))
                        
                                              
                    #ax.errorbar(range_, [np.mean(result[bt][cost_type][bs][k][el-2]) for el in range_], [np.std(result[bt][cost_type][bs][k][el-2]) for el in range_], 
                    #           color='#'+couleurs[k],marker=markers[k], label=str(WEIGHTNAMES[k]))
                    ax.grid(linestyle='-')
                    ax.set_ylim(bottom=0, top=1.1)
                    ax.set_title('Bin type {}, bin size {}, cost type {}'.format(bt, bs, cost_type))
                i+=1
    fig.suptitle("Divergence {} on {}, algo {}".format(div_name, 'real data', "mini-batch K-Means"))
    p.legend()
    if sh:
        p.show()
    else: 
        print 'image saved as {}'.format(os.path.join(folder, 'stability_calc.png'))
        p.savefig(os.path.join(folder, 'stability_calc.png'))


def importStability_monoParam(filenames,simulated=False, sh=False,
            folderList =['/media/lalil0u/New/workspace2/Tracking/resultData/sinkhornKMeans/results_PCA2']):
    
    algo_names = ['K-means', "Mini batch K-means"]
    range_=range(2,20)

    if not simulated:
        data_name = 'real data'
        output_file = '{}.png'.format(filenames[0][:-4])
    else:
        folder = '/media/lalil0u/New/workspace2/Tracking/resultData/simulated_traj/results_second'
        data_name = 'simulated data'
        output_file = '{}.png'.format(filenames[0][:-4])
    
    result={}
    min_ = 1
    #k=2
#    f=open(os.path.join(folder, 'seeds_num_quantile_s50_trans_w2_40000_0.pkl'))
#    tss=pickle.load(f)[-1]; f.close()
    for k,filename in enumerate(filenames):
        result[filename]=np.zeros(shape=(len(range_), 2))
        for l in range_:
            try:
                f=open(os.path.join(folderList[k], filename.format( l)))
                stability, inertia=pickle.load(f); f.close()
                
            except IOError:
                print 'pas {}'.format(l)
            else:
                result[filename][l-2]=[np.mean(stability), np.std(stability)]
                if np.min(stability)<min_: min_=np.min(stability)

    fig=p.figure(figsize=(24,13))
    ax=fig.add_subplot(1,1,1)
    for k,filename in enumerate(filenames):
        ax.errorbar(range_, result[filename][:,0], result[filename][:,1], color='#'+couleurs[k],marker=markers[k], label=filename)
        ax.grid(linestyle='-')
        ax.set_ylim(bottom=min_, top=1)
    #ax.set_title('Bin type {}, bin size {}, cost type {}'.format(bt, bs, cost_type))
    fig.suptitle('Comparaison on stability index')
    p.legend()
    if sh:
        p.show()
    else: 
        print 'image saved as {}'.format(os.path.join(folder, output_file))
        p.savefig(os.path.join(folder, output_file))

def importHistBatchKMeans_monoParam(filenames,simulated=False, sh=False,
            folderList =['/media/lalil0u/New/workspace2/Tracking/resultData/sinkhornKMeans/results_PCA1',
                         '/media/lalil0u/New/workspace2/Tracking/resultData/sinkhornKMeans/results_PCA2']):
    
    algo_names = ['K-means', "Mini batch K-means"]
    range_=range(2,20)

    if not simulated:
        data_name = 'real data'
        output_file = '{}.png'.format(filenames[0][:-4])
    else:
        folder = '/media/lalil0u/New/workspace2/Tracking/resultData/simulated_traj/results_second'
        data_name = 'simulated data'
        output_file = '{}.png'.format(filenames[0][:-4])
    
    result={}
    min_ = 1
    #k=2
#    f=open(os.path.join(folder, 'seeds_num_quantile_s50_trans_w2_40000_0.pkl'))
#    tss=pickle.load(f)[-1]; f.close()
    for k,filename in enumerate(filenames):
        result[filename]=np.zeros(shape=(len(range_), 2))
        for l in range_:
            try:
                f=open(os.path.join(folderList[k], filename.format( l)))
                load_=pickle.load(f); f.close()
    #            inertia=np.array(inertia)/float(tss)
    #            f=open(os.path.join(folder, filename.format( l)), 'w')
    #            pickle.dump([stability, inertia],f); f.close()
            except IOError:
                print 'pas {},{}'.format(l)
            else:
                if len(load_)==3:
                    inertia=load_[-1]
                    result[filename][l-2]=[inertia,0]
                    if inertia<min_: min_=inertia
                else:
                    inertia=load_[-1]
                    result[filename][l-2]=[np.mean(inertia), np.std(inertia)]
                    if np.min(inertia)<min_: min_=np.min(inertia)
    
    fig=p.figure(figsize=(24,13))
    ax=fig.add_subplot(1,1,1)
    for k,filename in enumerate(filenames):
        ax.errorbar(range_, result[filename][:,0], result[filename][:,1], color='#'+couleurs[k],marker=markers[k], label=filename)
        ax.grid(linestyle='-')
        ax.set_ylim(bottom=min_, top=1)
    #ax.set_title('Bin type {}, bin size {}, cost type {}'.format(bt, bs, cost_type))
    fig.suptitle('Comparaison on cohesion index')
    p.legend()
    if sh:
        p.show()
    else: 
        print 'image saved as {}'.format(os.path.join(folder, output_file))
        p.savefig(os.path.join(folder, output_file))

def importHistBatchKMeans(div_name,filename, algo=1, iteration=0, ddimensional = False, simulated=False, sh=False, kmax=15, max_division =False,
            folder = '/media/lalil0u/New/workspace2/Xb_screen/resultData/sinkhornKMeans/results_PCA_MAR1'):
#    filename+='_a{}_k{}_d{}_w{}_bt{}_bs{}_{}.pkl'
    filename=filename+'_a{}_k{}_d{}_w{}_bt{}_bs{}_ct{}_{}.pkl' if not ddimensional else 'dd_'+filename+'_a{}_k{}_d{}_w{}_bt{}_bs{}_ct{}_{}.pkl'
    weights = WEIGHTS if not ddimensional else ddWEIGHTS
#    elif file_type==-1:
#        #ici tout ce qui est en quantile a le cost_type value
#        filename+='_a{}_k{}_d{}_w{}_bt{}_bs{}_ik-means++_{}.pkl'
#    elif file_type==-2:
#        filename+='_a{}_k{}_d{}_w{}_bt{}_bs{}_{}.pkl'
    
    algo_names = ['K-means', "Mini batch K-means"]
    range_=range(2,kmax)

    if not simulated:
        data_name = 'real data'
        output_file = '{}_a{}_notSimPCA_MAR1.png'.format(div_name, algo)
    else:
        folder = '/media/lalil0u/New/workspace2/Tracking/resultData/simulated_traj/results_second'
        data_name = 'simulated data'
        output_file = '{}_a{}_Sim2.png'.format(div_name, algo)
    
    result={}
    min_ = 0;
    for bt in ['quantile', 'minmax']:
        result[bt]={}
        if bt=='quantile': cost_choices=['num', 'val']
        else: cost_choices=["num"]
        for cost_type in cost_choices:
            result[bt][cost_type]={}
            for bs in [10]:
                result[bt][cost_type][bs]={k:np.empty(shape=(len(range_),), dtype=list) for k in range(len(WEIGHTS))}
                for k in range(len(weights)):
                    print k
                    result[bt][cost_type][bs][k].fill([0]) 
                    for iter_ in range(iteration):
                        max_ = 1
                        for l in range_:
                            try:
                                f=open(os.path.join(folder, filename.format(algo, l,div_name[:5], k, bt[:3], bs, cost_type, iter_)))
                                _,_,inertia=pickle.load(f); f.close()
                            except IOError:
                                print 'pas {},{}, iter {}'.format(l,k, iter_)
                            else:
                                if max_division:    
                                    max_ = max(max_, inertia)
                                            
                        for l in range_:
                            try:
                                f=open(os.path.join(folder, filename.format(algo, l,div_name[:5], k, bt[:3], bs, cost_type, iter_)))
                                _,_,inertia=pickle.load(f); f.close()
                            except IOError:
                                print 'pas {},{}, iter {}'.format(l,k, iter_)
                            else:    
                                if result[bt][cost_type][bs][k][l-2] == [0]:
                                    result[bt][cost_type][bs][k][l-2]=[inertia/float(max_)]
                                else:
                                    result[bt][cost_type][bs][k][l-2].append(inertia/float(max_))

                            
    print result
    fig=p.figure(figsize=(24,13)); i=0
    for bt in ['quantile', 'minmax']:
        if bt=='quantile': cost_choices=['num', 'val']
        else: cost_choices=["num"]
        for cost_type in cost_choices:
            for bs in [10]:
                ax=fig.add_subplot(1,3,i)
                for k in range(len(weights)):#range(2,8):    
                    if result[bt][cost_type][bs][k][2]!=[0]:
#                        m= float(np.max( [np.max(result[bt][cost_type][bs][k][el-2]) for el in range_]))
                        m=1
                        if np.any(np.array([np.mean(result[bt][cost_type][bs][k][el-2]) for el in range_])==0):
                            ax.errorbar(range_, [np.mean(result[bt][cost_type][bs][k][el-2])/m for el in range_], [np.std(result[bt][cost_type][bs][k][el-2])/m**2 for el in range_], 
                                    color='#'+couleurs[k],fmt=markers[k], label=str(WEIGHTNAMES[k]))
                        else:
                            ax.errorbar(range_, [np.mean(result[bt][cost_type][bs][k][el-2])/m for el in range_], [np.std(result[bt][cost_type][bs][k][el-2])/m**2 for el in range_], 
                                   color='#'+couleurs[k],marker=markers[k], label=str(WEIGHTNAMES[k]))
                    else:
                        ax.errorbar(range_, [0 for zz in range_],[0 for zz in range_], color='#'+couleurs[k],fmt=markers[k],label=str(WEIGHTNAMES[k]))
    
                    ax.grid(linestyle='-')
                    ax.set_ylim(bottom=min_, top=1.1)
                    ax.set_title('Bin type {}, bin size {}, cost type {}'.format(bt, bs, cost_type))
                i+=1
    fig.suptitle("Divergence {} on {}, algo {}".format(div_name, data_name, algo_names[algo]))
    p.legend()
    if sh:
        p.show()
    else: 
        print 'image saved as {}'.format(os.path.join(folder, output_file))
        p.savefig(os.path.join(folder, output_file))
        
def compHistBatchKMeans(div_name,filename, algo=1, simulated=0, iteration=0, sh=False):
    filename1=filename+'_a{}_k{}_d{}_w{}_bt{}_bs{}_{}.pkl'
    filename2=filename+'_a{}_k{}_d{}_w{}_bt{}_bs{}_ik-means++_{}.pkl'
    algo_names = ['K-means', "Mini batch K-means"]
    #WEIGHTS=[0,1]
    range_=range(2,20)

    if not simulated:
        folder = '/media/lalil0u/New/workspace2/Tracking/resultData/sinkhornKMeans/results_second'
        data_name = 'real data'
        output_file = '{}_compGroundMetric_a{}_notSim.png'.format(div_name, algo)
    else:
        folder = '/media/lalil0u/New/workspace2/Tracking/resultData/simulated_traj/results_second'
        data_name = 'simulated data'
        output_file = '{}_compGroundMetric_a{}_Sim.png'.format(div_name, algo)
    
    dict_={i:{k:np.zeros(shape=(len(range_),)) for k in range(len(WEIGHTS))} for i in range(4)}
    i=0;min_ = 1
    bt='quantile'
    for filename in [filename1, filename2]:
        for bs in [10,50]:
            for k in range(len(WEIGHTS)):
                print k
                for l in range_:
                    try:
                        f=open(os.path.join(folder, filename.format(algo, l,div_name[:5], k, bt[:3], bs, iteration)))
                        _,_,inertia=pickle.load(f); f.close()
                    except IOError:
                        print 'pas {},{}'.format(l,k)
                    else:
                        dict_[i][k][l-2]=inertia
                        if inertia<min_: min_=inertia
            i+=1
                        
    fig=p.figure(figsize=(24,13)); i=0
    for filename in [filename1, filename2]:
        for bs in [10,50]:
            ax=fig.add_subplot(2,2,i)
            for k in range(len(WEIGHTS)):                                                                                                              
                ax.scatter(range_, dict_[i][k], color='#'+clustering.couleurs[k],marker=markers[k], label=str(WEIGHTNAMES[k]))
                ax.grid(linestyle='-')
                ax.set_ylim(bottom=min_, top=1)
                if filename==filename1:
                    ax.set_title('Bin size {}, ground metric = bin distance'.format(bs))
                else:
                    ax.set_title('Bin size {}, ground metric = bin value distance'.format(bs))
            i+=1
    fig.suptitle("Divergence {} on {}, algo {}, bin type quantile".format(div_name, data_name, algo_names[algo]))
    p.legend()
    if sh:
        p.show()
    else: 
        print 'image saved as {}'.format(os.path.join(folder, output_file))
        p.savefig(os.path.join(folder, output_file))

def importBatchKMeans(folder, baseName, sh=False):
    if type(baseName)==list:
        sil=[]; coh=[]
        for name in baseName:
            print name
            s, c = importBatchKMeans(folder, name)
            sil.append(s); coh.append(c)
        sil.append('Silhouette score'); coh.append("Intra-cluster cohesion")
        plotClustInd(folder, range(2,31),baseName, None,sh, sil, coh)
    else:
        sil=None; coh=None
        for k in range(50):
            #pdb.set_trace()
            try:
                f=open(os.path.join(folder, baseName+'{}.pkl'.format(k)), 'r')
                s, c = pickle.load(f)
                f.close()
                print k,
            except IOError:
                #print 'bou'
                pass
            else:
                sil=np.vstack((sil, np.array(s)[:,0])) if sil is not None else np.array(s)[:,0]
                coh=np.vstack((coh, np.array(c))) if coh is not None else np.array(c)
#                sil=np.vstack((sil, np.array(s)[:,0][:14])) if sil is not None else np.array(s)[:,0][:14]
 #               coh=np.vstack((coh, np.array(c)[:14])) if coh is not None else np.array(c)[:14]
                    
        return np.array(sil),np.array(coh)

def importCMeans(folder, baseName, sh=False):
    if type(baseName)==list:
        partition, FS, XB=[],[],[]
        for name in baseName:
            print name
            p, fs, xbi= importCMeans(folder, name)
            partition.append(p); FS.append(fs); XB.append(xbi)
        partition.append("Partition entropy coefficient"); FS.append('Fukuyama-Sugeno index')
        XB.append('Xie-Beni index');
        plotClustInd(folder, range(2,31),baseName, None,sh, XB, FS, partition)
    else:
        partition, FS, XB = None, None, None#partition-entropy coefficient, see Validity_survey p. 31
        for k in range(50):
            try:
                nom = baseName+'{}.pkl'.format(k)
                f=open(os.path.join(folder, nom), 'r')
                fs, xbi, p = pickle.load(f)
                f.close()
                print k,
            except:
                pass
            else:
                XB = np.vstack((XB, np.array(xbi))) if XB is not None else np.array(xbi)
                partition =np.vstack((partition, np.array(p))) if partition is not None else np.array(p)
                FS = np.vstack((FS, np.array(fs))) if FS is not None else np.array(fs)
    
        return partition, FS, XB

def importGaussianM(folder, baseName, sh=False):
    if type(baseName)==list:
        BIC, partition, FS, sc=[],[],[],[]
        for name in baseName:
            print name
            b, p, fs, score = importGaussianM(folder, name)
            BIC.append(b); partition.append(p); FS.append(fs); sc.append(score)
        BIC.append('Bayesian information criterion'); partition.append("Partition entropy coefficient"); FS.append('Fukuyama-Sugeno index')
        sc.append('Log-probability of data under model');
        plotClustInd(folder, range(2,31),baseName, None,sh, BIC, FS, sc)
    else:
        BIC, partition, FS, sc =None, None, None, None#partition-entropy coefficient, see Validity_survey p. 31
        for k in range(50):
            try:
                nom = baseName+'{}.pkl'.format(k)
                f=open(os.path.join(folder, nom), 'r')
                b, fs, score, p = pickle.load(f)
                f.close()
                print k,
            except:
                pass
            else:
                sc = np.vstack((sc, np.array(score)[:,0])) if sc is not None else np.array(score)[:,0]
                BIC = np.vstack((BIC, np.array(b))) if BIC is not None else np.array(b)
                partition =np.vstack((partition, np.array(p))) if partition is not None else np.array(p)
                FS = np.vstack((FS, np.array(fs))) if FS is not None else np.array(fs)
    
        return BIC, partition, FS, sc


def importSpecClust(folder, baseName, neighbours, sigma, show=False):
    if type(baseName)==list:
        eigen_gap = []; sil=[]; coh=[];BIC = []; FS = []; partition=[];
        for name in baseName:
            for si in sigma:
                #for n in neighbours:
                n=neighbours[0]
                e, b, fs, p, s, c = importSpecClust(folder, name, n, si)
                eigen_gap.append(e)
                BIC.append(b); partition.append(p); FS.append(fs)
                sil.append(s); coh.append(c)
        eigen_gap.append('Eigengap, sigma ={}'.format(si))
        sil.append('Silhouette score, sigma ={}'.format(si)); coh.append("Intra-cluster cohesion")
        BIC.append('Bayesian information criterion, sigma ={}'.format(si)); partition.append("Partition entropy coefficient"); FS.append('Fukuyama-Sugeno index')
        plotClustInd(folder, range(2,31),baseName, sigma,show, eigen_gap ,sil, coh)
        plotClustInd(folder, range(2,31),baseName, sigma,show, FS, partition, BIC)
    else:  
        eigen_gap=None
        BIC, partition, FS =None, None, None
        sil=None; coh=None
        for k in range(50):
            try:
                #print os.path.join(folder, baseName+'{}_{}_{}.pkl'.format(neighbours,sigma, k))
                f=open(os.path.join(folder, baseName+'{}_{}_{}.pkl'.format(neighbours,sigma, k)), 'r')
                e, b, fs, p, s, c, v = pickle.load(f)
                f.close()
            except IOError:
                pass
            else:
                #print np.array(e).shape
                eigen_gap=np.vstack((eigen_gap, np.array(e))) if eigen_gap is not None else np.array(e)
                BIC = np.vstack((BIC, np.array(b))) if BIC is not None else np.array(b)
                partition =np.vstack((partition, np.array(p))) if partition is not None else np.array(p)
                FS = np.vstack((FS, np.array(fs))) if FS is not None else np.array(fs)
                sil=np.vstack((sil, np.array(s))) if sil is not None else np.array(s)
                coh=np.vstack((coh, np.array(c))) if coh is not None else np.array(c)
    
        return eigen_gap, BIC, FS, partition, sil, coh

if __name__ == '__main__':
    description =\
'''
%prog - Parallel clustering
You can in particular set up the noise level
'''
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    
    parser.add_option('--action', type=str, default=None)
    parser.add_option('--outputname', type=str, default='results_whole_iter5_median_015')
    (options, args) = parser.parse_args()

    
    if options.action=='collectingTrajectories':
        f=open(os.path.join('../resultData/features_on_films/', options.outputname+'_hit_exp.pkl'))
        l=pickle.load(f)
        f.close()
        
        ll=strToTuple(l, os.listdir('/share/data20T/mitocheck/tracking_results'))
        
        print 'this is launched'
        _, r, _,  who,ctrlStatus, length, genes, sirna, time_length=histConcatenation('/share/data20T/mitocheck/tracking_results/', 
                                ll, '../data/mitocheck_siRNAs_target_genes_Ens75.txt', '../data/qc_export.txt', hist=False, perMovie=False)
        f=open(os.path.join('../resultData/features_on_films/', options.outputname+'hit_exp_data.pkl'), 'w')
        pickle.dump([r,  who,ctrlStatus, length, genes, sirna, time_length], f)
        f.close()
        
        from tracking.trajPack import featuresSaved
        from sklearn.decomposition import PCA
        
        r=np.hstack((r[:,:len(featuresNumeriques)], r[:,featuresSaved.index('mean persistence'), np.newaxis], r[:, featuresSaved.index('mean straight'), np.newaxis]))
        print r.shape

        pca=PCA(n_components=r.shape[1])
        nr=(r-np.mean(r,0))/np.std(r,0)
        pcaed=pca.fit_transform(nr)
        pcaed=pcaed/np.std(pcaed,0)
        f=open(os.path.join('../resultData/features_on_films/', options.outputname+'_pcaed.pkl'), 'w')
        pickle.dump((pca, pcaed),f); f.close()
    
#    f=open('../resultData/features_on_films/hit_experiments_siRNAhighconf_PCAed_data.pkl')
#    narr=pickle.load(f)
#    f.close()
    
#    print 'BatchKMeans'
#    silhouette_r, cohesion_r = BatchKmeans(narr, 2, 30, N=10)
#    f=open('../resultData/features_on_films/batchKMeans_hitexp_highconfsiRNAs_1.pkl', 'w')
#    pickle.dump([silhouette_r, cohesion_r], f); f.close()
    
#    parser.add_option("-f", "--folder", dest="folder", default='/cbio/donnees/aschoenauer/workspace2/Tracking/resultData/',
#                      help="The folder where the data is and everything saved")
#    
#    parser.add_option("-a", "--algo", dest="algo", default = 0, 
#                      help="0 for Isomap\
#                       1 for Gaussian mixtures\
#                       2 for K-means\
#                       3 for mini-batch K-means\
#                       4 for fuzzy C-means\
#                       5 for IBP")
#
#    parser.add_option("-n", "--n_iter", dest="n_iter", default = 0, 
#                      help="Nb de fois")
#    
#    parser.add_option("-g", "--neighbours", dest="neighbours", default = 5, type=int,
#                      help="Nb de voisins")
#    
#    parser.add_option("-d", "--data", dest="data", default = 0, 
#                      help="Quelles donnees")
#    
#    parser.add_option("-b", "--noise", dest="noise",type=float, default = 0.1, 
#                      help="Noise level")
#    
#    parser.add_option("-s", "--sigma", dest="sigma",type=float, default = 0.01, 
#                      help="Noise level")
#    
#    parser.add_option("--density", dest="dens",type=int, default = 5, 
#                      help="Max nb of neighbours during track")
#    parser.add_option("--covariance", dest="covar",type=str, default = 'full', 
#                      help="Covariance type for GMM algo")
#    parser.add_option("--fuzzy", dest="fuzzifier",type=int, default = 2, 
#                      help="Fuzzifier parameter for Fuzzy CMeans algo")
#    
#    parser.add_option("--numsampling", dest="num_samp",type=int, default = 10, 
#                      help="Number of iterations for Gibbs sampling to estimate parameters for IBP")
#
#    (options, args) = parser.parse_args()
#
##    folder = '/cbio/donnees/aschoenauer/data/tracking/migration'
## PUIS features with densities: folder = '/share/data/mitocheck/tracking_results'
##    allData, whoIsInData, length = concatenation(folder)
#    names = dict(zip(range(7), ['isomap', 'gaussianM', 'kMeans', 'batchKMeans', 'SpecClust', 'cMeans', 'IBP']))
#    baseNames = dict(zip(range(4), ['', 'd_p_trsfRR_', 'd_p_trsfRG_0.5out_', 'real_'])); baseName = baseNames[int(options.data)]
##FOR DMAX50<5
#    f_size = len(features1)
#    s = 500000    
#    if int(options.algo)==0 and baseName+names[int(options.algo)]+'{}_{}.pkl'.format(options.neighbours, options.n_iter) in os.listdir(options.folder):
#            print 'already treated {}'.format(baseName+names[int(options.algo)])
#            sys.exit()
#    elif int(options.algo)==4 and baseName+names[int(options.algo)]+'{}_{}_{}.pkl'.format(options.neighbours,options.sigma, options.n_iter) in os.listdir(options.folder):
#            print 'already treated {}'.format(baseName+names[int(options.algo)])
#            sys.exit()
#    elif baseName+names[int(options.algo)]+'{}.pkl'.format(options.n_iter) in os.listdir(options.folder):
#            print 'already treated {}'.format(baseName+names[int(options.algo)])
#            sys.exit()
#    elif int(options.data)==0:
#        print "Working with whole real data unplanned. Pgm exiting"
#        sys.exit()
##        f=open(os.path.join(options.folder, 'data_norm.pkl'), 'r')
##        r, X, Y, who,ctrlStatus, length = pickle.load(f); f.close()
##        dataCtrl = concatCtrl(r, ctrlStatus, length, len(ctrlStatus)-np.sum(ctrlStatus))
##        np.random.shuffle(dataCtrl)
##        ctrl = dataCtrl[:1000]################
##        dataPheno = returnPheno(r, ctrlStatus, length, np.sum(ctrlStatus))
##        np.random.shuffle(dataPheno)
##        dPheno = dataPheno[0]
##        for d in dataPheno: 
##            dPheno= np.vstack((dPheno, d))
##            if dPheno.shape[0]>9000:
##                break
##        currData = np.vstack((ctrl, dPheno[:9000]))###############
#
#    elif int(options.data)==1:
#        print "Working with random data (uniform), shape {} x {}".format(s,f_size)
#        f=open(os.path.join(options.folder, 'l123_minMax.pkl'), 'r')
#        mini, maxi = pickle.load(f); f.close()
#
#        r = np.random.random_sample(size=(s, f_size))*(maxi-mini)+mini
#        r = (r-np.mean(r, 0))/np.std(r, 0)
#        currData=r[:15000]
#        
#    elif int(options.data)==2:
#        print "Working with random data (Gaussian), noise level {}, shape {} x {}".format(options.noise,s,f_size)
#        f=open(os.path.join(options.folder, 'l123_mean.pkl'), 'r')
#        mean = pickle.load(f); f.close()
#        f=open(os.path.join(options.folder, 'l123_std.pkl'), 'r')
#        stdev = pickle.load(f); f.close()
#        f=open(os.path.join(options.folder, 'l123_minMax.pkl'), 'r')
#        mini, maxi = pickle.load(f); f.close()
#        
#        outliers_percentage = 0.2
#        r = np.random.normal(mean, stdev, size=(s, f_size))
#        
#        data_size = s if int(options.algo) in [2,3] else 15000
#            
#        outliers = int(np.random.rand(1)*outliers_percentage*data_size); left = data_size-outliers
#    
#        currData = np.vstack((r[:left], np.random.random_sample(size=(outliers, f_size))*(maxi-mini)+mini))
#        currData=(currData-np.mean(currData, 0))/np.std(currData, 0)
#        currData += np.random.normal(0, float(options.noise), size=(currData.shape[0], f_size))
#        
#    elif int(options.data)==3:
#        print "Working with all real data, only {} features, max nb of neighbours {}".format(f_size, options.dens)
#        f=open(os.path.join(options.folder, 'small_l123.pkl'), 'r')#'trsf_maxd50inf{}_data.pkl'.format(options.dens)), 'r')
#        zou= pickle.load(f); f.close()
#        r=zou-np.mean(zou, 0)
#        #print mot
#        #r=r[:,:-4]
#        if int(options.algo) not in [2,3]:
#            np.random.shuffle(r); currData = r[:15000]
##            dataCtrl = concatCtrl(r, ctrlStatus, length, len(ctrlStatus)-np.sum(ctrlStatus))
##            np.random.shuffle(dataCtrl)
##            ctrl = dataCtrl[:5000]################
##            dPheno = concatCtrl(r, 1-np.array(ctrlStatus), length, np.sum(ctrlStatus))
##    #pas fait pour les premieres experiences sur trsformed data
##            np.random.shuffle(dPheno)
##            
##            currData = np.vstack((ctrl, dPheno[:10000]))###############
#        
#    else:
#        print "Settings not understood"
#        sys.exit()
#        
#    neighbours = [int(options.neighbours)]
#    if int(options.algo) ==0:
#        print 'isomapping, neighbours = ', neighbours
#        corr_r, error_r = isomapping(currData, 2, 15, neighbours)
#        f=open(os.path.join(options.folder, baseName+'isomap{}_{}.pkl'.format(options.neighbours, options.n_iter)), 'w')
#        pickle.dump([corr_r, error_r], f); f.close()
#        
#    elif int(options.algo) ==1:
#        print 'gaussian mixturing, covariance type {}'.format(options.covar)
#        BIC_r, FS_r, score_r, partition_entropy_r = GaussianMixture(currData, 2, 30, covar = options.covar)
#        f=open(os.path.join(options.folder, baseName+'gaussianM{}_{}.pkl'.format(options.covar, options.n_iter)), 'w')
#        pickle.dump([BIC_r, FS_r, score_r, partition_entropy_r], f); f.close()
#    elif int(options.algo) ==2:
#        print "KMeans"
#        silhouette_r, cohesion_r = Kmeans(r, 2, 15, N=10)
#        f=open(os.path.join(options.folder, baseName+'kMeans{}.pkl'.format(options.n_iter)), 'w')
#        pickle.dump([silhouette_r, cohesion_r], f); f.close()
#    elif int(options.algo) ==3:
#        print 'BatchKMeans'
#        silhouette_r, cohesion_r = BatchKmeans(r, 2, 30, N=10)
#        f=open(os.path.join(options.folder, baseName+'batchKMeans{}.pkl'.format(options.n_iter)), 'w')
#        pickle.dump([silhouette_r, cohesion_r], f); f.close()
#        
#    elif int(options.algo) ==4:
#        print 'Spectral Clustering, neighbours = {}, sigma = {}'.format(options.neighbours, options.sigma)
#        r = homeMadeSpectralClust(currData, 2, 30, int(options.neighbours), options.sigma, show=False)
#        f=open(os.path.join(options.folder, baseName+'SpecClust{}_{}_{}.pkl'.format(options.neighbours,options.sigma, options.n_iter)), 'w')
#        pickle.dump(r, f); f.close()
#        
#    elif int(options.algo) ==5:
#        print 'Fuzzy C-Means, fuzzifier {}'.format(options.fuzzifier)
#        r = CMeans(currData, 2, 30, options.fuzzifier)
#        f=open(os.path.join(options.folder, baseName+'CMeans{}_{}.pkl'.format(options.fuzzifier, options.n_iter)), 'w')
#        pickle.dump(r, f); f.close()
#        
#    elif int(options.algo) ==6:
#        print 'Indian Buffet Process'
#        r = IBP_LinGaussianModel(currData, options.num_samp)
#        f=open(os.path.join(options.folder, baseName+'IBP{}_{}.pkl'.format(options.num_samp,options.n_iter)), 'w')
#        pickle.dump(r, f); f.close()
#        
#        
#        
        
        
        
        
        
        
        
        
        
        
        
        #CODE POUR LANCER LES ALG DE CLUSTERING SANS LES INFOS DE DENSITE
        
        
#        
#    names = dict(zip(range(5), ['isomap', 'gaussianM', 'kMeans', 'batchKMeans', 'SpecClust' ]))
#    baseNames = dict(zip(range(4), ['', 'pRR_', 'pRG_', 'pure_'])); baseName = baseNames[int(options.data)]
#    f_size = len(sureFeatures)
#    
#    if int(options.algo)==0 and baseName+names[int(options.algo)]+'{}_{}.pkl'.format(options.neighbours, options.n_iter) in os.listdir(options.folder):
#            print 'already treated {}'.format(baseName+names[int(options.algo)])
#            sys.exit()
#    elif int(options.algo)==4 and baseName+names[int(options.algo)]+'{}_{}_{}.pkl'.format(options.neighbours,options.sigma, options.n_iter) in os.listdir(options.folder):
#            print 'already treated {}'.format(baseName+names[int(options.algo)])
#            sys.exit()
#    elif baseName+names[int(options.algo)]+'{}.pkl'.format(options.n_iter) in os.listdir(options.folder):
#            print 'already treated {}'.format(baseName+names[int(options.algo)])
#            sys.exit()
#    elif int(options.data)==0:
#        print "Working with real data"
#        f=open(os.path.join(options.folder, 'data_norm.pkl'), 'r')
#        r, X, Y, who,ctrlStatus, length = pickle.load(f); f.close()
#        dataCtrl = concatCtrl(r, ctrlStatus, length, len(ctrlStatus)-np.sum(ctrlStatus))
#        np.random.shuffle(dataCtrl)
#        ctrl = dataCtrl[:1000]################
#        dataPheno = returnPheno(r, ctrlStatus, length, np.sum(ctrlStatus))
#        np.random.shuffle(dataPheno)
#        dPheno = dataPheno[0]
#        for d in dataPheno: 
#            dPheno= np.vstack((dPheno, d))
#            if dPheno.shape[0]>9000:
#                break
#        currData = np.vstack((ctrl, dPheno[:9000]))###############
#
#    elif int(options.data)==1:
#        print "Working with random data (uniform), shape 300 000 x {}".format(f_size)
#        f=open(os.path.join(options.folder, 'minMax.pkl'), 'r')
#        mini, maxi = pickle.load(f); f.close()
#        s = 300000
#        r = np.random.random_sample(size=(s, f_size))*(maxi[sureFeatures]-mini[sureFeatures])+mini[sureFeatures]
#        r = (r-np.mean(r, 0))/np.std(r, 0)
#        currData=r[:10000]
#        
#    elif int(options.data)==2:
#        print "Working with random data (Gaussian), noise level {}, shape 300 000 x 13".format(options.noise)
#        f=open(os.path.join(options.folder, 'mean.pkl'), 'r')
#        mean = pickle.load(f)[sureFeatures]; f.close()
#        f=open(os.path.join(options.folder, 'stdev.pkl'), 'r')
#        stdev = pickle.load(f)[sureFeatures]; f.close()
#        f=open(os.path.join(options.folder, 'minMax.pkl'), 'r')
#        mini, maxi = pickle.load(f); f.close()
#        s = 300000; outliers_percentage = 0.1
#        r = np.random.normal(mean, stdev, size=(s, f_size))
#        
#        data_size = s if int(options.algo) in [2,3] else 10000
#            
#        outliers = int(np.random.rand(1)*outliers_percentage*data_size); left = data_size-outliers
#    
#        currData = np.vstack((r[:left], np.random.random_sample(size=(outliers, f_size))*(maxi[sureFeatures]-mini[sureFeatures])+mini[sureFeatures]))
#        currData=(currData-np.mean(currData, 0))/np.std(currData, 0)
#        currData += np.random.normal(0, float(options.noise), size=(currData.shape[0], f_size))
#        
#    elif int(options.data)==3:
#        print "Working with all real data, only 13 features"
#        f=open(os.path.join(options.folder, 'pure_data_norm.pkl'), 'r')
#        r, X, Y, who,ctrlStatus, length = pickle.load(f); f.close()
#
#        dataCtrl = concatCtrl(r, ctrlStatus, length, len(ctrlStatus)-np.sum(ctrlStatus))
#        np.random.shuffle(dataCtrl)
#        ctrl = dataCtrl[:1000]################
#        dataPheno = returnPheno(r, ctrlStatus, length, np.sum(ctrlStatus))
#        np.random.shuffle(dataPheno)
#        dPheno = dataPheno[0]
#        for d in dataPheno: 
#            dPheno= np.vstack((dPheno, d))
#            if dPheno.shape[0]>9000:
#                break
#        currData = np.vstack((ctrl, dPheno[:9000]))###############
#        
#    else:
#        print "Settings not understood"
#        sys.exit()
    
