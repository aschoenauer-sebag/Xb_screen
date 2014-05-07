import numpy as np
import cPickle as pickle
import os, pdb, getpass, sys
sys.path.append("/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc")
if getpass.getuser()=='aschoenauer':
    import matplotlib as mpl
    mpl.use('Agg')
    show = False
elif getpass.getuser()=='lalil0u':
    import matplotlib as mpl
    show=True
import matplotlib.pyplot as p
from collections import Counter
from itertools import product
from scipy.stats import pearsonr
from optparse import OptionParser

from sklearn.cluster import MiniBatchKMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap
from sklearn.metrics import silhouette_score
from sklearn.mixture import GMM
from scipy.spatial.distance import cdist
from scipy.sparse.linalg import eigsh
from scipy.sparse import csc_matrix
from scipy.spatial import KDTree

from peach import FuzzyCMeans


from util.fileManagement import ArffReader
from util.kkmeans import KernelKMeans
from util.plots import couleurs, markers, makeColorRamp
from util.sandbox import dist, homeMadeGraphLaplacian
from nucleus import methods, params, extended_params, param_sets, pca_param, whiten_param, outFile

def getData(filename, pca=26, whiten=False):
    '''
    From the name of an Arff file, reads the data into it and give back the data matrix as well as the labels list
    Actions on the data :
    - shuffle
    - normalization
    - standardization
    - if asked, PCAed and possibly whitened
    '''
    reader = ArffReader(filename)
    data = None
    labels = []
    for el in reader.dctFeatureData:
        data=reader.dctFeatureData[el] if data is None else np.vstack((data, reader.dctFeatureData[el]))
        labels.extend([el for k in range(len(reader.dctFeatureData[el]))])
    
    data_matrix = np.hstack((data, np.array(labels)[:,np.newaxis]))
    np.random.shuffle(data_matrix)    
    data = np.asarray(data_matrix[:,:-1], dtype=np.float64)
    labels = data_matrix[:,-1]
    #Deleting 2 elements for which h8 features are NaN
    labels = np.delete(labels, np.where(np.isnan(data))[0],0)
    data = np.delete(data, np.where(np.isnan(data))[0],0)
    
    #We standardize the matrix prior to PCA to rescale the different features
    data=(data-np.mean(data, 0))/np.std(data,0)
    if type(pca)==int:
        #perform PCA
        print "Performing PCA, retaining {} components".format(pca)
        pca = PCA(n_components = pca, whiten =whiten)
        data = pca.fit_transform(data)
        
    print 'Number of elements ', data.shape[0]
    print 'Number of features ', data.shape[1]
    print 'Labels ', reader.dctClassNames.values()

    return reader.dctClassLabels.keys(), data, labels

def givingAllParameterSets(algo=None, l=None):
    if algo is not None:
        paramL=param_sets[algo]
    else:
        paramL = l
        
    print paramL
    param_set=[]; param_dict = []
    for x in paramL:
        if type(x)==int or type(x)==str:
            param_set.append([x])
        elif type(x)==list:
            param_set.append(x)
        elif type(x)==dict:
            subset=[]
            for k in x:
                subParamL = givingAllParameterSets(l=x[k])
                for ll in subParamL:
                    ll.insert(0, k)
                subset.extend(subParamL)
            param_dict.extend(subset)
    if param_dict == []:
        return [list(el) for el in  product(*param_set)]
    else:
        result0 = [list(el) for el in  product(*param_set)]
        result = []
        for el in result0:
            for le in param_dict:
                new = list(el)
                new.extend(le)
                result.append(new)
        return result
    
def savingParameterSets(algo, outFolder = '../resultData/nucleus/parameters'):
    if not os.path.isdir(outFolder):
        os.makedirs(outFolder)
    
    parameter_set = givingAllParameterSets(algo)
    for i, paramL in enumerate(parameter_set):
        d = dict(zip(extended_params[algo], paramL))
        f=open(os.path.join(outFolder, 'parameters_{}_{}.pkl'.format(algo, i)), 'w')
        pickle.dump(d, f); f.close()
        
    return

def writeJobArray(algo, paramFolder = '../resultData/nucleus/parameters'):
    jobNum= len(filter(lambda x: algo in x, os.listdir(paramFolder)))

    progFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc'
    scriptFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/scripts'
    path_command = """setenv PATH /cbio/donnees/nvaroquaux/.local/bin:${PATH}
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/cbio/donnees/nvaroquaux/.local/lib
    setenv LIBRARY_PATH /cbio/donnees/nvaroquaux/.local/lib
    setenv PYTHONPATH /cbio/donnees/aschoenauer/workspace2/cecog/pysrc:/cbio/donnees/aschoenauer/workspace2/Tracking/src
    setenv DRMAA_LIBRARY_PATH /opt/gridengine/lib/lx26-amd64/libdrmaa.so
    """
    pbsOutDir = '/cbio/donnees/aschoenauer/PBS/OUT'
    pbsErrDir = '/cbio/donnees/aschoenauer/PBS/ERR'
    pbsArrayEnvVar = 'SGE_TASK_ID'
    
    baseName = algo
    
    for k in range(jobNum):
        head = """#!/bin/sh
cd %s""" %progFolder

        cmd = ''
            # command to be executed on the cluster
        temp_cmd = """
python nucleus/clustering.py -a %s -o %i
"""
        temp_cmd %= (
                     algo,
                     k
                    )

        cmd += temp_cmd
        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(k+1))
        script_file = file(script_name, "w")
        script_file.write(head + cmd)
        script_file.close()

        # make the script executable (without this, the cluster node cannot call it)
        os.system('chmod a+x %s' % script_name)

        # write the main script
    array_script_name = '%s.sh' % os.path.join(scriptFolder, baseName)
    main_script_file = file(array_script_name, 'w')
    main_content = """#!/bin/sh
%s
#$ -o %s
#$ -e %s
%s$%s.sh
""" % (path_command,
       pbsOutDir,  
       pbsErrDir, 
       os.path.join(scriptFolder, baseName),
       pbsArrayEnvVar)

    main_script_file.write(main_content)
    main_script_file.close()
    os.system('chmod a+x %s' % array_script_name)

    # the submission commando is:
    #sub_cmd = 'qsub -o %s -e %s -t 1-%i %s' % (self.oBatchSettings.pbsOutDir,  
    #                                           self.oBatchSettings.pbsErrDir, 
    #                                           jobCount, array_script_name)
    sub_cmd = 'qsub -t 1-%i %s' % (jobNum, array_script_name)
    print sub_cmd
    return

def plotMoy(outFolder='../resultData/nucleus', show=False):
    f=p.figure(figsize=(24,17))
    
    for i, algo in enumerate(["gaussianmixture", 'mini-batchk-means', 'kernelk-means','spectralclustering']):
        file_=open(os.path.join(outFolder, outFile.format(algo)), 'r')
        resultDict = pickle.load(file_); file_.close()
    
        perClass = []; perCluster = []
        for parameter_set in resultDict:
            currLperClass=[]
            currLperCluster = []
            for morphoClass in resultDict[parameter_set]:
                currLperClass.append(np.mean(resultDict[parameter_set][morphoClass][0]))
                currLperCluster.append(np.mean(resultDict[parameter_set][morphoClass][1]))
            perClass.append(currLperClass)
            perCluster.append(currLperCluster)
            
        
        ax=f.add_subplot(2,2,i)
        for k in range(len(perClass)):
            ax.errorbar(np.mean(perCluster[k]), np.mean(perClass[k]), xerr = np.std(perCluster[k]), yerr=np.std(perClass[k]), 
                            fmt = markers[methods.index(algo)],
                           label = '{}'.format(algo))
        ax.grid(True)
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_title("Algorithm {}, on average on all morphological classes".format(algo))
        ax.set_xlabel("Nb of class elements in cluster/cluster size")
        ax.set_ylabel("Nb of class elements in cluster/class size")
    if show:
        p.show()
    else:
        p.savefig(os.path.join(outFolder, 'figure {}_MOY.png'.format(algo)))


class NucleusClustering():
    
    def __init__(self, inFolder, filename, algo,pca, whiten, outFolder='../resultData/nucleus',n_iteration =10, **kwargs):
        '''
        pca : None pour pas de PCA, nombre de composantes principales a garder sinon
        dans kwargs je mets les parametres de la methode de clustering
        n'oublions pas qu'il n'y a que 3,000 elements dans le training set donc tout va bien :))
        '''
        assert algo in methods, 'Given algorithm is not implemented'
        #assert sorted(params[algo]) == sorted(kwargs.keys()), 'Parameters for given algorithm are not given'
        self.file=os.path.join(inFolder, filename)
        self.algo = algo
        self.iterations= n_iteration
        self.outFolder = outFolder
        self.outFile = outFile.format(self.algo)
        
        self.parameters = kwargs
        self.pca = pca
        self.whiten= whiten
        
        try:
            self.k = self.parameters["n_clusters"]
        except KeyError:
            self.k = self.parameters["n_components"]
        
    def cluster(self,data):
        if self.algo == 'mini-batchk-means':
            return self.batch_kmeans(data)
        elif self.algo=='kernelk-means':
            return self.kernel_kmeans(data)
        elif self.algo=='c-means':
            return self.c_means(data)
        elif self.algo=='gaussianmixture':
            return self.gmm(data)
        elif self.algo=='isomapping':
            return self.isomap(data)
        elif self.algo=='spectralclustering':
            return self.spectral_clustering(data)
        
    def batch_kmeans(self, data, parameters = None):
        if parameters == None:
            parameters = self.parameters
        model = MiniBatchKMeans(**parameters)
        labels = model.fit_predict(data)
        
#        print 'Computing silhouette score'
#        silhouette = silhouette_score(data, labels, metric='euclidean')
#        
#        print 'Computing cohesion'
#        total_sum_squares = np.sum((data - np.mean(data, 0))**2)
#        centers = model.cluster_centers_
#        cohesion=0
#        for p in range(parameters['n_clusters']):
#            cohesion+=np.sum((data[np.where(labels == p)]-centers[p])**2)
#        cohesion = cohesion/float(total_sum_squares)
        
        return labels
    
    def kernel_kmeans(self,data):
        model = KernelKMeans(**self.parameters)

        labels = model.fit_predict(data)
        final_dist = model.final_dist
        
        print 'Computing silhouette score'
    #TODO on peut calculer le silhouette score mais il faut redefinir la distance qui est utilisee 
        #silhouette = silhouette(data, labels, N)
        
        print 'Computing cohesion'
        total_sum_squares = data.shape[0]*np.trace(model.gram_matrix)-np.sum(model.gram_matrix)
        cohesion=np.sum((final_dist)**2)/float(total_sum_squares)
        
        return labels
    
    def gmm(self, data):
        #partition-entropy coefficient, see Validity_survey p. 31
        #Fukuyama-Sugeno index, see same article p.33
        
        model = GMM(**self.parameters)#, thresh=0.01, min_covar=0.001, params='wmc', init_params='wmc')
        
        model.fit(data)
        prob_labels = model.predict_proba(data)+10**(-15)
        labels = model.predict(data)
        
        score = model.score(data)
        BIC = model.bic(data) 
        FS=0
        for p in range(self.k):
            zz = data[np.where(labels==p)]
            if len(zz)==0:
                continue
            center = np.mean(zz, 0)
            FS+=np.sum(prob_labels[:, p]*(np.sum((data-center)**2, 1)- np.linalg.norm(center)**2))
               
            #Partition entropy coefficient normalized by its max (not done for real and RR data)
        if self.k>1:
            partition_entropy=-1/float(data.shape[0]*np.log2(self.k))* np.sum(prob_labels * np.log2(prob_labels))
            
        return labels
    
    def spectral_clustering(self, data):
        model = SpectralClustering(**self.parameters)
        labels = model.fit_predict(data)
    #possibilite de calculer cohesion et silhouette score si on utilise kmeans. Discretize kezako ?
        return labels
    
    
    def homemade_spectral_clustering(self, data):
        '''
        Affinity(x,y)=exp(-sigma*||x-y||^2)
        there is a pbl smw I think
        '''
    
        print "Computing KDTree on data"
        tree = KDTree(data)
        W = np.zeros(shape=(data.shape[0], data.shape[0]), dtype=np.float)
        ## Put the values in the matrix according to k-NN approach + rbf similarity measure
        #ms si j'ai bien compris on peut aussi faire juste l'un ou juste l'autre, il faut regarder
 
        print 'Computing affinity matrix'
        for n in range(data.shape[0]):
            for k in range(n+1):
                #now using complete affinity matrix instead of kNN
                W[n, k]=dist(data[k], data[n], self.parameters["sigma"])
                W[k,n]=W[n,k]
#            dou, i=tree.query(data[n], k=neighbours_nb)
#            for couple in filter(lambda x: x[1]<np.inf, zip(i,dou)):
#                W[n,couple[0]] = dist(data[n], data[couple[0]], sig)
#                #we use kNN approach and not mutual kNN approach
#                W[couple[0],n] = W[n,couple[0]]

        Lrw = homeMadeGraphLaplacian(W, normed=True)
            
        if np.any(np.diagonal(Lrw)==0):
            raise
#            outliersGraves = np.where(np.diagonal(Lrw)==0)
#            dataBis = np.delete(data, outliersGraves, 0)
#            Wbis = np.delete(np.delete(W, outliersGraves, 0), outliersGraves, 1)
#            f=open("outliersGraves{}.pkl".format(time.clock()), 'w'); pickle.dump(data[outliersGraves],f); f.close()
#            print "Important outliers that are only connected to themselves, I take them off"
#            return homeMadeSpectralClust(dataBis, cluster_nb_min, cluster_nb_max, neighbours_nb, sig, show=show, affinity = Wbis)
    
        print 'Diagonalizing normalized Laplacian'; Lrw=csc_matrix(Lrw)
        vals_rw, vecs_rw = eigsh(Lrw, k=self.k+10, which = 'LM', tol=0.01, sigma=0, return_eigenvectors=True)
        
        #chaque ligne de u est donc l'embedding d'une ligne de data. Maintenant on peut faire KMeans dessus. On peut aussi plotter sur les deux premieres composantes
        #pour voir ce que ca donne
        U = vecs_rw[:,:self.k]
        eigengap=vals_rw[self.k] - vals_rw[self.k-1]
        
#        #let's look at the Gaussian mixture result
#        BIC, FS, score, partition_entropy = GaussianMixture(U, cluster_nb, cluster_nb)
#        BIC_r.extend(BIC); FS_r.extend(FS); PEC_r.extend(partition_entropy)
        
#TODO relire spectral clustering, is new_data == U ou U[:,:self.k] ??
        return self.batch_kmeans(U, **dict(zip(params["mini-batch k-means"], [17, 1000, 500, 1000, 'k-means++', 5]))) 
            
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
    
    def isomap(self, data):
        
        print 'Isomap neighbours :', self.parameters["n_neighbors"]
        print 'Isomap components, ie final number of coordinates :', self.k
        m = Isomap(**self.parameters)#eigen_solver='auto', tol=0, path_method='auto', neighbors_algorithm='kd_tree')
        x = m.fit_transform(data)
        error=m.reconstruction_error() 
        geod_d = m.dist_matrix_.flatten()
        new_euclid_d = cdist(x, x, metric='euclidean').flatten()
        corr=1- pearsonr(geod_d, new_euclid_d)[0]**2
        
        new_data = x
        return self.batch_kmeans(new_data, parameters = dict(zip(params["mini-batch k-means"], [17, 1000, 500, 1000, 'k-means++', 5])))
    
    def c_means(self, data):
        k=self.parameters['n_clusters']
        m=self.parameters["fuzzyfier"]
        
        XB_courant=[]; FS_courant=[]; partition_entropy_courant = []
        labels = []
        for iteration in range(self.parameters["n_iter"]):
            print iteration,
            mu=np.zeros(shape=(data.shape[0], k))
            #initialization: matrix Mu must be a stochastic matrix with no class with everyone or no one
            
            mu[:,0:int(k/2)]=0.5*np.random.rand(data.shape[0], int(k/2))
            if k%2==0: mu[:,k/2:]=0.5-mu[:,0:k/2]
            else:
                mu[:,int(k/2):-1]=0.65*(0.5-mu[:,0:int(k/2)])
                mu[:,-1]=1-np.sum(mu[:,:-1], 1)
            np.random.shuffle(mu); np.random.shuffle(mu.transpose())

            model = FuzzyCMeans(data, mu, m)
            centers = model(emax=1.e-10, imax=self.parameters["max_iter"])
            prob_labels = model.mu+10**(-15); labels.append(prob_labels)
            
            tree=KDTree(centers)
            min_center_distance = min([max(tree.query(centers[l], 2)[0]) for l in range(k)])
            
            FS=0; XB = 0
            for p in range(k):
                center = centers[p]
                FS+=np.sum((prob_labels[:, p]**m)*(np.sum((data-center)**2, 1)- np.linalg.norm(center)**2))
                XB+=np.sum((prob_labels[:, p]**2)*np.sum((data-center)**2, 1))
            XB/=(data.shape[0]*min_center_distance)
            XB_courant.append(XB); FS_courant.append(FS) 
            partition_entropy_courant.append(-1/float(data.shape[0]*np.log2(k))* np.sum(prob_labels * np.log2(prob_labels)))
            
        #sortie de la boucle d'iterations et choix des meilleurs indicateurs
        print np.argmin(FS_courant), np.argmin(XB_courant), np.argmin(partition_entropy_courant)
        print partition_entropy_courant
        return labels[np.argmin(XB_courant)]
        
        
    def __call__(self):
        self.classNames, data, trueLabels = getData(self.file, pca = self.pca, whiten =self.whiten)
        result = {el:([], []) for el in self.classNames}
        for iter_ in range(self.iterations):
            try:
                predictedLabels = self.cluster(data)
            except:
                sys.stderr.write('Problem iter {} algo {}, pca {}, whiten {}, params {}'.format(iter_, self.algo,self.pca, self.whiten, self.parameters))
                continue
        
    #        print Counter(predictedLabels)
            if len(predictedLabels.shape)>1:
                return predictedLabels
            if True:#self.k == len(self.classNames):        
                confrontation = np.array(zip(trueLabels, predictedLabels))
                confusion_matrix = None
                for el in self.classNames:
                    arr = np.array([len(np.where(confrontation[np.where(confrontation[:,0]==el)][:,1]==str(j))[0]) for j in range(self.k)])
                    confusion_matrix = arr if confusion_matrix is None else np.vstack((confusion_matrix, arr))
    
                confusion_matrix = np.array(confusion_matrix)
                for i, el in enumerate(self.classNames):
                    print '{:<30}{}'.format( el, np.argmax(confusion_matrix[i]))
                    result[el][0].append(
                                np.max(confusion_matrix[self.classNames.index(el)]/float(np.sum(confusion_matrix[self.classNames.index(el)]))))
                    result[el][1].append(
                                np.max(confusion_matrix[self.classNames.index(el)])/float(np.sum(confusion_matrix[:,np.argmax(confusion_matrix[self.classNames.index(el)])]))
                                )
                
                
        self.parameters['pca'] = self.pca
        self.parameters['whiten']= self.whiten
        
        result = {tuple(zip(self.parameters.keys(), self.parameters.values())):result}
        
        if self.outFile in os.listdir(self.outFolder):
            f=open(os.path.join(self.outFolder, self.outFile), 'r')
            d=pickle.load(f); f.close()
            d.update(result)
        else:
            d=result
        f=open(os.path.join(self.outFolder, self.outFile), 'w')
        pickle.dump(d, f); f.close()
        self.plot(d, show)
        return
        
    def plot(self, resultDict=None, show=False):
        if resultDict == None:
            f=open(os.path.join(self.outFolder, self.outFile))
            resultDict= pickle.load(f); f.close()
            
        if type(resultDict.keys()[0])==str:
            #meaning we are dealing with many alg results
            pass
        elif type(resultDict.keys()[0])==tuple:
            f=p.figure(figsize=(24,17))
            ax=f.add_subplot(111)
            colors = makeColorRamp(17, hex_output=True)
            for param_set in resultDict:
                for morpho in resultDict[param_set]:
                    r = resultDict[param_set][morpho]
                    ax.errorbar(np.mean(r[1]), np.mean(r[0]), xerr = np.std(r[1]), yerr=np.std(r[0]), 
                            color =colors[self.classNames.index(morpho)], fmt = markers[methods.index(self.algo)],
                           label = '{} {}'.format(self.algo, morpho))
            ax.grid(True)
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)
            ax.set_title("Algorithm {}".format(self.algo))
            ax.set_xlabel("Nb of class elements in cluster/cluster size")
            ax.set_ylabel("Nb of class elements in cluster/class size")
            ax.legend()
            if show:
                p.show()
            else:
                p.savefig(os.path.join(self.outFolder, 'figure {}.png'.format(self.algo)))
                                
                
if __name__ == '__main__':
    description =\
'''
%prog - Nucleus clustering
'''
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)

    parser.add_option("-a", "--algo", dest="algo", type=str,
                      help="Clustering algorithm")
    parser.add_option('-o', '--o', dest = 'i', type=int, help = 'Option file')
    
    (options, args) = parser.parse_args()
    inFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/resultData/nucleus/parameters'
    f=open(os.path.join(inFolder, 'parameters_{}_{}.pkl'.format(options.algo, options.i)))
    clustering_parameters = pickle.load(f); f.close()
    
    print clustering_parameters
    
    for el in product(pca_param, whiten_param):
        print 'pca, whiten', el
        n = NucleusClustering('../data/challenge_clustering/MitoCheck_10x/data', 'all_features.arff', options.algo, 
                              el[0], el[1],  **clustering_parameters)
        n()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
