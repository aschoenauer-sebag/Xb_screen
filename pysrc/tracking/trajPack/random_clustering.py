from optparse import OptionParser
import os,sys
import cPickle as pickle
import numpy as np

from clustering import BatchKmeans

if __name__ == '__main__':
    parser = OptionParser(usage="usage: %prog [options]",
                         description='description')
    
    parser.add_option("-f", "--folder", dest="folder", default='../resultData/ECCB',
                      help="The folder where the data is and everything saved")
    
    parser.add_option("-a", "--algo", dest="algo", default = 3, 
                      help="0 for Isomap\
                       1 for Gaussian mixtures\
                       2 for K-means\
                       3 for mini-batch K-means\
                       4 for fuzzy C-means\
                       5 for IBP")
    parser.add_option("-n", "--n_iter", dest="n_iter", default = 0, 
                      help="Nb de fois")
    
    parser.add_option("-d", "--data", dest="data", default = 0, 
                      help="Quelles donnees")
    
    (options, args) = parser.parse_args()

#    folder = '/cbio/donnees/aschoenauer/data/tracking/migration'
# PUIS features with densities: folder = '/share/data/mitocheck/tracking_results'
#    allData, whoIsInData, length = concatenation(folder)
    #names = dict(zip(range(7), ['isomap', 'gaussianM', 'kMeans', 'batchKMeans', 'SpecClust', 'cMeans', 'IBP']))
    baseNames = dict(zip(range(4), ['', 'RR_', 'RG_', 'real_'])); baseName = baseNames[int(options.data)]
    f_size = 7
    s = 1000000    

    if int(options.data)==1:
        print "Working with random data (uniform), shape {} x {}".format(s,f_size)
        f=open(os.path.join(options.folder, 'minMax.pkl'), 'r')
        mini, maxi = pickle.load(f); f.close()

        r = np.random.random_sample(size=(s, f_size))*(maxi-mini)+mini
        r = (r-np.mean(r, 0))/np.std(r, 0)
        
    elif int(options.data)==2:
        print "Working with random data (Gaussian), noise level {}, shape {} x {}".format(0,s,f_size)
        mean = np.zeros(shape=(f_size,))
        stdev = np.ones(shape=(f_size,))
#        f=open(os.path.join(options.folder, 'minMax.pkl'), 'r')
#        mini, maxi = pickle.load(f); f.close()
        
        currData = np.random.normal(mean, stdev, size=(s, f_size))
        r=(currData-np.mean(currData, 0))/np.std(currData, 0)
        
#    elif int(options.data)==3:
#        print "Working with all real data, only {} features, max nb of neighbours {}".format(f_size, options.dens)
#        f=open(os.path.join(options.folder, 'small_l123.pkl'), 'r')#'trsf_maxd50inf{}_data.pkl'.format(options.dens)), 'r')
#        zou= pickle.load(f); f.close()
#        r=zou-np.mean(zou, 0)
#        #print mot
#        #r=r[:,:-4]
#        if int(options.algo) not in [2,3]:
#            np.random.shuffle(r); currData = r[:15000]
#            dataCtrl = concatCtrl(r, ctrlStatus, length, len(ctrlStatus)-np.sum(ctrlStatus))
#            np.random.shuffle(dataCtrl)
#            ctrl = dataCtrl[:5000]################
#            dPheno = concatCtrl(r, 1-np.array(ctrlStatus), length, np.sum(ctrlStatus))
#    #pas fait pour les premieres experiences sur trsformed data
#            np.random.shuffle(dPheno)
#            
#            currData = np.vstack((ctrl, dPheno[:10000]))###############
        
    else:
        print "Settings not understood"
        sys.exit()
        
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
    if int(options.algo) ==3:
        print 'BatchKMeans'
        silhouette_r, cohesion_r = BatchKmeans(r, 2, 30, N=10)
        f=open(os.path.join(options.folder, baseName+'batchKMeans{}.pkl'.format(options.n_iter)), 'w')
        pickle.dump([silhouette_r, cohesion_r], f); f.close()
        
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