import cPickle as pickle
import numpy as np
import os, pdb
#import matplotlib
#matplotlib.use('Agg')
import fileManagement
from trajPack import clustering
'''
File for scripts when opening ipython
'''

#
#f=open('../data/lib123.pkl'); l=pickle.load(f); f.close()
#done = fileManagement.countingDone(l)
#
#means=[]
#f=open('../means_ctrl.pkl'); means.append(pickle.load(f)); f.close()
#
#for k in range(10):
#    print k
#    r, _,_,_,_,_,_,_=clustering.histConcatenation('/share/data20T/mitocheck/tracking_results/', done[1000*k:1000*(k+1)], '../data/mitocheck_siRNAs_target_genes_Ens72.txt', '../data/qc_export.txt')
#    means.append((r.shape, np.mean(r,0)))
#    
#f=open('../all_means.pkl', 'w')
#pickle.dump(means, f); f.close()


#from sandbox import histLogTrsforming
#f=open('../resultData/simulated_traj/hist_tabFeatures_WO.pkl')
#tabFeatures, _, histNC, _, _, _ = pickle.load(f)
#r2 = histLogTrsforming(tabFeatures)  
#print r2.shape, 'not normalized'
#r2=r2[:,:-1]; r2=(r2-np.mean(r2,0))/np.std(r2,0)
#print 'computing bins'
#histogrammeMatrix = transportation.computingBins(histNC, mat_hist_sizes[0])
#print 'saving'
#f=open('../resultData/simulated_traj/histogramsNtotSim.pkl', 'w')
#pickle.dump(np.hstack((r2, histogrammeMatrix)), f); f.close()
