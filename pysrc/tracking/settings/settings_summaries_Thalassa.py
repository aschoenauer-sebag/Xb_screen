from tracking.trajPack import featuresNumeriques

#settings for aschoenauer@cbio

#folder for trajectory features data
data_folder = '/share/data20T/mitocheck/tracking_results'

#folder for output data
result_folder = '/share/data20T/mitocheck/Alice/summaries'
summary_filename = "summary_iter{}_{}_{}.pkl"
clustering_filename = 'centers_iter{}.pkl'
hit_filename = 'hits_KS.pkl'
pval_filename = 'pval_L_iter{}_{}.pkl'
label_filename = 'labels_iter{}_{}.pkl'
figure_name = "comparaison_{}_{}.png" 
threshold_filename = 'thresholds_significancy_5_{}.pkl' 

mitocheck_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/mitocheck_siRNAs_target_genes_Ens72.txt'

quality_control_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/qc_export.txt'

#fraction on which to check clustering stability
fraction = 0.8
num_iterations_stability = 10

#number of clusters to try
k_min = 3
k_max = 11 

#number of trajectories to represent a cluster
n_representatives = 10

#number of numeric features
nb_feat_num = len(featuresNumeriques)

#list to take into account or not pca, and whitening
pcaParameters = [(0,0), (1,0), (1,1)]
#number of principal components that are retained to get approx 92% of the variance
nb_composantes = 6

#Fisher test or Kolmogorv-Smirnov ?
Fisher = True

#do you want to calculate all p-values from all parameter sets even if they have already been calculated?
redo = False