from tracking.trajPack import featuresNumeriques

#settings for aschoenauer@cbio

#folder for trajectory features data
data_folder = '/share/data20T/mitocheck/tracking_results'

#folder for output data
result_folder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/resultData/summaries'
summary_filename = "summary_{}_{}.pkl"
clustering_filename = 'centers.pkl'
hit_filename = 'hits.pkl'

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