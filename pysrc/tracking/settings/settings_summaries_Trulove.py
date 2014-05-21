from tracking.trajPack import featuresNumeriques

#settings for lalil0u@Trulove

#folder for trajectory features data
data_folder = '/media/lalil0u/New/data/test'

#folder for output data
result_folder = '/media/lalil0u/New/workspace2/Xb_screen/resultData/summaries'
summary_filename = "summary_{}_{}.pkl"
clustering_filename = 'centers.pkl'
hit_filename = 'hits_KS.pkl'
pval_filename = 'pval_KS_{}.pkl'
figure_name = "comparaison_p-values.png" 

mitocheck_file = '/media/lalil0u/New/workspace2/Xb_screen/data/mapping/mitocheck_siRNAs_target_genes_Ens72.txt'

quality_control_file = '/media/lalil0u/New/workspace2/Xb_screen/data/mapping/qc_export.txt'

#fraction on which to check clustering stability
fraction = 0.8
num_iterations_stability = 6

#number of clusters to try
k_min = 3
k_max = 11 

#number of trajectories to represent a cluster
n_representatives = 10

#number of numeric features
nb_feat_num = len(featuresNumeriques)