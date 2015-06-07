
mitocheck_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/mitocheck_siRNAs_target_genes_Ens75.txt'
qc_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/qc_export.txt'

#folders
outputFolder = '../resultData/cell_cycle/galeries'
hdf5Folder = '/share/data20T/mitocheck/Alice/results'
rawDataFolder='/share/data20T/mitocheck/compressed_data'
trackingFolder = '/share/data20T/mitocheck/tracking_results'
#mitocheck_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/mitocheck_siRNAs_target_genes_Ens75.txt'
#quality_control_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/qc_export.txt'

##filenames
imageFilename = "--W{}--P00001--T{:>05}-"
trackingFilename='traj_noF_densities_w{}.hdf5.pkl'
outputImage="{}crop_P{}_W{}_{}_t{}_id{}.png"
outputFile="cell_cycle_cens_{}.pkl"
#modelFilename = "best_model_split.pkl"
#outputTrainingFilename = "thrivisions_featureMatrix_training.pkl"
#outputTestingFilename="thrivisions_featureMatrix_testing.pkl"
#outputPredictingFilename = "thripred_{}_{}.pkl"

#Say if you also want trajectories that start with a mitosis and end with the end of the movie
not_ending_track=False

min_=1
max_=2**8
margin=20
new_h5=False
renumber=True
XMAX = 1344 
YMAX = 1024

#Objective to get total intensity
#objective = {'name':'total intensity','function':lambda x:np.multiply(x['roisize'], x['n2_avg']), 'features':['roisize', 'n2_avg']}
#Objective to get nucleus size
objective='roisize'