import getpass
if getpass.getuser()=='aschoenauer':
    mitocheck_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/mitocheck_siRNAs_target_genes_Ens75.txt'
    qc_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/qc_export.txt'
    
    #folders
    outputFolder = '../resultData/cell_cycle/galeries'
    hdf5Folder = '/share/data20T/mitocheck/Alice/results'
    rawDataFolder='/share/data20T/mitocheck/compressed_data'
    trackingFolder = '/share/data20T/mitocheck/tracking_results'
    
    new_h5=False
    
    trackingFilename='traj_noF_densities_w{}.hdf5.pkl'
    imageFilename = "--W{}--P00001--T{:>05}-"
    min_=1
    max_=2**8
    
if getpass.getuser()=='lalil0u':
    mitocheck_file = '/media/lalil0u/New/workspace2/Xb_screen/data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt'
    qc_file = '/media/lalil0u/New/workspace2/Xb_screen/data/mapping_2014/qc_export.txt'
    
    #folders
    outputFolder = '../resultData/cell_cycle/galeries'
    hdf5Folder = '/media/lalil0u/New/projects/H2B_PCNA'
    rawDataFolder='/media/lalil0u/New/data/H2B_PCNA'
    trackingFolder = '/media/lalil0u/New/projects/H2B_PCNA/results'
    new_h5=True
    exist_QC=False
    
    trackingFilename='traj_intQC_w{}.ch5.pkl'
    imageFilename = "__P{}__T{:>05}__Cgfp"
    
    min_=60
    max_=800

##filename for output image
outputImage="{}crop_P{}_W{}_{}_t{}_id{}.png"

#modelFilename = "best_model_split.pkl"
#outputTrainingFilename = "thrivisions_featureMatrix_training.pkl"
#outputTestingFilename="thrivisions_featureMatrix_testing.pkl"
#outputPredictingFilename = "thripred_{}_{}.pkl"

#Say if you also want trajectories that start with a mitosis and end with the end of the movie. Put to True for that
not_ending_track=True
if not_ending_track:
    outputFile="cell_cycle_cens_{}.pkl"
else:
    outputFile="cell_cycle_info_{}.pkl"

#When doing galeries, say if you want to recompute the complete/incomplete tracks
redo=False

margin=20

renumber=True
XMAX = 1344 
YMAX = 1024

#Objective to get total intensity
#objective = {'name':'total intensity','function':lambda x:np.multiply(x['roisize'], x['n2_avg']), 'features':['roisize', 'n2_avg']}
#Objective to get nucleus size
objective='roisize'