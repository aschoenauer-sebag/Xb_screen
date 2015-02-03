
qc_file='/media/lalil0u/New/workspace2/Xb_screen/data/xb_manual_qc.pkl'
qc_file2='/media/lalil0u/New/workspace2/Xb_screen/data/xb_focus_qc.pkl'

#folders
outputFolder = '/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions__settings2/thrivisions'
rawDataFolder='/media/lalil0u/TOSHIBA EXT/images_A/'
hdf5Folder = '/media/lalil0u/New/projects/Xb_screen/plates__all_features_2bis'
trackingFolder = '/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions__settings2/'
#mitocheck_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/mitocheck_siRNAs_target_genes_Ens75.txt'
#quality_control_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/qc_export.txt'

##filenames
imageFilename = "--W{}--P0001_t{:>05}_c00002"
trackingFilename='traj_intQC_{}.hdf5.pkl'
trackingFilename2='traj_intQC_w{}.hdf5.pkl'
outputImage="crop_P{}_W{}_t{}_{}_id{}.png"
outputFile="thrivisions_{}_{}.pkl"

min_=50
max_=1000
margin=20
new_h5=True
renumber=False
XMAX=1920
YMAX=1440