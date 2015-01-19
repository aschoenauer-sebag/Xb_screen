result_folder = '../../../projects/Xb_screen/dry_lab_results/track_predictions__settings2/features_on_films'
outputFile = 'distances_tw_{}.pkl'
ctrl_exp_filename = 'ctrl_exp_{}_{}.pkl'
data_folder = '../../../projects/Xb_screen/dry_lab_results/track_predictions__settings2/'
visual_qc = '../data/xb_manual_qc.pkl'
flou_qc='../data/xb_focus_qc.pkl'

##filename for trajectory features
filename='features_intQC_{}_01.pkl'
filename_twindow="features_intQC_{{}}_01_t{time_window}.pkl"

##filename format for the wells
#format_='{:>05}_01'

#Option to say that we're only dealing with numerical features for the moment
histDataAsWell=False

#Option to say if we want to redo experiments anyway
redo=True
#to tell it's xb screen data
xb_screen=True