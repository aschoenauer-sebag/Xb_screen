result_folder = '../../../projects/Xb_screen/dry_lab_results/track_predictions__settings2/features_on_films'
outputFile = 'distances_XbScr_5CtrlC_{}.pkl'
ctrl_exp_filename = 'ctrl_exp_{}_{}.pkl'
data_folder = '../../../projects/Xb_screen/dry_lab_results/'
visual_qc = '../data/xb_manual_qc.pkl'
flou_qc='../data/xb_focus_qc.pkl'

##filename for trajectory features
filename='features_intQC_{}_01.pkl'

##filename format for the wells
#format_='{:>05}_01'

#Option to say that we're only dealing with numerical features for the moment
histDataAsWell=False

#Option to say if we want to redo experiments anyway
redo=True
#to tell it's not the xb screen
xb_screen=True