import getpass, os

if getpass.getuser()=='aschoenauer':
    outputFolder = "/share/data20T/mitocheck/tracking_results"
elif getpass.getuser()=='lalil0u':
    outputFolder = "/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions__settings2"
    dataFolder="/media/lalil0u/New/projects/Xb_screen/plates__all_features_2bis"
    
traj_filename = 'traj_intQC_w{}.pkl'
feature_filename='features_intQC_{}.pkl'

training = False
predict = True
repeat=True
filtering_fusion=False

time_windows=[(0, 48), (0,24), (12, 36), (24,48),(36, 60)] 