import getpass, os

if getpass.getuser()=='aschoenauer':
    outputFolder = "/cbio/donnees/aschoenauer/projects/Xb_screen/dry_lab_results/track_predictions__settings2"
    dataFolder="/cbio/donnees/aschoenauer/Xb_screen/plates__all_features_2bis"
elif getpass.getuser()=='lalil0u':
    outputFolder = "/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions__settings2"
    dataFolder="/media/lalil0u/New/projects/Xb_screen/plates__all_features_2bis"
    
traj_filename = 'traj_intQC_w{}.pkl'
feature_filename='features_intQC_{}_t{}.pkl'

predict = True
repeat=False
filtering_fusion=False

#Indicating in hours the time windows under study
time_windows=[(0, 48), (0,24), (12, 36), (24,48),(36, 60)] 
#Indicating the number of frames per hour
time_lapse=4