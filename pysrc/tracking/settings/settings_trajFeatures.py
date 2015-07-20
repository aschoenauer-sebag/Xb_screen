import getpass, os

if getpass.getuser()=='aschoenauer':
    outputFolder = "/share/data20T/mitocheck/tracking_results"
elif getpass.getuser()=='lalil0u':
    outputFolder = "/media/lalil0u/New/projects/H2B_PCNA/results"
    dataFolder="/media/lalil0u/New/projects/H2B_PCNA"
    
loadingFolder='../prediction'

intensity_qc_file=None
traj_filename = 'traj_intQC_w{}.pkl'
all_trajectory_xml=False

feature_filename='features_intQC_{}.pkl'

training = False
predict = True
repeat=True
filtering_fusion=False

time_windows=[(0, 48), (0,24), (12, 36), (24,48),(36, 60)] 

separating_function=lambda x: x.split('_')[1]

removeTempFiles=True

redo=False