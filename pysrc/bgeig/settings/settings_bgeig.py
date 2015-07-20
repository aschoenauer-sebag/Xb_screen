import getpass, os

outputFolder = "/media/lalil0u/New/projects/Geiger/results"
dataFolder="/media/lalil0u/New/projects/Geiger/nuclei"
    
loadingFolder='../prediction'
    
intensity_qc_file=None
traj_filename = 'traj_intQC_w{}.pkl'
all_trajectory_xml=True

feature_filename='features_intQC_{}.pkl'

training = False
predict = True
repeat=True
filtering_fusion=False

time_windows=[(0, 48), (0,24), (12, 36), (24,48),(36, 60)] 

separating_function=lambda x: int(x.split('_')[1])

removeTempFiles=True

redo=False