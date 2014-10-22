To apply the workflow to simulated trajectories
i. Produce the simulated trajectories
		p=simulator.PlateSimulator(settings_filename='tracking/simulator/simulation_settings_2014_01_16.py')
		p()
		
	AND the extra controls that are needed
		
ii. Extract the features of those trajectories
		%run tracking/trajPack/features_script -b simFeat -d ../resultData/simulated_traj/simres/plates/ --simulated 1
		
iii. If never done before, generate the false Quality Control file
		simulator.generateQC()
		
iv. Do the jobs for the simulated traj
		%run tracking/trajPack/simulated
