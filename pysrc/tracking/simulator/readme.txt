To apply the workflow to simulated trajectories
i. Produce the simulated trajectories AND the extra controls that are needed
		p=simulator.PlateSimulator(settings_filename='tracking/settings/settings_simulator_14_10_20.py')
		p()
		
ii. Extract the features of those trajectories
		%run tracking/trajPack/features_script -b simFeat -d ../resultData/simulated_traj/simres/plates/ --simulated 1
		
iii. If never done before, generate the false Quality Control file and change the name of feature files
		simulator.generateQC()
		
iv. Do the jobs for the simulated traj, controls first
		%run tracking/trajPack/feature_cell_extraction_script -b simDistances --div_name KS --testCtrl 1 --simulated 1  
		
v. For the experiments
	for k in range(5):
		%run tracking/trajPack/feature_cell_extraction_script -b simDistances --div_name KS --simulated 1 --iter $k --simulated 1 --siRNA ../data/siRNA_simulated.pkl
		
vi. Extract collected distances
	r=feature_cell_extraction.collectingDistances('dist_sim_CTRL.pkl', 
		'../resultData/simulated_traj/simres', qc_filename='../data/qc_simulated.txt',mapping_filename= None, testCtrl=True, redo=True,long_version=False, key_name='dist_sim')
		
	r=feature_cell_extraction.collectingDistances('dist_sim.pkl', 
		'../resultData/simulated_traj/simres', qc_filename='../data/qc_simulated.txt',mapping_filename= None, testCtrl=False, redo=True,long_version=False, key_name='dist_sim')
vii. Compute hits
	r=feature_cell_extraction.multipleHitDistances('../resultData/simulated_traj/simres','dist_sim', 
	qc_filename='../data/qc_simulated.txt', mapping_filename=None, filename='all_dist_sim', combination='max', redo=False, trad=False, without_mean_persistence=True)
	
viii. Compare with planned hits

		
		
To apply the workflow to real trajectories
i. Check if there's hdf5/traj/features left to compute
	
	f=open('../data/expL_targeted_Mitocheck_2014.pkl')
	import cPickle as pickle
	l=pickle.load(f)
	f.close()
	from util.listFileManagement import countingDone, appendingControl, expSi, strToTuple
	ctrl=appendingControl(l)
	l.extend(ctrl)
	raw,tracking,feat=countingDone(l, featlistonly=False, featureTabName='hist_tabFeatures_{}.pkl')
	
	** HDF5 left
		f=open('../mitocheck_hdf5_left.pkl', 'w')
		pickle.dump(raw,f)
		f.close()
		
		[cd ~/workspace2/cecog/pysrc/scripts/EMBL/cluster]
		python pbs_script_generation.py -b correction_mitocheck -p ../mitocheck_hdf5_left.pkl 
	
	** FEATURES left	
		f=open('../mitocheck_feat_left.pkl', 'w')
		pickle.dump(feat,f)
		f.close()
		
		python tracking/trajPack/features_script.py -b correction_feat --feat ../mitocheck_feat_left.pkl
		
ii. Workflow jobs, controls first
		%run tracking/trajPack/feature_cell_extraction_script -b wholeDistances --div_name KS --testCtrl 1  
		
iii. For the experiments
	for k in range(5):
		%run tracking/trajPack/feature_cell_extraction_script -b wholeDistances --div_name KS --iter $k --siRNA ../data/siRNA_targeted_Mitocheck_2014.pkl



