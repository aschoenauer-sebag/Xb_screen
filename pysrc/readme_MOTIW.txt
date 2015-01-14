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
	r=feature_cell_extraction.collectingDistances('all_dist_sim_CTRL.pkl', 
		'../resultData/simulated_traj/simres', qc_filename='../data/qc_simulated.txt',mapping_filename= None, testCtrl=True, redo=True,long_version=False, key_name='dist_sim')
		
	r=feature_cell_extraction.collectingDistances('all_dist_sim.pkl', 
		'../resultData/simulated_traj/simres', qc_filename='../data/qc_simulated.txt',mapping_filename= None, testCtrl=False, redo=True,long_version=False, key_name='dist_sim')

vii. Compute hits
	empirical_qval,siRNAL, exp_hit, siRNA_HC, exp_of_highconfsiRNAs, gene_highconf=feature_cell_extraction.multipleHitDistances('../resultData/simulated_traj/simres','dist_sim', 
		qc_filename='../data/qc_simulated.txt', mapping_filename=None, filename='all_dist_sim', combination='max', redo=False, trad=False, without_mean_persistence=True,
		save=True)
	
viii. STUDY RESULTS
	a. Compare with planned hits: this will give accuracy and precision in comparison with the groundtruth
		simulator.evalWorkflowOutput(exp_hit,siRNA_HC)
		
	b. Trajectory clusters:
		* Collect the data + PCA
			from util import sandbox
			sandbox.generic_single_script('sim_collectTraj', 'python tracking/trajPack/clustering.py --action collectingTrajectories --simulated 1 --outputname all_dist_sim')
		
		**. Cluster number
			from tracking.trajPack import clustering_script
			clustering_script.generationScript('sim_KM',outputname='all_dist_sim', simulated=True)
		
	c. Other approach: 
		simulator.evalTrajectoryClustering(....)
		
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
		
iv. Extract collected distances
	r=feature_cell_extraction.collectingDistances('all_distances_whole_CTRL.pkl', 
		'../resultData/features_on_films', qc_filename='../data/qc_export.txt',
		mapping_filename='../data/mitocheck_siRNAs_target_genes_Ens75.txt', testCtrl=True, redo=True,long_version=False, key_name='distances_whole_5CtrlC')
		
	r=feature_cell_extraction.collectingDistances('all_distances_whole.pkl', 
		'../resultData/simulated_traj/simres', qc_filename='../data/qc_export.txt',
		mapping_filename='../data/mitocheck_siRNAs_target_genes_Ens75.txt', testCtrl=False, redo=True,long_version=False, key_name='distances_whole_5CtrlC')
		
v. Compute hits
	empirical_qval,siRNAL, exp_hit, siRNA_HC, exp_of_highconfsiRNAs, gene_highconf=feature_cell_extraction.multipleHitDistances('../resultData/features_on_films','distances_whole_5CtrlC', 
		qc_filename='../data/qc_export.txt',
		mapping_filename='../data/mitocheck_siRNAs_target_genes_Ens75.txt', 
		filename='all_distances_whole', combination='max', redo=False, trad=True, without_mean_persistence=True,save=True, threshold=0.02)
		
vi. Study results
	from util import sandbox
	sandbox.generic_single_script('real_collectTraj', 'python tracking/trajPack/clustering.py --action collectingTrajectories --outputname all_distances_whole')
	
	PUIS
	
	from tracking.trajPack import clustering_script
	clustering_script.generationScript('real_distances_KM',outputname='all_distances_whole', simulated=False)
	
vii. Cluster trajectories
	from tracking.trajPack import exploiting_clustering
	exploitingKMeans(data,ctrldata, length,ctrlStatus,8, who, genes,sirna, iteration=False, Td=True, show=True)
	f=open('../resultData/features_on_films/labels.pkl', 'w')
	pickle.dump(labels,f); f.close()
	
	
	
	



