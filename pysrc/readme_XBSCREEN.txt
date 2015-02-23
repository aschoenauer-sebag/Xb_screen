1.Preliminary steps:
	i. Get the data in *.czi files from the microscope
	ii. Do the csv files containing the information for the database regarding wells, clone name, %serum, medium, xenobiotics and doses
		in projects/Xb_screen/protocols_etal/plate_setups
	iii. Do the image extraction in Windows using Zen (Zeiss software)
	
2. Renaming the images to be read by CellCognition
	from util import listFileManagement
	listFileManagement.renameFromZeiss(inputFolder=[folder with the images], plate='[plate name]')

3. Extract features with CellCognition, settings file CellCognition projects/Xb_screen/settings/settings_all_features_2.conf

Rq (06/02/2015) : now possible to do it on the cluster. 
	Settings file for using CellCognition on the cluster pysrc/analyzer/settings/pbs_Xb_screen1.py
	Command python analyzer/Cecog_script_generation.py -b analyzer/settings/pbs_Xb_screen1.py -p [filename with relevant wells]
	eg for all wells of the screen that passed the QC : data/all_QC_wells_XBSC.pkl
	
4. Enter the info in the db, compute the images with regard to cell count, number of out of focus objects, etc
	from analyzer import interface
	html=interface.HTMLGenerator(settings_file='analyzer/settings/settings_xb1.py')
	html(plateL=['[plate name]'], doQC=True, fillDBOnly=False)
	
5. Copy the data to Olympia and save the images and CZI files to CBIO

6. Compute the cell trajectories
	wells=range(1,[well_number])
	for w in wells:
		ww='{:>05}_01.hdf5'.format(w)
		%run tracking/trajPack/parallel_trajFeatures -p [plate name] -w $ww -d /media/lalil0u/New/projects/Xb_screen/plates__all_features_2bis -c 0
		%run tracking/trajPack/parallel_trajFeatures -p [plate name] -w $ww -d /media/lalil0u/New/projects/Xb_screen/plates__all_features_2bis -c 1 --ff 0 -n features_intQC_{{}}.pkl
		
7. Launch hit detection step
	i. Generate scripts for feature_cell_extraction
	
	
	ii. Generate scripts for phenotype analysis
	from analyzer import phenotype_analysis
	phenotype_analysis.script()
	
	
8. Looking at PCAed features using ctrls
	from analyzer import quality_control, xb_analysis
	exp_list=quality_control.usable_XBSC('DMSO', 0)
	exp_list.extend(quality_control.usable_XBSC('Nonane', 0))
	
	r_D, who,ctrlStatus, length_D, xb, others, time_length_D=xb_analysis.xbConcatenation(folder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results/', exp_list=exp_list)
	r=np.hstack((r_D[:,:len(featuresNumeriques)], r_D[:,featuresSaved.index('mean straight')][:,np.newaxis]))
	nr=(r-mean)/std
	pca=PCA(n_components=15, whiten=False)
	pnr=pca.fit_transform(nr)
	all_r, who,ctrlStatus, length, xb, others, time_length=xb_analysis.xbConcatenation(folder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results/')
	all_r=np.hstack((all_r[:,:len(featuresNumeriques)], all_r[:,featuresSaved.index('mean straight')][:,np.newaxis]))
	
	pnall_r=pca.transform((all_r-mean)/std)
	
#puis pour regarder les plots : on peut considérer que l'on a des infos sur quatre composantes pour chaque puits, et regarder ce que cela donne au niveau des conditions.
 
	r={k:np.array([pnall_r[np.sum(length[:j]):np.sum(length[:j+1]),k] for j in range(len(length))]) for k in range(4)}
	phenotype_analysis.plotResults(None, r, who, dose_list, xb, None, None, features=True)
	
#Ne donne rien de très concluant. Donc autre idée, regarder comment les trajectoires se répartissent dans les clusters trouvés avec Mitocheck. DONC :

#En fait je l'avais déjà fait : c'est à peu près la fonction xb_analysis.comparTrajectoriesMitoXB().......
	
	f=open('../resultData/features_on_films/all_distances_whole_dataonly.pkl')
	arr=pickle.load(f);f.close()
	data=arr[0]; data=np.hstack((data[:,:len(featuresNumeriques)], data[:,featuresSaved.index('mean straight')][:,np.newaxis]))
	
	pca=PCA(n_components=15, whiten=False)
	mean=np.mean(data,0)
	std=np.std(data,0)
	ndata=(data-mean)/std
	pca.fit(ndata)
	
	all_r, who,ctrlStatus, length, xb, others, time_length=xb_analysis.xbConcatenation(folder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results/')
	all_r=np.hstack((all_r[:,:len(featuresNumeriques)], all_r[:,featuresSaved.index('mean straight')][:,np.newaxis]))
	#to put everybody in the same referencial
	pnall_r=pca.transform((all_r-np.mean(all_r,0))/np.std(all_r,0))
	model=MiniBatchKMeans(n_clusters=8, batch_size = 2000, init='k-means++',n_init=1000,max_iter=1000, max_no_improvement=100, compute_labels = True)
	model.fit((pndata/np.std(pndata,0))[:,:7])
	std_p=np.std(pndata,0)
	xb_screen_labels=model.predict((pnall_r/std_p)[:,:7])

	#Sauvegardes:			
	f=open('../../../projects/Xb_screen/dry_lab_results/track_predictions__settings2/used_MotIW_pca.pkl', 'w')
	pickle.dump((pca, mean, std, std_p),f)
	f.close()
	
	f=open('../../../projects/Xb_screen/dry_lab_results/track_predictions__settings2/xb_screen_labels.pkl', 'w')
	pickle.dump(xb_screen_labels, f)
	f.close()
	percentages=exploiting_clustering.filmCharacterization(xb_screen_labels, length, ctrlStatus, others, colors=None, show=False,plotcov=False)
	
#Maintenant on a de nouveau pour chaque puits huit caractéristiques, qui sont les pourcentages de trajectoires dans chacun des clusters. On peut donc regarder cela
#de la même manière avec des boxplots
	phenotype_analysis.plotResults(None, r, who, dose_list, xb, None, None, features=True)
	
#To compute Fisher's pvalues on trajectory distributions in Mitocheck trajectory clusters
	r1,r2=xb_analysis.comparClusterDistributions(labels, compound_list, who, ctrlStatus, length, dose_list, n_cluster=8)
	
#PLATE NORMALIZATION (instead of negative control normalization)

	file='projects/Xb_screen/dry_lab_results/track_predictions__settings2/xb_screen_Fisherspval_plate_normalization.pkl'
	result, xb, who, length,conditions = pickle.load(open(file, 'r'))
	
	
	
	
	
	
	
	
	

	