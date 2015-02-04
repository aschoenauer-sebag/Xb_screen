1.Preliminary steps:
	i. Get the data in *.czi files from the microscope
	ii. Do the csv files containing the information for the database regarding wells, clone name, %serum, medium, xenobiotics and doses
		in projects/Xb_screen/protocols_etal/plate_setups
	iii. Do the image extraction in Windows using Zen (Zeiss software)
	
2. Renaming the images to be read by CellCognition
	from util import listFileManagement
	listFileManagement.renameFromZeiss(inputFolder=[folder with the images], plate='[plate name]')

3. Extract features with CellCognition, setting file projects/Xb_screen/settings/settings_all_features_2.conf
	
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
	i. Generate scripts 
	
	[missing info here]
	
	
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
	
	
	