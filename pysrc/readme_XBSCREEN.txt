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
		