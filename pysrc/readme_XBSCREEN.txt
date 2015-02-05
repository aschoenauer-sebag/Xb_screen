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
	
	
	making movie number 2 LT0100_03--050
converting images ... 
generating movie ... 
mencoder "mf://../resultData/Munich_plots/mean_straight/temp/*.jpg" -mf fps=3 -o ../resultData/Munich_plots/mean_straight/LT0100_03/LT0100_03--050.mov -ovc xvid -oac copy -xvidencopts fixed_quant=2.5
movie generated: ../resultData/Munich_plots/mean_straight/LT0100_03/LT0100_03--050.mov
MEncoder 1.1-4.1.2 (C) 2000-2012 MPlayer Team

WARNING: OUTPUT FILE FORMAT IS _AVI_. See -of help.
success: format: 16  data: 0x0 - 0x0
MF file format detected.
[mf] search expr: ../resultData/Munich_plots/mean_straight/temp/*.jpg
============ Sorry, this file format is not recognized/supported =============
=== If this file is an AVI, ASF or MPEG stream, please contact the author! ===
Cannot open demuxer.

Exiting...
mencoder "mf://../resultData/Munich_plots/feature_movies_2/mean_straight/temp/*.png" -mf fps=3 -o ../resultData/Munich_plots/feature_movies_2/mean_straight/meanstraight_LT0100_03--050.mov -ovc xvid -oac copy -xvidencopts fixed_quant=2.5
movie generated: ../resultData/Munich_plots/feature_movies_2/mean_straight/meanstraight_LT0100_03--050.mov
MEncoder 1.1-4.1.2 (C) 2000-2012 MPlayer Team

WARNING: OUTPUT FILE FORMAT IS _AVI_. See -of help.
success: format: 16  data: 0x0 - 0x0
MF file format detected.
[mf] search expr: ../resultData/Munich_plots/feature_movies_2/mean_straight/temp/*.png
[mf] number of files: 93 (744)
[demux_mf] file type was not set! trying 'type=png'...
VIDEO:  [MPNG]  0x0  24bpp  3.000 fps    0.0 kbps ( 0.0 kbyte/s)
[V] filefmt:16  fourcc:0x474E504D  size:0x0  fps:3.000  ftime:=0.3333
xvid: using library version 1.3.2 (build xvid-1.3.2)
Opening video filter: [expand osd=1]
Expand: -1 x -1, -1 ; -1, osd: 1, aspect: 0.000000, round: 1
==========================================================================
Opening video decoder: [ffmpeg] FFmpeg's libavcodec codec family
libavcodec version 54.23.100 (internal)
Selected video codec: [ffpng] vfm: ffmpeg (FFmpeg PNG)
==========================================================================
Could not find matching colorspace - retrying with -vf scale...
Opening video filter: [scale]
Movie-Aspect is undefined - no prescaling applied.
[swscaler @ 0xdfa500]BICUBIC scaler, from rgb24 to yuv420p using MMX2
videocodec: XviD (1344x1024 fourcc=44495658 [XVID])
xvid: par=0/0 (vga11), displayed=1344x1024, sampled=1344x1024
xvid: Fixed Quant Rate Control -- quantizer=5/2=2.50
New_Face failed. Maybe the font path is wrong.
Please supply the text font file (~/.mplayer/subfont.ttf).
subtitle font: load_sub_face failed.
New_Face failed. Maybe the font path is wrong.
Please supply the text font file (~/.mplayer/subfont.ttf).
subtitle font: load_sub_face failed.
Writing header...
ODML: vprp aspect is 16384:12483.
Writing header...
	