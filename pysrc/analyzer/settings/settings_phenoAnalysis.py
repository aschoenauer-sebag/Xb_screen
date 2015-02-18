#SETTINGS FOR lalil0u@Trulove

###DIRECTORY SETTINGS
#where the images are
raw_data_dir = "/media/lalil0u/TOSHIBA EXT/images_A"#"/media/lalil0u/XB SCREEN/ok cbio/images"#"/media/lalil0u/New/data/Xb_screen/Images_Cecog"#
base_result_dir = '/media/lalil0u/New/projects/Xb_screen'

#Where to load processed results
loadingFolder = os.path.join(base_result_dir, 'dry_lab_results')
#Where to save results
savingFolder=os.path.join(loadingFolder, 'phenotype_analysis')
outputFile='phenoAnalysis_plateNorm_{}_{}.pkl'

#Local regression parameters
h=10
deg=3
#Default phenotypes of interest
pheno_list=['Anaphase_ch1', 'Apoptosis_ch1', 'Folded_ch1', 'Interphase_ch1', 'Metaphase_ch1',
            'Polylobbed_ch1', 'Prometaphase_ch1', 'WMicronuclei_ch1']
#Regarding normalization
norm='plate'
n_n=37

