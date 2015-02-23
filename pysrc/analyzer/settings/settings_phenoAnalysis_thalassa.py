#SETTINGS FOR lalil0u@Trulove

###DIRECTORY SETTINGS
#where the images are
base_result_dir ='/cbio/donnees/aschoenauer/projects/Xb_screen'

#Where plate setups are
plate_setups_folder=os.path.join(base_result_dir, 'protocols_etal/plate_setups')

#Where to load processed results
loadingFolder = os.path.join(base_result_dir, 'dry_lab_results')

ok_wells_asLIST = '../data/xb_OK_wells_asLIST.pkl'
ok_wells_asDICT= '../data/xb_OK_wells_asDICT.pkl'
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

