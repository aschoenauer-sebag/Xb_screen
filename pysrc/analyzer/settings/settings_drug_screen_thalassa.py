#SETTINGS FOR aschoenauer@cbio.ensmp.fr sur thalassa

###DIRECTORY SETTINGS
#where the images are /share/data40T/
raw_data_dir = "/share/data40T/Thomas/drug_screen/ORDERED_DRUGS_screen/data"

base_result_dir = '/cbio/donnees/aschoenauer/projects/drug_screen'
base_html = "/cbio/donnees/aschoenauer/public_html"

#where hdf5 files are
raw_result_dir = "/share/data40T/aschoenauer/drug_screen"
#intensity_qc_filename='../data/xb_intensity_qc.pkl'

#Where to save processed results
result_dir = os.path.join(base_result_dir, 'results')

#Where to save data to be put online
media_dir = os.path.join(base_html, 'interface_screen/plates/static')
plot_dir = os.path.join(media_dir, 'plots')
movie_dir = os.path.join(media_dir, "movies")

# Plate setups directory
confDir = os.path.join(base_result_dir, 'plate_setups')

#Renaming plates for filenames and visualization on the interface
def newPlateName(plate,dateFormat = "%m%y%d"):
    import datetime
    date = datetime.datetime.strptime(plate[2:6]+plate[-2:], dateFormat)
    return date.strftime('%d')+date.strftime('%m')+date.strftime('%y')

#information about well folder naming: where should the well number be indicated in the well image folder?
whereWell = lambda x: x[:3]
    
# if Zeiss plate setup did not include some columns, indicate it here
missing_cols = {'11414':(1,2)}

###FEATURES OF INTEREST
#Plate features and channel of extraction
featuresOfInterest = ['Focus']
featureChannels = [0]
#Well features
well_features = ["cell_count", 'Focus_ch1']

###DATA BASE SETTINGS
#Plate name
name = 'Drug screen'
xbscreen = False

#Date format
date_format = '%m%y%d'

###OTHER SETTINGS
###decide if the first well has number zero or 1
startAtZero = False
### is there more than one channel ?
secondaryChannel =False
### do you want to count empty wells according to the plate setup ?
countEmpty = False

density_plot_settings = {
    'min_count': 20,
    'max_count': 600,
    'min_class': 0,
    'max_class': 0.6,
    'min_circularity': 0.1,
    'max_circularity': 0.8,
    'min_proliferation': 1, 
    'max_proliferation': 3.0,
    'min_death': 0.1, 
    'max_death': 2.5
}

well_plot_settings={
                    'cell_count':(50, 700),
                    'circularity':(0, 0.75),
                    'Focus_ch1':(0,1),
                    'Nuclear morphologies':(0,1)
}

TRANSLATION_WHOLENAMED = {
    'max_Interphase': 'Interphase',
    'max_Large': 'Large',
    'max_Shape1': 'Binuclear',
    'max_Shape3': 'Polylobed',
    'max_Grape': 'Grape',
    'max_Elongated': 'Elongated',
    'max_Metaphase': 'Metaphase',
    'max_Anaphase': 'Anaphase',
    'max_MetaphaseAlignment': 'MAP',
    'max_Prometaphase': 'Prometaphase',
    'max_ADCCM': 'ADCCM',
    'max_Apoptosis': 'Cell Death',
    'max_Hole': 'Hole',
    'max_Folded': 'Folded',
    'max_SmallIrregular': 'Small Irregular',
    'max_Artefact': 'Artefact',
    'max_UndefinedCondensed': 'Condensed',
    'maxSum_Dynamic': 'Dynamic changes',
    'maxSum_MitosisPhenotype': 'Mitotic Delay',
    'proliferationDiff': 'Proliferation Difference',
    'proliferation': 'Proliferation',
    'subTrack_dist_mean_norm': 'Migration: Distance', 
    'frameToFrame_max': 'Migration: Speed'
    }
