#SETTINGS FOR lalil0u@Trulove

###DIRECTORY SETTINGS
raw_data_dir = "/media/lalil0u/New/data/Xb_screen/Screen_trials"

base_result_dir = '/media/lalil0u/New/projects/Xb_screen'
raw_result_dir = os.path.join(base_result_dir, 'plates')
result_dir = os.path.join(base_result_dir, 'results')
plot_dir = os.path.join(result_dir, 'plots')
movie_dir = os.path.join(result_dir, "movies")

# Plate setups directory
confDir = os.path.join(base_result_dir, 'plate_setups')

###DEFAULT PLATE
plate = '11414'
# if Zeiss plate setup did not include some columns, indicate it here
missing_cols = {'11414':(1,2)}

###FEATURES OF INTEREST
#Plate features and channel of extraction
featuresOfInterest = ["circularity"]
featureChannels = [1]
#Well features
well_features = ["cell_count", 'circularity']

###DATA BASE SETTINGS
#Plate name
name = 'Trial'

#Date format
date_format = '%d%m%y'

###OTHER SETTINGS
###decide if the first well has number zero or 1
startAtZero = False
### is there more than one channel ?
secondaryChannel =True
### do you want to count empty wells according to the plate setup ?
countEmpty = False
### do you want to redo videos that have already been extracted ?
redoMovies = False

density_plot_settings = {
    'min_count': 20,
    'max_count': 600,
    'min_circularity': 0.1,
    'max_circularity': 0.7,
    'min_proliferation': 0.1, 
    'max_proliferation': 2.5,
    'min_death': 0.1, 
    'max_death': 3.0
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
    
# color dictionnary for the time curve plots
COLORD = {'Interphase': 'chocolate3',
          'Large': 'darkgoldenrod1',
          'Elongated': 'chocolate4',
          'StrangeInterphase': 'cornsilk3',
          
          'Shape': 'blue',
          'Shape1': 'blue',
          'Shape3': 'deepskyblue',
          'Grape': 'dodgerblue3',

          #'MitosisPhenotype': 'green',
          'MitosisPhenotype': 'chartreuse',          
          'Metaphase': 'forestgreen',
          'Anaphase': 'darkolivegreen1',
          'MetaphaseAlignment': 'chartreuse',
          'Prometaphase': 'palegreen2',
          'ADCCM': 'lightgreen',
          
          'Apoptosis': 'red',
          
          'Dynamic': 'magenta',
          'Hole': 'orchid1',
          'Folded': 'purple',
          'SmallIrregular': 'magenta3',

          'Artefact': 'grey',
          'UndefinedCondensed': 'red4',
          
          'proliferationDiff': 'cyan2',
          'proliferation': 'cyan3',
          
          'subTrack_dist_mean_norm': 'yellow', 
          'frameToFrame_max': 'orange'

          }