#SETTINGS FOR aschoenauer@cbio.ensmp.fr sur thalassa

# directory settings
#raw_data_dir = "/media/lalil0u/New/data/Xb_screen/Screen_trials"
base_result_dir = '/cbio/donnees/aschoenauer/projects/Xb_screen'

raw_result_dir = "/share/data20T/mitocheck/Alice/Xb_screen/results"
result_dir = os.path.join(base_result_dir, 'html_results')
# Plate Configuration
confDir = os.path.join(base_result_dir, 'plate_setups')

featuresOfInterest = ["circularity"]
featureChannels = [1]

# moviedir = /g/mitocheck/Thomas/results_drug_screen/LT0900_01/movies
#result_dir = '/g/mitocheck/Thomas/results_drug_screen'
local_plot_dir = os.path.join(result_dir, 'plots')
html_dir = result_dir
#movie_raw_dir = os.path.join(

# full filename of result pickle file
id_result_filename = os.path.join(result_dir, 'id_results.pickle')

# plate lists
plate = ['11414']
## if Zeiss plate setup did not include some columns, indicate it here
missing_cols = {'11414':(1,2)}

#Different settings
###decide if the first well has number zero or 1
startAtZero = False

### is there explotation of data to do on more than one channel ?
secondaryChannel =True

### do you want to count empty wells according to the plate setup ?
countEmpty = False

### not effect
saveResult = False
join_subwells = True
multiple_positions = True
#position_regex = re.compile('--W(?P<Well>\d+)_P(?P<Position>\d+)--')

#INFOS THAT SHOULD BE IN THE CSV FILE
# control positions : 
negCtrlPosL = []# ['%03i' % i for i in range(301, 312)]
#negCtrlPosL = ['001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '028', '029', '030', '031', '032', '033', '034', '035', '036', '048', '049', '072', '073', '096', '097', '120', '121', '144', '145', '168', '169', '192', '193', '216', '217', '240', '241', '264', '265', '288', '289', '312', '313', '336', '337', '360', '361', '362', '363', '364', '365', '366', '367', '368', '369', '370', '371', '372', '373', '374', '375', '376', '377', '378', '379', '380', '381', '382', '383', '384']
# labtek layout
#nb_row = 12
#nb_col = 12
#
#gene_key = 'drug'
#sirna_key= 'concentration'
#
#key_level1 = gene_key
#key_level2 = sirna_key

controlD = {'negative': negCtrlPosL}

labeledPosD = {'d': negCtrlPosL}

# prediction suffix
prediction_suffix='_prediction_track.dat'


density_plot_settings = {
    'min_count': 20,
    'max_count': 600,
    'min_proliferation': 0.1, 
    'max_proliferation': 5.0
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