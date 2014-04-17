
# for the pilot screen

# directory settings
baseDir = '/data/twalter/data/results_drug_pilot'
result_dir = os.path.join(baseDir, 'flat_files_structure')
local_plot_dir = '/g/mitocheck/Thomas/drugs_html'
html_dir = local_plot_dir

# full filename of result pickle file
id_result_filename = os.path.join(baseDir, 'results', 'id_results.pickle')

# labtek lists
labtekBaseL = ['0699']
labtekL = ['LT0699_01', 'LT0699_02', 'LT0699_03', 'LT0699_04']

# control positions
negCtrlPosL = ['001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '028', '029', '030', '031', '032', '033', '034', '035', '036', '048', '049', '072', '073', '096', '097', '120', '121', '144', '145', '168', '169', '192', '193', '216', '217', '240', '241', '264', '265', '288', '289', '312', '313', '336', '337', '360', '361', '362', '363', '364', '365', '366', '367', '368', '369', '370', '371', '372', '373', '374', '375', '376', '377', '378', '379', '380', '381', '382', '383', '384']

controlD = {'negative': negCtrlPosL}

labeledPosD = {'d': negCtrlPosL}

# labtekConfiguration
confDir = os.path.join(baseDir, 'configuration')
confFilename = 'labtekConfiguration.pickle'

# prediction suffix
prediction_suffix='_prediction_clean.dat'

# labtek layout
nb_row = 32
nb_col = 12

gene_key = 'drug'
sirna_key= 'concentration'

key_level1 = gene_key
key_level2 = sirna_key

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