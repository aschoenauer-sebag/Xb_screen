
# for the pilot screen

# directory settings
baseDir = '/data/twalter/data/results_drug_pilot'
result_dir = os.path.join(baseDir, 'flat_files_structure')
local_plot_dir = '/g/mitocheck/Thomas/drugs_html'
html_dir = local_plot_dir

# full filename of result pickle file
id_result_filename = os.path.join(baseDir, 'results', 'id_results.pickle')

# labtek lists
labtekBaseL = ['1000']
labtekL = ['LT1000_01', 'LT1000_02']

# control positions
negCtrlPosL = ['012', '013', '014', '015', '016', '017', '018', '019', '020', '021', '022', '034', '035', '036', '037', '038', '039', '040', '041', '042', '043', '044', '078', '079', '080', '081', '082', '083', '084', '085', '086', '087', '088']
emptyL = ['100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '166', '167', '168', '169', '170', '171', '172', '173', '174', '175', '176', '188', '189', '190', '191', '192', '193', '194', '195', '196', '197', '198', '210', '211', '212', '213', '214', '215', '216', '217', '218', '219', '220', '232', '233', '234', '235', '236', '237', '238', '239', '240', '241', '242', '254', '255', '256', '257', '258', '259', '260', '261', '262', '263', '264', '276', '277', '278', '279', '280', '281', '282', '283', '284', '285', '286', '298', '299', '300', '301', '302', '303', '304', '305', '306', '307', '308']

controlD = {'negative': negCtrlPosL, 'empty': emptyL}

labeledPosD = {'d': negCtrlPosL,
               '-': emptyL}

# labtekConfiguration
confDir = os.path.join(baseDir, 'configuration')
confFilename = 'labtekConfiguration.pickle'

# prediction suffix
prediction_suffix='_prediction_track.dat'

# labtek layout
nb_row = 14
nb_col = 22

gene_key = 'drug'
sirna_key= 'concentration'

key_level1 = gene_key
key_level2 = sirna_key

density_plot_settings = {
    'min_count': 15,
    'max_count': 800,
    'min_proliferation': 0.0, 
    'max_proliferation': 4.0
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