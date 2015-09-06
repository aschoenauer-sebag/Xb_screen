#VARS FOR DRUG SCREEN PROJECT

import cPickle as pickle
from collections import Counter

couleurs=['green', 'red', 'yellow', 'blue', 'purple', 'orange', 'black', 'cyan']
couleurs.append("#9e0142")
couleurs.append("#5e4fa2")
DRUGS=['Acyclovir',
     'Adenosine3',
     'Aminopurine6benzyl',
     'Anisomycin',
     'Azacytidine5',
     'Camptothecine(S,+)',
     'Daunorubicinhydrochloride',
     'Dexamethasoneacetate',
     'Doxorubicinhydrochloride',
     'Epiandrosterone',
     'Etoposide',
     'Hesperidin',
     'Idoxuridine',
     'JNJ7706621',
     'MLN8054',
     'Methotrexate',
     'Nocodazole',
     'Paclitaxel',
     'R763',
     'Ribavirin',
     'Sulfaguanidine',
     'Sulfathiazole',
     'Thalidomide',
     'VX680',
     'Zidovudine,AZT']

CLASSES=['Interphase',
     'Large',
     'Elongated',
     'Binucleated',
     'Polylobed',
     'Grape',
     'Metaphase',
     'Anaphase',
     'MetaphaseAlignment',
     'Prometaphase',
     'ADCCM',
     'Apoptosis',
     'Hole',
     'Folded',
     'SmallIrregular']

KNOWN_TARGETS = {'Anisomycin': ['RPL10L', 'RPL13A', 'RPL23',  'RPL15',
                          'RPL19','RPL23A','RSL24D1','RPL26L1','RPL8','RPL37','RPL3','RPL11','NHP2L1'],
           'Azacytidine5':['DNMT1'],
           'Camptothecine(S,+)':['TOP1'],
            'Daunorubicinhydrochloride':['TOP2A', 'TOP2B'],
            'Doxorubicinhydrochloride':['TOP2A', 'TOP2B'],
           'Etoposide':['TOP2A', 'TOP2B'],
           'JNJ7706621' : ['AURKA', 'AURKB', 'CDK1', 'CDK2', 'CDK3', 'CDK4', 'CDK6'],
           'MLN8054':['AURKA', 'AURKB'],
           'Nocodazole':['HPGDS', 
                         'TUBB2A',
                        'TUBB4A',
                        'TUBB4B',
                        'TUBB6',],
           'Paclitaxel':['TUBB1', 'BCL2','NR1I2','MAPT','MAP4','MAP2'],
           'VX680':['AURKA', 'AURKB'],
           }
DISTANCES = {'N_pheno_score':'Normalized\n phenotypic score',
 'nature':'Phenotypic\n trajectory',
 'ttransport_MAX':'Max time \n Sinkhorn div.',
 'U_pheno_score':'Phenotypic score',
 'ttransport_INT':'Sum of time\n Sinkhorn div.',
 'transport':'Global\n Sinkhorn div.'}


f=open('/media/lalil0u/New/projects/drug_screen/results/DS_pheno_scores.pkl')
_,_,_, exposure_=pickle.load(f); f.close()
PASSED_QC_COND=Counter(exposure_)