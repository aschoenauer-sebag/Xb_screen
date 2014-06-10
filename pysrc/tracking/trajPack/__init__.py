

#MOMENTS : nb de puissances des distances parcourues que l'on va regarder pour calculer movement type
moments = range(1,7)

#DENSITIES : a cb de pixels va-t-on regarder pour la densite maximale autour d'une cellule, sur toute sa trajectoire ? 
densities =[50.0]#, 60.0, 80.0, 100.0]

#RAYON : rayon des boules pour l'entropie d'une trajectoire
rayon = [0,14]#1,2,3,4, 5, 7, 9, 11, 14]

#WINDOW : longueur d'une fenetre temporelle pour regression et RSS (autoroute de Thomas)
windows = [10]#[6, 8, 10]

#TIMEDEPT : nb de frames ou l'on va regarder persistance, signed turning angle et turning angle
TIMEDEPTH=1
#FEATURES : max
featuresListArray= [#'FLATmoments', 
                    'convex hull area','norm convex hull area', 
                    'corrected straightness index', 
                    'diffusion coefficient', 'direction', 
                    'effective space length', 'effective speed', 
                    'entropy1', 'entropy2', 'entropy3', 'entropy4', 'entropy5', 'entropy6', 'entropy7', 'entropy8', 
                    'entropy9', 'entropy10', 'largest move', 'mean instantaneous acceleration', 
                    'mean instantaneous speed', 'stdev instantaneous speed',
                    'mean squared displacement', 'movement type', 
                    'persistence1', 'stdev persistence1',
                    'persistence2', 'persistence3', 'persistence4', 'persistence5', 
                    'RSS window 6', 'RSS window 8', 'RSS window 10',
                    'signed turning angle1', 'signed turning angle2', 'signed turning angle3', 'signed turning angle4', 'signed turning angle5', 
                    'sinuosity',  'time length', 'total space length', 
                    'turning angle1', 'turning angle2', 'turning angle3', 'turning angle4', 'turning angle5', 'app', 'disp',
                    'mean_d50.0', 'max_d50.0','std_d50.0', 'min_d50.0',  
                    'mean_d60.0', 'std_d60.0', 'min_d60.0', 'max_d60.0', 
                    'mean_d80.0', 'std_d80.0', 'min_d80.0', 'max_d80.0', 
                    'mean_d100.0', 'std_d100.0', 'min_d100.0', 'max_d100.0'
                    ]

#FEATURES utilisees pour MoGDIW
sdFeatures1=[0, 1, 2, 3, 5, 6, 7, 16, 17, 18, 19, 20, 21, 23, 24, 31, 32, 37, 39, 40] #34, 36, 37]#indices sans avoir rien enleve
features1=['convex hull area', 'norm convex hull area',
       'corrected straightness index', 'diffusion coefficient',
       'effective space length', 'effective speed', 
       'entropy1', 'entropy2',#rayons 0, 14 ci-dessous 
       'largest move', 
       'mean instantaneous acceleration', 'mean instantaneous speed', 'stdev instantaneous speed',
       'mean squared displacement', 
       'persistence1', 'stdev persistence1',
       'RSS window 10', 'signed turning angle1', 'sinuosity',
       'total space length', 'turning angle1']
logTrsf = [1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0]

#FEATURES pour melanger histogrammes et features numeriques
featuresSaved=[
               #i. features that i'll use in the analysis
               'ball number1','ball number2',
               'norm convex hull area',
               'corrected straightness index', 
               'diffusion coefficient',
               'effective space length', 'effective speed', 
               'entropy1', 'entropy2',#rayons 0, 14 ci-dessus 
               'largest move', 
               'mean squared displacement',
               'movement type', 
               'signed turning angle', 
               'diffusion adequation',
       
       #ii.features that will not be used in the analysis
       'mvt type adequation',
       'mean acceleration',  'mean displacements',  'mean persistence',  'mean straight',  'mean turningangle',#we want this for analyzing the feature roughly
       'convex hull area',
       'correlation speed and persistence',
       'slope SP',
       'total space length']

featuresNumeriques=['ball number1','ball number2', 
                    'norm convex hull area',
                    'corrected straightness index', 
                    'diffusion coefficient',
                    'effective space length', 'effective speed', 
                    'entropy1', 'entropy2',#rayons 0, 14 ci-dessus 
                    'largest move', 
                    'mean squared displacement',
                    'movement type', 
                    'signed turning angle',
                    'diffusion adequation'
                    ]

featuresHisto=['acceleration', 'displacements', 'persistence', 'straight', 'turningangle']

#transformation ou pas pour les featuresNumeriques. Verification refaite apres modification du calcul des features le 6/2/14 pour ne pas prendre
#les frames avant/apres merge/split
histLogTrsf = [1,1,0,
               1, 1, 1,
               0, 0, 1,
               1, 1, 1, 1, 0]
histLogTrsf_meanHistFeat=[0,1,0,1,0]

#RSS window 6, 8, 10 : oui oui oui

d_ctrl={'LT0007_26': ['00074_01', '00015_01'],
         'LT0013_42': ['00026_01', '00304_01'], 
         'LT0012_29': ['00352_01', '00304_01'], 
         'LT0024_41': ['00074_01', '00063_01'], 
         'LT0039_45': [], 'LT0008_26': ['00015_01'], 
         'LT0038_27': [], 'LT0041_35': [], 
         'LT0010_27': ['00026_01', '00304_01'], 
         'LT0037_52': ['00026_01', '00315_01']}

