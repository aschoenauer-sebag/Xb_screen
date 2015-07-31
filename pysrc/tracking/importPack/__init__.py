FEATURE_NUMBER = 239
k=5#nb de nearest neighbours
dmax=45#distances pour le nearest neighbour
#dist = 30#pour la connexion des tracklets #Alice le 30 juillet, est-ce qu'on utilise cette variable ?
d_bf_move= 50
d_bf_more = 30
XMAX = 1024#attention dans les fichiers hdf5 xmax vaut 1344 et ymax 1024
YMAX = 1344
featuresToKeep = [221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238]


EVENTS = []
EVENTS.append("move")
EVENTS.append("appear")
EVENTS.append("disappear")
EVENTS.append("merge")
EVENTS.append("split")

POSSIBLE_MATCH = []
POSSIBLE_MATCH.append((1,0,1,1,1))
POSSIBLE_MATCH.append((1,0,1,1,1))
POSSIBLE_MATCH.append(None)
POSSIBLE_MATCH.append((1,0,1,1,1))
POSSIBLE_MATCH.append((1,0,1,1,1))

SIZES = []
SIZES.append((1,1))
SIZES.append((0,1))
SIZES.append((1,0))
SIZES.append(("M",1))
SIZES.append((1,"M"))

nb_big_folds=5
nb_folds= 5
c_list = [2.5]#[1,1.5,2]#[0,0.5, 1, 1.5, 2, 2.5, 3]
#bestC=[?, 1.5, ?, 1.5 ,]
#MORPHOLOGY FEATURES
l1=['roisize','perimeter','dist_max','dist_min','dist_ratio', 'circularity', 'irregularity', 'irregularity2', 'n_wdist']
l2=['moment_I1', 'moment_I2', 'moment_I3', 'moment_I4', 'moment_I5', 'moment_I6', 'moment_I7', 'eccentricity', 'gyration_radius', 'gyration_ratio', 'ellip_major_axis', 'ellip_minor_axis',
    'ellip_axis_ratio', 'princ_gyration_x', 'princ_gyration_y', 'princ_gyration_ratio', 'skewness_x', 'skewness_y']
l3=['ch_cc', 'ch_thresh_cc' , 'ch_area_ratio', 'ch_rugosity', 'ch_max_val_0', 'ch_max_val_1', 'ch_max_val_2', 'ch_acd', 'ch_mean_area', 'ch_variance_area']
l4=['granu_close_area_1', 'granu_close_area_2', 'granu_close_area_3', 'granu_close_area_5', 'granu_close_area_7', 'granu_close_volume_1', 'granu_close_volume_2', 'granu_close_volume_3', 'granu_close_volume_5','granu_close_volume_7', 'granu_open_area_1','granu_open_area_2','granu_open_area_3','granu_open_area_5','granu_open_area_7','granu_open_volume_1','granu_open_volume_2','granu_open_volume_3', 'granu_open_volume_5','granu_open_volume_7']
l5=['dyn_distance_nb_max', 'dyn_distance_radius_0', 'dyn_distance_radius_1', 'dyn_distance_radius_2', 'dyn_distance_radius_3']

#index correspondants dans les features CeCog
l_indexes =[236, 232, 11, 12, 13, 10, 165, 166, 230, 215, 216, 217, 218, 219, 220, 221, 19, 43, 44, 21, 22, 20, 234, 235, 233, 237, 238, 2, 8, 1, 7, 3, 4, 5, 0, 6, 9, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 14, 15, 16, 17, 18]