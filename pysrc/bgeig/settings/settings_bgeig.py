import getpass, os

outputFolder = "/media/lalil0u/New/projects/Geiger/results"
dataFolder="/media/lalil0u/New/projects/Geiger/nuclei"
allDataFolder="/media/lalil0u/New/projects/Geiger/cytoplasmic_membrane"
    
loadingFolder='../prediction'
    
intensity_qc_file=None
traj_filename = 'traj_intQC_w{}.pkl'
all_trajectory_xml=True

feature_filename='features_intQC_{}.pkl'
outputFile = 'bgeig_features_{}.pkl' 
figname='evolution_features_{}_{}.png'

training = False
predict = True
repeat=True
filtering_fusion=False

time_windows=[(0, 48), (0,24), (12, 36), (24,48),(36, 60)] 

separating_function=lambda x: int(x.split('_')[1])

removeTempFiles=True

redo=False

dict_corresp_nuclei_cyto='/media/lalil0u/New/projects/Geiger/results/dict_nuclei_cytoplasm_object_number.pkl'

objective = [
    #AREAS
             {'name':'cell_area',
              'channel':'primary__primary3',
              'feature':'roisize'
              },
             {'name':'nucleus_area',
              'channel':'secondary__primary3',
              'feature':'roisize'
              },
    #AXIAL RATIOS   
             {'name':'cell_axial_ratio',
              'channel':'primary__primary3',
              'feature':'ellip_axis_ratio'
              },
             {'name':'nuclear_axial_ratio',
              'channel':'secondary__primary3',
              'feature':'ellip_axis_ratio'
              },
    #TEXTURE AUTOCORR. Would be on the whitefield image but we have too much noise. So we need to stay in the fluorescent channel
#              {'name':'cell_autocorr_1',
#               'channel':'tertiary__expanded',
#               'feature':'granu_close_volume_1',
#               'function':(lambda x:x['granu_close_volume_1'][:-1]-x['granu_close_volume_1'][1:])
#               },
#              {'name':'cell_autocorr_2',
#               'channel':'tertiary__expanded',
#               'feature':'granu_open_volume_1',
#               'function':(lambda x:x['granu_open_volume_1'][:-1]-x['granu_open_volume_1'][1:])
#               },
#              {'name':'cell_autocorr_3',
#               'channel':'tertiary__expanded',
#               'feature':'granu_close_volume_3',
#               'function':(lambda x:x['granu_close_volume_3'][:-1]-x['granu_close_volume_3'][1:])
#               },
#              {'name':'cell_autocorr_4',
#               'channel':'tertiary__expanded',
#               'feature':'granu_open_volume_3',
#               'function':(lambda x:x['granu_open_volume_3'][:-1]-x['granu_open_volume_3'][1:])
#               },
    #NUCLEAR POSITION WITH RESPECT TO NUCLEUS
             {'name':'nuclear_pos_min',
              'feature':['bounding_box','center'],
              'function': lambda x: np.minimum(
                                           np.minimum(np.fabs(x['bounding_box']['left']-x['center']['x'])[:,0], 
                                                      np.fabs(x['bounding_box']['right']-x['center']['x'])[:,0]),
                                           np.minimum(np.fabs(x['bounding_box']['top']-x['center']['y'])[:,0],  
                                                      np.fabs(x['bounding_box']['bottom']-x['center']['y'])[:,0])
                                           )
              },
             {'name':'nuclear_pos_max',
              'feature':['bounding_box','center'],
              'function': lambda x: np.maximum(
                                           np.maximum(np.fabs(x['bounding_box']['left']-x['center']['x'])[:,0], 
                                                      np.fabs(x['bounding_box']['right']-x['center']['x'])[:,0]),
                                           np.maximum(np.fabs(x['bounding_box']['top']-x['center']['y'])[:,0],  
                                                      np.fabs(x['bounding_box']['bottom']-x['center']['y'])[:,0])
                                           )
              },
             {'name':'nuclear_pos_ratio',
              'feature':['nuclear_pos_max','nuclear_pos_min'],
              'function': lambda x: x['nuclear_pos_min']/float(x['nuclear_pos_max'])
              },
             ]