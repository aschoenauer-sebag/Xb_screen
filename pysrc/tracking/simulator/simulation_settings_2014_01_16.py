import os
#out_folder = '/media/lalil0u/New/workspace2/Tracking/resultData/simulated_traj'
#out_folder = '/cbio/donnees/aschoenauer/workspace2/Tracking/resultData/simulated_traj'
#plot_folder = os.path.join(out_folder, 'plots')
out_folder = '/Users/twalter/data/migration_sim/data/output'
plot_folder = '/Users/twalter/data/migration_sim/data/plots'
real_data_folder = '/Users/twalter/data/migration_sim/data/real'

adapt_axis = True
grid = True

#L = {
#    'mean': 21,
#    'stdev': 5
#    }
#
length_file = os.path.join(real_data_folder, 'lengthDistribution.pkl')
nb_traj_file = os.path.join(real_data_folder, 'nb_tracks.pickle')

L = {
    'distribution': 'Uniform',
    'min': 12,
    'max': 40,
    }

simulator_settings = {
'movement_normal': {'N': 1,
                    'radius': {'distribution': 'Normal', 'mean': 15, 'stdev': 3},
                    'angle': {'distribution': 'Uniform', 'min': 0, 'max': 360},
                    },
'movement_fast': {'N': 1,
                  'radius': {'distribution': 'Normal', 'mean': 25, 'stdev': 8},
                  'angle': {'distribution': 'Uniform', 'min': 0, 'max': 360},
                  },
'directional': {'N': 1,
                'angle_0': {'distribution': 'Uniform', 'min': 0, 'max': 360}, # mean value
                'angle': {'distribution': 'Normal', 'stdev': 25.0},
                'radius': {'distribution': 'Normal', 'mean': 15, 'stdev': 3},
                },
'directional_bended': {'N': 20,
                       'angle_0': {'distribution': 'Uniform', 'min': 0, 'max': 360}, # start value
                       'angle': {'distribution': 'Normal_Moving_Mean', 'stdev': 25.0},
                       'angle_delta': {'distribution': 'Normal', 'mean': 50.0, 'stdev': 1.0},
                       'radius': {'distribution': 'Normal', 'mean': 15, 'stdev': 3},
                },

# trajectory consisting in two different speeds. 
#'radius_switch': {'N': 20,
#                  'radius': {'distribution': 'Normal_two_states',
#                             'radius1': {'distribution': 'Normal', 'mean': 15, 'stdev': 3},
#                             'radius2': {'distribution': 'Normal', 'mean': 35, 'stdev': 3},
#                             'trans_count_0_1': {'min': 5, 'max': 10},
#                             'trans_count_1_0': {'min': 1, 'max': 5},},
#                  'angle': {'distribution': 'Uniform', 'min': 0, 'max': 360},                             
#                  },

# indicates directed movement in +- 180 degrees.
'angle_switch': {'N': 1,
                 'radius': {'distribution': 'Normal', 'mean': 15, 'stdev': 3},
                 'angle_0': {'distribution': 'Uniform', 'min': 0, 'max': 360}, # mean value
                 'angle': {'distribution': 'Normal', 'stdev': 25.0, 
                           'switch': {'trans_count_0_1': {'min': 10, 'max': 15},
                                      'trans_count_1_0': {'min': 5, 'max': 12}}}
                  },                                    

# mode_switch indicates the presence of 2 movements, one directed, the other random, with different speeds respectively.
'mode_switch': {'N': 1,
                'common_switch': {'trans_count_0_1': {'min': 8, 'max': 15},
                                  'trans_count_1_0': {'min': 3, 'max': 5}},     
                'radius1': {'distribution': 'Normal', 'mean': 15, 'stdev': 3},
                'radius2': {'distribution': 'Normal', 'mean': 35, 'stdev': 3},                
                'angle1': {'distribution': 'Uniform', 'min': 0, 'max': 360},
                'angle2': {'distribution': 'Normal', 'stdev': 25.0},
                },
}


probabilities = {   
                 'normal': {'movement_normal': 0.92,
                            'movement_fast': 0.05,
                            'directional_bended': 0.01,
                            'angle_switch': 0.01,
                            'mode_switch': 0.01,
                            },
                 'directed1': {'movement_normal': 0.64,
                               'movement_fast': 0.05,
                               'directional_bended': 0.15,
                               'angle_switch': 0.15,
                               'mode_switch': 0.01,
                               },
                 'directed2': {'movement_normal': 0.7,
                               'movement_fast': 0.05,
                               'directional_bended': 0.2,
                               'angle_switch': 0.04,
                               'mode_switch': 0.01,
                               },
                 'fast': {'movement_normal': 0.6,
                          'movement_fast': 0.3,
                          'directional_bended': 0.01,
                          'angle_switch': 0.01,
                          'mode_switch': 0.08,
                          },
                 'switch_fast': {'movement_normal': 0.7,
                                 'movement_fast': 0.1,
                                 'directional_bended': 0.01,
                                 'angle_switch': 0.01,
                                 'mode_switch': 0.18,
                                 },  
                 }



