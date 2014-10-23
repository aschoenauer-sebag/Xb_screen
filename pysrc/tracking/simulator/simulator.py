import os, re, sys, time
import cPickle as pickle
import numpy as np
from collections import Counter
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
#from PyPack.fHacktrack2 import ensTraj, trajectoire
#from trajPack.trajFeatures import histogramPreparationFromTracklets
from util.settings import Settings

import pdb
import plotter_stats

from util.sandbox import accuracy_precision

def evalWorkflowOutput(exp_hit,siRNA_hit,folder='../resultData/simulated_traj/simres',  num_replicates=[1,2,3]):
#Evaluating the workflow at the level of the experiment
    
    annotations = sorted(os.listdir(os.path.join(folder, 'annotations')))
    truth=[]; types={}; normals=0
    for annotation in annotations:
        f=open(os.path.join(folder, "annotations", annotation))
        ann=pickle.load(f)
        f.close()
        
        plate = annotation.split('_')[-1][:6]
        for replicate in num_replicates:
            try:
                truth.extend(['{}_{:>02}--{:>03}'.format(plate,replicate, w) for w in ann[replicate-1] if ann[replicate-1][w] not in ["control", "normal"]])
                types.update({'{}_{:>02}--{:>03}'.format(plate,replicate, w):ann[replicate-1][w] for w in ann[replicate-1]})
            except KeyError:
                print plate
                continue
            else:
                normals+=len([w for w in ann[0] if ann[0][w]=='normal'])
    accuracy_precision(normals, truth, exp_hit)

    #retourne les faux negatifs et les faux positifs
    print Counter((types[el] for el in truth if el not in exp_hit)), Counter((types[el] for el in exp_hit if el not in truth))

#Evaluating the workflow at the level of the siRNA

    f=open(os.path.join(folder,"siRNA_hit_truth.pkl"), 'r')
    siRNA_hit_truth=pickle.load(f); f.close()
    
    f=open('../data/siRNA_simulated.pkl', 'r')
    normal_siRNAs=len(pickle.load(f)); f.close()
    
    accuracy_precision(normal_siRNAs, siRNA_hit_truth, siRNA_hit)
    
    return [(el, types[el]) for el in truth if el not in exp_hit], [(el, types[el]) for el in exp_hit if el not in truth]
    

def generateQCFile(num_plates=None, num_replicates=[1,2,3], rename=False,
                   folder = '../resultData/simulated_traj/simres'):
    
    dataFolder = os.path.join(folder, 'plates')
    plates=filter(lambda x: 'LT' in x, os.listdir(dataFolder))
    
    f=open('../data/qc_simulated.txt', 'w')
    f.write('Id\tsirnaId\ttype\tlabtekQC\tmanualSpotQC\tautoSpotQC\taccepted\tremark\n')
    siRNAid=1
    siRNA_hits=[]
    if num_plates is None:
        num_plates = range(1, len(plates)/3+1)
    
    for plate in num_plates:
        print plate, ' starting at siRNA ', siRNAid
        l=len(os.listdir(os.path.join(dataFolder, 'LT{:>04}_01'.format(plate))))
        if l <50:
            wells = [15, 26, 63, 74, 304, 315, 352]
        else:
            wells = range(1,385)
            file_=open(os.path.join(folder, 'annotations', 'annotation_LT{:>04}.pickle'.format(plate)))
            ann=pickle.load(file_); file_.close()
        for well in wells:
            if well in [15, 26, 63, 74, 304, 315, 352]:
                experiment_type='scrambled'
                s='103860'
            else:
                experiment_type='marker'
                s=siRNAid
                if ann[0][well] not in ['control', "normal"]:
                    siRNA_hits.append(s)
                siRNAid+=1
                
            for replicate in num_replicates:
                f.write('LT{:>04}_{:>02}--{:>03}\t{}\t{}\tpass\tpass\tpass\tTrue\tok\n'.format(plate, replicate, well, s, experiment_type) )                
                if rename:
                    try:
                        os.rename(os.path.join(dataFolder, 'LT{:>04}_{:>02}'.format(plate, replicate), 'hist_tabFeatures_{:>03}.pkl'.format(well) ), 
                                  os.path.join(dataFolder, 'LT{:>04}_{:>02}'.format(plate, replicate), 'hist_tabFeatures_{:>05}_01.pkl'.format(well)))
                    except:
                        print 'LT{:>04}_{:>02}'.format(plate, replicate), well
                
                
    f.close()
    
    print 'Maximum siRNA', siRNAid
    print 'Saving siRNA list in ', '../data/siRNA_simulated.pkl'
    f=open('../data/siRNA_simulated.pkl', 'w')
    pickle.dump(range(1, siRNAid), f); f.close()
    
    print 'Saving siRNA hit list in ', '../resultData/simulated_traj/simres/siRNA_hit_truth.pkl'
    f=open('../resultData/simulated_traj/simres/siRNA_hit_truth.pkl', 'w')
    pickle.dump(siRNA_hits, f); f.close()
    
    
    return 1
    

class BasicCountMachine(object):
    # transition_conditions = {source_state: (target_state, count_number)}
    def __init__(self, transition_conditions):
        self.transition_conditions = transition_conditions
        
    def get_states(self, length, start_state=None):
        if start_state is None:
            states = sorted(self.transition_conditions.keys())
            start_state = states[0]

        states = [start_state]
        in_state_count = 0
        current_state = start_state
        for abscount in range(length-1):
            in_state_count += 1
            if in_state_count >= self.transition_conditions[current_state][1]:
                current_state = self.transition_conditions[current_state][0]
                in_state_count = 0
            states.append(current_state)

        return np.array(states)

    def get_state_lengths(self, states):
        state_lengths = {}
        current_state = states[0]
        i = 0
        k = 0
        for s in states:
            if s == current_state:
                k += 1
            else:
                state_lengths[i] = {'state': current_state, 'l': k}
                k = 1
                i += 1
                current_state = s
        state_lengths[i] = {'state': current_state, 'l': k}
        
        return state_lengths
    
class PlateSimulator(object):
    def __init__(self, settings_filename=None, settings=None):
        if settings is None and settings_filename is None:
            raise ValueError("Either a settings object or a settings filename has to be given.")
        if not settings is None:
            self.settings = settings
        elif not settings_filename is None:
            self.settings = Settings(os.path.abspath(settings_filename), dctGlobals=globals())
        self.trajectory_simulator = TrajectorySimulator(settings=self.settings)
        self.movement_types = None
        
    def __call__(self, control_only=False):
        # sum should be 377 (7 negative controls)
        #parameter mean_variation_alpha for varying mean speed from one replicate to another
        if not control_only:
            hit_vec = {
                       'normal': 350, 
                       'directed1': 4,
                       'directed2': 8, 
                       'fast': 8,
                       'switch_fast': 7,
                       }
            res, annotation = self.simulate_plates_with_replicates(hit_vec, low_nb=100, high_nb=150, nb_min_thresh=20, 
                                                                   mean_variation_alpha=0.0)
            self.export_plate_simulation(res, annotation)
    
            hit_vec = {
                       'normal': 320, 
                       'directed1': 14,
                       'directed2': 6, 
                       'fast': 33,
                       'switch_fast': 4,
                       }        
            res, annotation = self.simulate_plates_with_replicates(hit_vec, low_nb=100, high_nb=150, nb_min_thresh=20, 
                                                                   mean_variation_alpha=0.5)
            self.export_plate_simulation(res, annotation)
    
            hit_vec = {
                       'normal': 357, 
                       'directed1': 1,
                       'directed2': 3, 
                       'fast': 4,
                       'switch_fast': 12,
                       }        
            res, annotation = self.simulate_plates_with_replicates(hit_vec, low_nb=100, high_nb=150, nb_min_thresh=20, 
                                                                   mean_variation_alpha=1.0)
            self.export_plate_simulation(res, annotation)
        
        for k in range(150):
            res, annotation = self.simulate_plates_with_replicates(control_only=True, low_nb=100, high_nb=150, nb_min_thresh=20, 
                                                               mean_variation_alpha=0.0)
            self.export_plate_simulation(res, annotation)
        for k in range(200):
            res, annotation = self.simulate_plates_with_replicates(control_only=True, low_nb=100, high_nb=150, nb_min_thresh=20, 
                                                               mean_variation_alpha=0.5)
            self.export_plate_simulation(res, annotation)
        for k in range(50):
            res, annotation = self.simulate_plates_with_replicates(control_only=True, low_nb=100, high_nb=150, nb_min_thresh=20, 
                                                               mean_variation_alpha=1)
            self.export_plate_simulation(res, annotation)

        return
    
    def get_standard_simulator_settings(self, mean_variation_alpha=0.0):

        simulator_settings = {

                              'movement_normal': {'N': 0,
                                                  'radius': {'distribution': 'Normal', 'mean': 0.26, 'stdev': 0.09},
                                                  'angle': {'distribution': 'Uniform', 'min': 0, 'max': 360},
                                                  },
                              'movement_fast': {'N': 0,
                                                'radius': {'distribution': 'Normal', 'mean': 0.4, 'stdev': 0.09},
                                                'angle': {'distribution': 'Uniform', 'min': 0, 'max': 360},
                                                },
                              'directional_bended': {'N': 0,
                                                     'angle_0': {'distribution': 'Uniform', 'min': 0, 'max': 360}, # start value
                                                     'angle': {'distribution': 'Normal_Moving_Mean', 'stdev': 25.0},
                                                     'angle_delta': {'distribution': 'Normal', 'mean': 50.0, 'stdev': 1.0},
                                                     'radius': {'distribution': 'Normal', 'mean': 0.26, 'stdev': 0.09},
                                                     },                               
                              'angle_switch': {# indicates directed movement in +- 180 degrees.
                                               'N': 0, 
                                               'radius': {'distribution': 'Normal', 'mean': 0.26, 'stdev': 0.09},
                                               'angle_0': {'distribution': 'Uniform', 'min': 0, 'max': 360}, # start angle
                                               'angle': {'distribution': 'Normal', 'stdev': 25.0, 
                                                         'switch': {'trans_count_0_1': {'min': 10, 'max': 15},
                                                                    'trans_count_1_0': {'min': 5, 'max': 12}}}
                                               },                                                                  
                              'mode_switch': {# mode_switch indicates the presence of 2 movements, 
                                              # one directed, the other random, with different speeds respectively.
                                              'N': 0, 
                                              'common_switch': {'trans_count_0_1': {'min': 8, 'max': 15},
                                                                'trans_count_1_0': {'min': 3, 'max': 5}},     
                                              'radius1': {'distribution': 'Normal', 'mean': 0.26, 'stdev': 0.09},
                                              'radius2': {'distribution': 'Normal', 'mean': 0.4, 'stdev': 0.09},                
                                              'angle1': {'distribution': 'Uniform', 'min': 0, 'max': 360},
                                              'angle2': {'distribution': 'Normal', 'stdev': 25.0},
                                              },
                              }
        self.movement_types = simulator_settings.keys()

        #pdb.set_trace()
        simulator_settings['movement_normal']['radius']['mean'] += mean_variation_alpha * np.random.randn() * 0.01 * simulator_settings['movement_normal']['radius']['stdev']
        simulator_settings['movement_fast']['radius']['mean'] += mean_variation_alpha * np.random.randn() * 0.01 * simulator_settings['movement_normal']['radius']['stdev']
        simulator_settings['directional_bended']['radius']['mean'] += mean_variation_alpha * np.random.randn() * 0.01 * simulator_settings['movement_normal']['radius']['stdev']
        simulator_settings['angle_switch']['radius']['mean'] += mean_variation_alpha * np.random.randn() * 0.01 * simulator_settings['movement_normal']['radius']['stdev']
        simulator_settings['mode_switch']['radius1']['mean'] += mean_variation_alpha * np.random.randn() * 0.01 * simulator_settings['movement_normal']['radius']['stdev']
        simulator_settings['mode_switch']['radius2']['mean'] += mean_variation_alpha * np.random.randn() * 0.01 * simulator_settings['movement_normal']['radius']['stdev']

        return simulator_settings
    
    def export_plate_simulation(self, res, annotation):
        simres_folder = os.path.join(self.settings.out_folder, 'simres')
        
        out_folder = os.path.join(simres_folder, 'plates')
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
                    
        # find the plates already simulated
        plates = list(set([x.split('_')[0] for x in 
                           filter(lambda y: y.rfind('_') >=0, os.listdir(out_folder))]))
        max_number = len(plates)
            
        # write annotation (hit categories and normal spots)
        annotation_folder = os.path.join(simres_folder, 'annotations')
        if not os.path.exists(annotation_folder):
            os.makedirs(annotation_folder)
        
        fp = open(os.path.join(annotation_folder, 'annotation_LT%04i.pickle' % (max_number +1)), 'w')
        pickle.dump(annotation, fp)
        fp.close()

        for plate_index in res.keys():
            plate_folder = os.path.join(out_folder, 'LT%04i_%02i' % ((max_number +1), (plate_index+1)))
            if not os.path.exists(plate_folder):
                print 'making %s' % plate_folder
                os.makedirs(plate_folder)
            for pos in res[plate_index].keys():
                output_struct = []
                k = 1
                for movement_type in res[plate_index][pos].keys():
                    for trajectory in res[plate_index][pos][movement_type]:
                        curr_traj={}
                        trajectory_length = len(trajectory)
                        offset = np.random.random_integers(1, 100-trajectory_length)
                        for j in range(trajectory_length): 
                            curr_traj[(offset + j, k)] = trajectory[j]
                        output_struct.append(curr_traj)
                # writes one pickle for each position
                fp = open(os.path.join(plate_folder, 'LT%04i_%02i--%03i.pickle' % ((max_number+1), (plate_index+1), pos)), 'w')
                pickle.dump(output_struct, fp)
                fp.close()
                
        return
    
    def simulate_plates_with_replicates(self, hit_vec=None,
                                        low_nb=100, high_nb=200, nb_min_thresh=20, nb_type='from_file', 
                                        mean_variation_alpha=0.0,
                                        control_only=False):
            
        res = {}
    #In particular those are the same negative controls as for the mitocheck plates
        neg_controls = [15, 26, 63, 74, 304, 315, 352]
        exp_spots_o = filter(lambda x: x not in neg_controls, range(1, 385))
        exp_spots = np.random.permutation(exp_spots_o)
        
        #annotation = dict(zip(range(1, 385), ['neg']))
        annotation = {}   
        for plate_index in range(3): 
            res[plate_index] = {}
            annotation[plate_index] = {}
            # this corresponds to a reset of simulator settings with N=1 for all movement types
            # This also sets the self.movement_types
            simulator_settings = self.get_standard_simulator_settings(mean_variation_alpha=mean_variation_alpha)
            
            if nb_type == 'uniform': 
                # the total number of trajectories is given by the uniform distribution
                # between low_nb and high_nb.            
                total_nb_traj = dict(zip(neg_controls + exp_spots.tolist(), 
                                         np.random.random_integers(low=low_nb, high=high_nb, size=384)))
            elif nb_type=='normal':
                # the total number of trajectories is given by the normal distribution
                # with mean uniformly distributed between low_nb and high_nb
                mean_nb = np.random.random_integers(low=low_nb, high=high_nb)
                sigma_nb = np.random.random_integers(10, 40)
                total_nb_traj = dict(zip(neg_controls + exp_spots.tolist(),
                                         np.maximum(sigma_nb * np.random.randn(384) + mean_nb, nb_min_thresh).astype(int)))
            elif nb_type=='from_file':
                fp = open(self.settings.nb_traj_file, 'r')
                temp = pickle.load(fp)
                fp.close()
                default_values = np.array(temp)[np.random.permutation(len(temp))][:384]
#COMMENTE EN ATTENDANT D'AVOIR LE FICHIER DE THOMAS
#                # choose a plate at random
#                real_plates = temp.keys()
#                nb_pos = 0
#                while(nb_pos < 50): 
#                    chosen_plate = real_plates[np.random.random_integers(len(temp.keys()))]
#                    chosen_positions = sorted(temp[chosen_plate].keys())
#                    nb_pos = len(chosen_positions)
#                    print 'read plate : %s %i' % (chosen_plate, len(chosen_positions))
#
#                values = dict(zip(chosen_positions, 
#                                  [temp[chosen_plate][x] for x in chosen_positions]))
#                
#                maxval = max(values.values())
#                minval = min(values.values())
#                default_values = np.random.random_integers(low=minval, high=maxval, size=384)
#                default_values[np.array(chosen_positions)-1] = [temp[chosen_plate][x] for x in chosen_positions]
                total_nb_traj = dict(zip(range(1, 385), default_values))
                
            # simulation for the negative controls
            for pos in neg_controls:          
                # get the counts for the different movement types
                counts = self.get_trajectory_counts('normal', total_nb_traj[pos])
    
                # assign for each movement_type the correct number of trajectories
                for mt in self.movement_types:
                    simulator_settings[mt]['N'] = counts[mt]
                res[plate_index][pos] = self.trajectory_simulator.make_trajectories(simulator_settings)
                annotation[plate_index][pos] = 'control'
            if control_only:
                continue
            # simulation for other spots
            k=0
            for exp_type in hit_vec:
                # nb_spots is the number of experiments for this hit category
                nb_spots = hit_vec[exp_type]
    
                # find the positions (at random without replacement) for this hit category.
                positions = exp_spots[k:min(k+nb_spots, len(exp_spots))]
                k = k + nb_spots
                
                for pos in positions:                 
                    # reset of the settings
                    simulator_settings = self.get_standard_simulator_settings(mean_variation_alpha=mean_variation_alpha)
    
                    counts = self.get_trajectory_counts(exp_type, total_nb_traj[pos])
                    # assign for each movement_type the correct number of trajectories
                    for mt in self.movement_types:
                        simulator_settings[mt]['N'] = counts[mt]
                    res[plate_index][pos] = self.trajectory_simulator.make_trajectories(simulator_settings)
                    annotation[plate_index][pos] = exp_type

        return res, annotation
    
    def get_trajectory_counts(self, exp_type, N):
        counts = {}
        pvec = self.settings.probabilities[exp_type]
        movement_types = pvec.keys()
        pvals = [pvec[x] for x in movement_types]
        
        # the counts are estimated from a multinomial distribtion
        # we draw N times according to the probabilities given in the settings file.
        counts = dict(zip(movement_types, np.random.multinomial(N, pvals)))
        
        return counts
    

class RealData(object):
    def __init__(self):
        self.data_folder = '/share/data20T/mitocheck/tracking_results'
        self.featuresSaved=[
               #i. features that i'll use in the analysis
               'ball number1',
               'ball number2',
               'norm convex hull area',
               'corrected straightness index',
               'diffusion coefficient',
               'effective space length', 
               'effective speed',
               'entropy1', 
               'entropy2',#rayons 0, 14 ci-dessus
               'largest move',
               'mean squared displacement',
               'movement type',
               'signed turning angle',
               'diffusion adequation',
      
               #ii.features that will not be used in the analysis
               'mvt type adequation',
               'mean acceleration',  
               'mean displacements',  
               'mean persistence',  
               'mean straight',  
               'mean turningangle',
               'convex hull area',
               'correlation speed and persistence',
               'slope SP',
               'total space length']
        return
    
    def get_plates(self):
        
        plates = filter(lambda x: len(x.split('--')[0].split('_')) >= 2 and int(x.split('--')[0].split('_')[0][2:]) < 200, 
                        os.listdir(self.data_folder))
        return plates


    def get_number_trajectories(self):
        
        plates = self.get_plates()
        
        nb_trajectories = {}
        
        for pl in plates:
            start_time = time.time()
            nb_trajectories[pl] = {}
            plate_folder = os.path.join(self.data_folder, pl)
            feature_files = filter(lambda x: x[:len('hist_tabFeatures')]=='hist_tabFeatures', 
                                   os.listdir(plate_folder))
            for filename in feature_files:
                #hist2_tabFeatures_00026_01
                pos = int(filename.split('_')[2])
                fp = open(os.path.join(plate_folder, filename), 'r')
                X = pickle.load(fp)
                fp.close()
                nb_trajectories[pl][pos] = len(X[1])
            diff_time = time.time() - start_time
            print 'elapsed time: %02i:%02i:%02i' % ((diff_time / 3600), ((diff_time % 3600) / 60), (diff_time % 60))
            
        return nb_trajectories
    
    def write_number_of_trajectories(self, filename):
        nb_trajectories = self.get_number_trajectories()
        fp = open(filename, 'w')
        pickle.dump(nb_trajectories, fp)
        fp.close()
        return

    def get_avg_feature_distribution(self, feature_name):
        if feature_name in self.featuresSaved:
            feature_index = self.featuresSaved.index(feature_name)
        else:
            raise ValueError("feature not calculated")

        plates = self.get_plates()
        features = {}
        for pl in plates:
            features[pl] = {}
            plate_folder = os.path.join(self.data_folder, pl)
            feature_files = filter(lambda x: x[:len('hist_tabFeatures')]=='hist_tabFeatures', 
                                   os.listdir(plate_folder))
            for filename in feature_files:
                pos = int(filename.split('_')[2])
                fp = open(os.path.join(plate_folder, filename), 'r')
                temp = pickle.load(fp)
                fp.close()
                X = temp[0]
                features[pl][pos] = np.mean(X[:,feature_index])

        return features

    def get_full_feature_distribution(self, feature_name):
        if feature_name in self.featuresSaved:
            feature_index = self.featuresSaved.index(feature_name)
        else:
            raise ValueError("feature not calculated")

        plates = self.get_plates()
        features = {}
        for pl in plates:
            features[pl] = {}
            plate_folder = os.path.join(self.data_folder, pl)
            feature_files = filter(lambda x: x[:len('hist_tabFeatures')]=='hist_tabFeatures', 
                                   os.listdir(plate_folder))
            for filename in feature_files:
                pos = int(filename.split('_')[2])
                fp = open(os.path.join(plate_folder, filename), 'r')
                temp = pickle.load(fp)
                fp.close()
                X = temp[0]
                features[pl][pos] = X[:,feature_index]

        return features
    
    def write_feature_avg(self, feature_name, out_folder):
        features = self.get_avg_feature_distribution(feature_name)
        filename = os.path.join(out_folder, '%s.pickle' % feature_name)
        fp = open(filename, 'w')
        pickle.dump(features, fp)
        fp.close()
        return

    def write_feature_full(self, feature_name, out_folder):
        features = self.get_full_feature_distribution(feature_name)
        filename = os.path.join(out_folder, 'full_%s.pickle' % feature_name)
        fp = open(filename, 'w')
        pickle.dump(features, fp)
        fp.close()
        return

    def analyze_nb_traj(self):
        fp = open('/Users/twalter/data/migration_sim/data/real/nb_tracks.pickle', 'r')
        nb_data = pickle.load(fp)
        fp.close()

        ####################### number of spots
        nb_of_spots = [len(nb_data[plate]) for plate in nb_data.keys()]
        print 'total: %i' % sum(nb_of_spots)
        print 'average: %f' % (sum(nb_of_spots) / float(len(nb_of_spots)))
        
        for plate in nb_data:
            print '%s\t%i' % (plate[:9], len(nb_data[plate]))
            
        return
    
    def analyze_mean_displacement(self):
        fp = open('/Users/twalter/data/migration_sim/data/real/mean displacements.pickle', 'r')
        featurevals = pickle.load(fp)
        fp.close()

        # plot histogram
        Xmat = []
        for plate in featurevals.keys():
            Xmat.extend(featurevals[plate][x] for x in featurevals[plate].keys())
        X = np.array(Xmat)    
        
        plot_folder = '/Users/twalter/data/migration_sim/data/plots'
        
        h = plotter_stats.Histogram()
        filename = os.path.join(plot_folder, 'histo', 'mean_displacement.png')
        h(X, filename, colorvec='blue',
          title='mean displacement (spot means)',
          bins=60)

        print 'mean displacement: mean=%f\t stdev=%f' % (np.mean(X), np.std(X))
        print 'max displacement: %f' % np.max(X)
        print 'min displacement: %f' % np.min(X)
        p = np.percentile(X, [50, 90, 95, 99])
        print 'percentiles: 50: %f\t90: %f\t95: %f\t99: %f' % tuple(p)

        ###############################################################
        ###################### FULL DISTRIBUTION ######################
        ###############################################################
        
        fp = open('/Users/twalter/data/migration_sim/data/real/full_mean displacements.pickle', 'r')
        featurevals = pickle.load(fp)
        fp.close()

        # full mean displacement
        # plot histogram
        Xmat = []
        for plate in featurevals.keys():
            for x in featurevals[plate].keys(): 
                Xmat.extend(featurevals[plate][x])
        X = np.array(Xmat)
        print X.shape
        
        plot_folder = '/Users/twalter/data/migration_sim/data/plots'
        
        h = plotter_stats.Histogram()
        filename = os.path.join(plot_folder, 'histo', 'full_mean_displacement.png')
        print
        print '*******************************************'
        print 'full mean displacement: mean=%f\t stdev=%f' % (np.mean(X), np.std(X))
        print 'full max displacement: %f' % np.max(X)
        print 'full min displacement: %f' % np.min(X)
        p = np.percentile(X, [50, 90, 95, 99])
        print 'full percentiles: 50: %f\t90: %f\t95: %f\t99: %f' % tuple(p)

        h(X, filename, colorvec='blue', 
          title='mean displacement (single cell values)',
          bins=200, xlim=(0, p[-1]))


        return
    
    def __call__(self, feature, feature_stat='avg'):
        out_folder = '/cbio/donnees/twalter/data/migration/features_tracks'
        if feature == 'nb_traj':
            filename = os.path.join(out_folder, 'nb_tracks.pickle')
            self.write_number_of_trajectories(filename)
        elif feature_stat=='avg': 
            self.write_feature_avg(feature, out_folder)
        else:
            self.write_feature_full(feature, out_folder)

        print 'DONE'
        return
    
            
class TrajectorySimulator(object):
    def __init__(self, settings_filename=None, settings=None):
        if settings is None and settings_filename is None:
            raise ValueError("Either a settings object or a settings filename has to be given.")
        if not settings is None:
            self.settings = settings
        elif not settings_filename is None:
            self.settings = Settings(os.path.abspath(settings_filename), dctGlobals=globals())

    def rand(self, alpha):        
        return 2*np.pi * alpha / 360
    
    def make_trajectories(self, simulator_settings=None, 
                          traj_length_from_file=True):

        if simulator_settings is None:
            simumlator_settings = self.settings.simulator_settings
            
        trajectories = {}
        for movement_type in simulator_settings:
            
            trajectories[movement_type] = []
            
            # number of trajectories
            N =  simulator_settings[movement_type]['N']

            if N<=0:
                continue
            
            # trajectory length based on a realistic lengths sample
            if traj_length_from_file:
                f=open(self.settings.length_file)
                lengths=pickle.load(f); f.close()
                seeds=np.random.mtrand.permutation(len(lengths))[:N]
                lvec=lengths[seeds]
            else:
                if self.settings.L['distribution'] == 'Normal': 
                    lvec = np.random.normal(self.settings.L['mean'], self.settings.L['stdev'], N).astype(int)
                elif self.settings.L['distribution'] == 'Uniform': 
                    lvec = np.random.uniform(self.settings.L['min'], self.settings.L['max'], N).astype(int)

            for i,l in enumerate(lvec):
                l = l-1

                if movement_type == 'mode_switch':
                    transition_count_0_1 = np.random.randint(simulator_settings[movement_type]['common_switch']['trans_count_0_1']['min'],
                                                             simulator_settings[movement_type]['common_switch']['trans_count_0_1']['max'])
                    transition_count_1_0 = np.random.randint(simulator_settings[movement_type]['common_switch']['trans_count_1_0']['min'],
                                                             simulator_settings[movement_type]['common_switch']['trans_count_1_0']['max'])
                    tc = {0: (1, transition_count_0_1),
                          1: (0, transition_count_1_0), }
                    bc = BasicCountMachine(tc)
                    states = bc.get_states(l)
                    state_lengths = bc.get_state_lengths(states)
                    angle_start = np.random.randint(0,360)
                    rvec = []
                    phivec = []
                    for k in sorted(state_lengths.keys()):
                        state = state_lengths[k]['state']
                        sl = state_lengths[k]['l']
                        if state == 0:
                            r = np.abs(np.random.normal(simulator_settings[movement_type]['radius1']['mean'],simulator_settings[movement_type]['radius1']['stdev'], sl))
                            phi = np.random.uniform(simulator_settings[movement_type]['angle1']['min'],
                                                    simulator_settings[movement_type]['angle1']['max'], sl)
                            angle_start = phi[-1]
                        else:
                            r = np.abs(np.random.normal(simulator_settings[movement_type]['radius2']['mean'], simulator_settings[movement_type]['radius2']['stdev'], sl))
                            phi = np.random.normal(angle_start, simulator_settings[movement_type]['angle2']['stdev'], sl)
                        rvec.extend(r)
                        phivec.extend(phi)

                else:                    
                    # draw radius
                    if simulator_settings[movement_type]['radius']['distribution'] == 'Normal':
                        
                        rvec = np.abs(np.random.normal(simulator_settings[movement_type]['radius']['mean'],
                                                       simulator_settings[movement_type]['radius']['stdev'],
                                                       l))
                    elif simulator_settings[movement_type]['radius']['distribution'] == 'Normals_two_states': 
                         
                        r1 = np.abs(np.random.normal(simulator_settings[movement_type]['radius']['radius1']['mean'],
                                                     simulator_settings[movement_type]['radius']['radius1']['stdev'], l))
                        r2 = np.abs(np.random.normal(simulator_settings[movement_type]['radius']['radius2']['mean'],
                                                     simulator_settings[movement_type]['radius']['radius2']['stdev'], l))
                        transition_count_0_1 = np.random.randint(simulator_settings[movement_type]['radius']['trans_count_0_1']['min'],
                                                                 simulator_settings[movement_type]['radius']['trans_count_0_1']['max'])
                        transition_count_1_0 = np.random.randint(simulator_settings[movement_type]['radius']['trans_count_1_0']['min'],
                                                                 simulator_settings[movement_type]['radius']['trans_count_1_0']['max'])

                        tc = {0: (1, transition_count_0_1),
                              1: (0, transition_count_1_0), }
                        bc = BasicCountMachine(tc)
                        states = bc.get_states(l)
                        rvec = states * r1 + (1-states) * r2
                        
	                # draw angle mean (one mean for the entire trajectory)
                    if 'angle_0' in simulator_settings[movement_type]:
                        if simulator_settings[movement_type]['angle_0']['distribution'] == 'Uniform':  
                            angle_mean = np.random.uniform(simulator_settings[movement_type]['angle_0']['min'],
                                                           simulator_settings[movement_type]['angle_0']['max'], 1)[0]
                    
                    # draw angle
                    if simulator_settings[movement_type]['angle']['distribution'] == 'Uniform':  
                        phivec = np.random.uniform(simulator_settings[movement_type]['angle']['min'],
                                                   simulator_settings[movement_type]['angle']['max'], l)
                    elif simulator_settings[movement_type]['angle']['distribution'] == 'Normal':
                        phivec = np.random.normal(angle_mean,
                                                  simulator_settings[movement_type]['angle']['stdev'], l)
                    elif simulator_settings[movement_type]['angle']['distribution'] == 'Normal_Moving_Mean':
                        deltaphivec = np.random.normal(simulator_settings[movement_type]['angle_delta']['mean'],
                                                       simulator_settings[movement_type]['angle_delta']['stdev'], l)                                        
                        phivec = np.random.normal(0,
                                                  simulator_settings[movement_type]['angle']['stdev'], l)                    
                        phivec += (2 * np.random.binomial(1, .5) - 1) * np.cumsum(deltaphivec)
                    else:
                        raise ValueError("Unknown angle distribution")
                    
                    # make a switch (180 degree turn)
                    if 'switch' in simulator_settings[movement_type]['angle']:
                        transition_count_0_1 = np.random.randint(simulator_settings[movement_type]['angle']['switch']['trans_count_0_1']['min'],
                                                                 simulator_settings[movement_type]['angle']['switch']['trans_count_0_1']['max'])
                        transition_count_1_0 = np.random.randint(simulator_settings[movement_type]['angle']['switch']['trans_count_1_0']['min'],
                                                                 simulator_settings[movement_type]['angle']['switch']['trans_count_1_0']['max'])
    
                        tc = {0: (1, transition_count_0_1),
                              1: (0, transition_count_1_0)}
                        
                        bc = BasicCountMachine(tc)
                        states = bc.get_states(l)
                        phivec += 180 * states
	                    
	           	                    
                # transform angle (to radiant and between 0 and 2 PI)
                for i in range(len(phivec)):
                    phivec[i] = self.rand(phivec[i])
                    while phivec[i] < 0:
                        phivec[i] += 2*np.pi
                    while phivec[i] > 2*np.pi:
                        phivec[i] -= 2*np.pi
            
                # make trajectory
                delta_x = [r*np.cos(phi) for r, phi in zip(rvec, phivec)]
                delta_y = [r*np.sin(phi) for r, phi in zip(rvec, phivec)]
                xvec = np.cumsum([0] + delta_x) # we assume that all trajectories start at (0, 0)
                yvec = np.cumsum([0] + delta_y) 
                trajectory = zip(xvec, yvec)
                trajectories[movement_type].append(trajectory)

        return trajectories

    def export_trajectories(self, trajectories, filename):
        filename = os.path.join(self.settings.out_folder, '%s.pkl' % filename)
        plate="PL0"; well="WO"
        dict_={plate:{well:ensTraj()}}
        cour = dict_[plate][well]
        num=0
        longueurFilm=0
        for movement_type in trajectories:
            print movement_type
            for i in range(len(trajectories[movement_type])):
                print i,
                traj=trajectoire(num, id=0)
                frame=0
                for x,y in trajectories[movement_type][i]:
                    traj.ajoutPoint(x,y, frame, 0)
                    frame+=1
                traj.density=[4 for k in range(frame)]
                cour.add(traj)
                num+=1
                longueurFilm=max(longueurFilm, frame)
        movie_length={plate:{well:longueurFilm}}
        connexions= {plate: {well: {}}}
        fp = open(filename, 'w')
        pickle.dump(dict(zip(['tracklets dictionary', 'connexions between tracklets', 'movie_length'], [dict_, connexions, movie_length])), fp)
        fp.close()
        return

    def plot_trajectories(self, trajectories):
        adapt_axis = self.settings.adapt_axis
        grid = self.settings.grid
        
        if adapt_axis:
            maxx = max([max([max([x[0] for x in trajectories[movement_type][i]]) for i in range(len(trajectories[movement_type]))])
                        for movement_type in trajectories])
            minx = min([min([min([x[0] for x in trajectories[movement_type][i]]) for i in range(len(trajectories[movement_type]))])
                        for movement_type in trajectories])
            maxy = max([max([max([x[1] for x in trajectories[movement_type][i]]) for i in range(len(trajectories[movement_type]))])
                        for movement_type in trajectories])
            miny = min([min([min([x[1] for x in trajectories[movement_type][i]]) for i in range(len(trajectories[movement_type]))])
                        for movement_type in trajectories])
    
        for movement_type in trajectories:
            for i in range(len(trajectories[movement_type])):
                filename = os.path.join(self.settings.plot_folder, '%s--%i.png' % (movement_type, i))
                xvec = [x[0] for x in trajectories[movement_type][i]]
                yvec = [x[1] for x in trajectories[movement_type][i]]

                fig = plt.figure(1, figsize=(20,10))
                ax = plt.subplot(1,1,1)
                ax.plot(xvec, yvec, color="blue", linewidth=1.0, marker='o')
                for k in range(len(xvec)):
                    ax.text(xvec[k], yvec[k], '{}'.format(k))

                plt.title('Trajectory %i for movement type %s (L: %i)' % (i, movement_type, len(xvec)))
                plt.xlabel('x')
                plt.ylabel('y')

                if adapt_axis:
                    ax.axis([minx, maxx, miny, maxy])
                      
                if grid:
                    ax.grid(b=True, which='major', linewidth=1.5)
    
                fig.savefig(filename)
                plt.close(1)

        return
    
#    def extractFeatures(self, filename, verbose):
#        fp = open(os.path.join(self.settings.out_folder, "{}.pkl".format(filename)))
#        data=pickle.load(fp); fp.close()
#        d,c, movie_length = data['tracklets dictionary'], data['connexions between tracklets'], data['movie_length']
#        res = histogramPreparationFromTracklets(d, c, self.settings.out_folder, True, verbose, movie_length)
        
    def __call__(self, filename, verbose):
        print 'doing trajectories'
        trajectories = self.make_trajectories()
        print 'exporting trajectories'
        self.export_trajectories(trajectories, filename)
        self.plot_trajectories(trajectories)
        print 'extracting trajectories features'
        self.extractFeatures(filename, verbose)
        return
    
    #SCRIPT UTILISE POUR FAIRE LES FEATURES
#    simulateur = TrajectorySimulator(settings_filename='simulator/simulation_settings_2014_01_16.py')
#    simulateur('simulated_trajectories', 0)
#
#    f=open('../resultData/simulated_traj/hist_tabFeatures_WO.pkl')
#    tabFeatures, _, histNC, _, _, _ = pickle.load(f)
#    r2 = histLogTrsforming(tabFeatures)  
#    print r2.shape, 'not normalized'
#    r2=r2[:,:-1]; r2=(r2-np.mean(r2,0))/np.std(r2,0)
#    print 'computing bins'
#    histogrammeMatrix = computingBins(histNC, mat_hist_sizes[0])
#    print 'saving'
#    f=open('../resultData/simulated_traj/histogramsNtotSim.pkl', 'w')
#    pickle.dump(np.hstack((r2, histogrammeMatrix)), f); f.close()

                

            
                 
