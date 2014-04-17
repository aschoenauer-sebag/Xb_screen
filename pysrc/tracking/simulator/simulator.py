import os, re, sys, time
import cPickle as pickle
import numpy as np

#import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
from PyPack.fHacktrack2 import ensTraj, trajectoire
from trajPack.trajFeatures import histogramPreparationFromTracklets
from util.settings import Settings

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
    
    def make_trajectories(self):

        trajectories = {}
        for movement_type in self.settings.simulator_settings:
            trajectories[movement_type] = []
            
            # number of trajectories
            N =  self.settings.simulator_settings[movement_type]['N']

            # trajectory length based on a realistic lengths sample
            f=open(self.settings.length_file)
            lengths=pickle.load(f); f.close()
            seeds=np.random.mtrand.permutation(len(lengths))[:N]
            lvec=lengths[seeds]
#            if self.settings.L['distribution'] == 'Normal': 
#                lvec = np.random.normal(self.settings.L['mean'], self.settings.L['stdev'], N).astype(int)
#            elif self.settings.L['distribution'] == 'Uniform': 
#                lvec = np.random.uniform(self.settings.L['min'], self.settings.L['max'], N).astype(int)

            for i,l in enumerate(lvec):
                l = l-1

                if movement_type == 'mode_switch':
                    transition_count_0_1 = np.random.randint(self.settings.simulator_settings[movement_type]['common_switch']['trans_count_0_1']['min'],
                                                             self.settings.simulator_settings[movement_type]['common_switch']['trans_count_0_1']['max'])
                    transition_count_1_0 = np.random.randint(self.settings.simulator_settings[movement_type]['common_switch']['trans_count_1_0']['min'],
                                                             self.settings.simulator_settings[movement_type]['common_switch']['trans_count_1_0']['max'])
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
                            r = np.random.normal(self.settings.simulator_settings[movement_type]['radius1']['mean'],
                                                 self.settings.simulator_settings[movement_type]['radius1']['stdev'], sl)                        
                            phi = np.random.uniform(self.settings.simulator_settings[movement_type]['angle1']['min'],
                                                    self.settings.simulator_settings[movement_type]['angle1']['max'], sl)
                            angle_start = phi[-1]
                        else:
                            r = np.random.normal(self.settings.simulator_settings[movement_type]['radius2']['mean'],
                                                 self.settings.simulator_settings[movement_type]['radius2']['stdev'], sl)                        
                            phi = np.random.normal(angle_start,
                                                   self.settings.simulator_settings[movement_type]['angle2']['stdev'], sl)
                        rvec.extend(r)
                        phivec.extend(phi)
                        #print i, l, len(rvec), len(phivec)
                else:
                
	                # draw radius
	                if self.settings.simulator_settings[movement_type]['radius']['distribution'] == 'Normal':
	                    rvec = np.random.normal(self.settings.simulator_settings[movement_type]['radius']['mean'],
	                                            self.settings.simulator_settings[movement_type]['radius']['stdev'],
	                                            l)
	                elif self.settings.simulator_settings[movement_type]['radius']['distribution'] == 'Normals_two_states':
	                    
	                    # r1
	                    r1 = np.random.normal(self.settings.simulator_settings[movement_type]['radius']['radius1']['mean'],
	                                          self.settings.simulator_settings[movement_type]['radius']['radius1']['stdev'], l)
	                    r2 = np.random.normal(self.settings.simulator_settings[movement_type]['radius']['radius2']['mean'],
	                                          self.settings.simulator_settings[movement_type]['radius']['radius2']['stdev'], l)
	
	                    transition_count_0_1 = np.random.randint(self.settings.simulator_settings[movement_type]['radius']['trans_count_0_1']['min'],
	                                                             self.settings.simulator_settings[movement_type]['radius']['trans_count_0_1']['max'])
	                    transition_count_1_0 = np.random.randint(self.settings.simulator_settings[movement_type]['radius']['trans_count_1_0']['min'],
	                                                             self.settings.simulator_settings[movement_type]['radius']['trans_count_1_0']['max'])
	
	                    tc = {0: (1, transition_count_0_1),
	                          1: (0, transition_count_1_0), }
	                    
	                    bc = BasicCountMachine(tc)
	                    states = bc.get_states(l)
	                    rvec = states * r1 + (1-states) * r2
	                    
	                # draw angle mean (one mean for the entire trajectory)
	                if 'angle_0' in self.settings.simulator_settings[movement_type]:
	                    if self.settings.simulator_settings[movement_type]['angle_0']['distribution'] == 'Uniform':  
	                        angle_mean = np.random.uniform(self.settings.simulator_settings[movement_type]['angle_0']['min'],
	                                                       self.settings.simulator_settings[movement_type]['angle_0']['max'], 1)[0]
	                
	                # draw angle
	                if self.settings.simulator_settings[movement_type]['angle']['distribution'] == 'Uniform':  
	                    phivec = np.random.uniform(self.settings.simulator_settings[movement_type]['angle']['min'],
	                                               self.settings.simulator_settings[movement_type]['angle']['max'], l)
	                elif self.settings.simulator_settings[movement_type]['angle']['distribution'] == 'Normal':
	                    phivec = np.random.normal(angle_mean,
	                                              self.settings.simulator_settings[movement_type]['angle']['stdev'], l)
	                elif self.settings.simulator_settings[movement_type]['angle']['distribution'] == 'Normal_Moving_Mean':
	                    deltaphivec = np.random.normal(self.settings.simulator_settings[movement_type]['angle_delta']['mean'],
	                                                   self.settings.simulator_settings[movement_type]['angle_delta']['stdev'], l)                                        
	                    phivec = np.random.normal(0,
	                                              self.settings.simulator_settings[movement_type]['angle']['stdev'], l)                    
	                    phivec += (2 * np.random.binomial(1, .5) - 1) * np.cumsum(deltaphivec)
	                else:
	                    raise ValueError("Unknown angle distribution")
	                
	                # make a switch (180 degree turn)
	                if 'switch' in self.settings.simulator_settings[movement_type]['angle']:
	                    transition_count_0_1 = np.random.randint(self.settings.simulator_settings[movement_type]['angle']['switch']['trans_count_0_1']['min'],
	                                                             self.settings.simulator_settings[movement_type]['angle']['switch']['trans_count_0_1']['max'])
	                    transition_count_1_0 = np.random.randint(self.settings.simulator_settings[movement_type]['angle']['switch']['trans_count_1_0']['min'],
	                                                             self.settings.simulator_settings[movement_type]['angle']['switch']['trans_count_1_0']['max'])
	
	                    tc = {0: (1, transition_count_0_1),
	                          1: (0, transition_count_1_0)}
	                    
	                    bc = BasicCountMachine(tc)
	                    states = bc.get_states(l)
	                    phivec += 180 * states
	                    
	                # 
	                    
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
    
    def extractFeatures(self, filename, verbose):
        fp = open(os.path.join(self.settings.out_folder, "{}.pkl".format(filename)))
        data=pickle.load(fp); fp.close()
        d,c, movie_length = data['tracklets dictionary'], data['connexions between tracklets'], data['movie_length']
        res = histogramPreparationFromTracklets(d, c, self.settings.out_folder, True, verbose, movie_length)
        
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

                

            
                 
