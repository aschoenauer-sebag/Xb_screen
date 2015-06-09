import os
import cPickle as pickle
from optparse import OptionParser
import sys, pdb
import numpy as np

from tracking.importPack.imp import importRawSegFromHDF5, frameLots

from util import settings

class TrackPrediction(object):
    def __init__(self, plate, w, settings,# loadingFolder, dataFolder, outputFolder, training = False, first=True, 
                 new_cecog_files=False):#, intensity_qc_dict=None, separating_function=None):
        self.plate=plate
        self.well = w
        
        self.settings=settings
        self.new_cecog_files = new_cecog_files
        
        self.dataFolder=os.path.join(self.settings.dataFolder, plate, 'hdf5')
        
    def _getIntensityQC(self):
        if 'LT' in self.plate:
            intensity_qc_dict=None
        else:
            w=int(self.well.split('_')[0])
            print "Loading manual and out of focus qc files"
            f=open('../data/xb_manual_qc.pkl', 'r')
            d1=pickle.load(f); f.close()
            if self.plate in d1 and w in d1[self.plate]:
                print "This well failed the manual quality control."
                return
            f=open('../data/xb_focus_qc.pkl', 'r')
            d2=pickle.load(f); f.close()
            if self.plate in d2 and w in d2[self.plate]:
                print "This well failed the focus/cell count quality control."
                return
            
            print 'Loading intensity qc file'
            intensity_qc_file=open('../data/xb_intensity_qc.pkl', 'r')
            intensity_qc_dict=pickle.load(intensity_qc_file); intensity_qc_file.close()    
            intensity_qc_dict=intensity_qc_dict[self.plate] if self.plate in intensity_qc_dict else None
            
        return intensity_qc_dict
    
    def _cleanNaNs(self,newFrameLot, feature_deleted_in_model, tabF, first):
        '''
        Checking which features are NaNs
        '''
        #global FEATURE_NUMBER
        features_to_del = np.where(np.isnan(tabF))[1]
        
        toZeros = filter(lambda x: x not in feature_deleted_in_model, features_to_del)
        if toZeros !=[]:
            sys.stderr.write("Attention attention, plate {}, some features here have NaN entries, and this was not the case in the training set. They are put to 0".format(self.plate))
            newFrameLot.zeros(toZeros)
    
        newFrameLot.clean(feature_deleted_in_model)
        if first:
            FEATURE_NUMBER =tabF.shape[1]-len(feature_deleted_in_model)
            
        return newFrameLot, FEATURE_NUMBER
    
    def _gettingRaw(self, filename, filenameT, well, secondary=False,name_primary_channel='primary__primary3', frames_to_skip=None):
        print "images loading : plate = "+self.plate+",well = "+well    
        tabF = None
        try:
            frameLotC, tabF = importRawSegFromHDF5(filename, self.plate, well, secondary=secondary, name_primary_channel=name_primary_channel, frames_to_skip=frames_to_skip,
                                                       separating_function=self.settings.separating_function)
        except ValueError:
            sys.stderr.write( sys.exc_info()[1])
            sys.stderr.write("File {} containing data for plate {}, well {} does not contain all necessary data".format(filename, self.plate, well))
            return None, None
    
        if filenameT is not None:
            print "training set loading : filename ="+filenameT
            frameLotC.addTraining2(filename, filenameT)
            print "current training set content :"
            print frameLotC.statisticsTraining2()
        
        return frameLotC, tabF

    def gettingSolu(self,first=True):
        tabF = None

        #loading features to delete    
        fichier = open(os.path.join(self.settings.loadingFolder,"featuresToDelete.pkl"), 'r')
        features_deleted_in_model = pickle.load(fichier)
        fichier.close()
        
        newFrameLot = None
        
        #loading information from h5 file
        filename = os.path.join(self.dataFolder, self.well)
        name_primary_channel='primary__primary' if not self.new_cecog_files else 'primary__primary3'
        well=self.well.split('.')[0]    
        
        #preparing training set loading if necessary
        if self.settings.training:
            filenameT = '/media/lalil0u/New/workspace2/Tracking/data/trainingset/PL'+self.plate+"___P"+well+"___T00000.xml"
        else:
            filenameT = None
            
        #preparing loading QC if there is one on a per frame basis
        intensity_qc_dict = self._getIntensityQC()
        
        if intensity_qc_dict is not None and int(well) in intensity_qc_dict and intensity_qc_dict[int(well)].shape[0]>0:
            frames_to_skip=intensity_qc_dict[int(well)]
            for el in frames_to_skip:
                if el+1 in frames_to_skip:
                    ind=np.where(frames_to_skip==el)
                    print "Issue with {}, {} in frames_to_skip".format(el, el+1)
                    pdb.set_trace()
                    #frames_to_skip=np.delete(frames_to_skip, [ind, ind+1])
        else:
            frames_to_skip=None
            
        #Now getting framelots
        frameLotC, tabFC = self._gettingRaw(filename, filenameT, well, name_primary_channel=name_primary_channel, frames_to_skip=frames_to_skip)
         
        if newFrameLot == None:
            newFrameLot = frameLotC 
        else: newFrameLot.addFrameLot(frameLotC)
        tabF = tabFC if tabF == None else np.vstack((tabF, tabFC))
    
        #cleaning from NaNs
        return self._cleanNaNs(newFrameLot, features_deleted_in_model, tabF, first)
        

if __name__ == '__main__':
    verbose=0
    description =\
'''
%prog - Parallel cell tracking
Performs trajectory extraction and trajectory feature extraction on data as contained in hdf5 files produced by 
CellCognition.

Input:
- plate, well: experiment of interest
- dataFolder: folder containing hdf5 files
- choice: if the goal is to extract trajectories or trajectory features
- name: generic filename for trajectory features
- repeat: if existing files should be overwritten

Note: there is a strange thing with the way iPython parses options. If you say
>>>ww="00020_01.hdf5"
>>>%run tracking/trajPack/parallel_trajFeatures -p 271114 -w $ww -d /media/lalil0u/New/projects/Xb_screen/plates__all_features_2 -c 1 --ff 0 -n features_intQC_{}.pkl
THEN it doesn't replace the first $ww with

'''
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    
    parser.add_option("-f", "--settings_file", dest="settings_file", default='tracking/settings/settings_trajFeatures.py',
                      help="Settings_file")

    parser.add_option("-p", "--plate", dest="plate",
                      help="The plate which you are interested in")
    
    parser.add_option("-w", "--well", dest="well",
                      help="The well which you are interested in")
    
    parser.add_option("-c", "--choice", dest="choice", default = False, 
                      help="False to build trajectories and true to compute features from existing trajectories")

    parser.add_option("-s", "--simulated", dest="simulated", default = 0, type=int, 
                      help="Use of simulated trajectories or no")
    parser.add_option("--cecog_file", dest="cecog_file", default = False, type=int, 
                       help="True for new type, False for old")    
    
    (options, args) = parser.parse_args()
    
    if options.well==None:
        print "You need to specify which well to treat. Pgm exiting"
        sys.exit()
    if type(options.choice)!=bool: options.choice=int(options.choice)
    
    settings = settings.Settings(options.settings_file, globals())
    outputFolder = os.path.join(settings.outputFolder, options.plate)
    fi=settings.traj_filename.format(options.well)
#     
#     if options.time_window is not None:
#         fi_trajfeatures = settings.feature_filename.format(options.well[:-5], options.time_window)
# #        time_window=time_windows[options.time_window]
#     else:
#         fi_trajfeatures = settings.feature_filename.format(options.well[:-5],'N')
# #        time_window=None
        
    if options.simulated:
        training=True
        outputFolder = os.path.join('../resultData/simulated_traj/simres/plates', options.plate)
        fi='{}--{}.pickle'.format(options.plate, options.well)
    
    try:
        os.mkdir(outputFolder)
    except OSError:
        print "Folder ", outputFolder, 'already created'
    
    if not settings.repeat:
        if not options.choice and fi in os.listdir(outputFolder):
            print "Trajectories already generated"
            sys.exit()
        elif options.choice and fi_trajfeatures in os.listdir(outputFolder):
            print "Trajectories features already calculated"
            sys.exit()
    
    print options.cecog_file
    if not options.choice: 
        print '### \n # \n ###\n We are going to predict trajectories for plate {}, well {}'.format(options.plate, options.well)
#FOR PREDICting DATA
        #i. Compute frameLots and save it
        settings.outputFolder=outputFolder
        predictor = TrackPrediction(options.plate, options.well, settings,new_cecog_files=bool(options.cecog_file))
        totalFrameLots, FEATURE_NUMBER = predictor.gettingSolu()
        del predictor
        
        if totalFrameLots is not None:
            if not os.path.isdir(os.path.join(outputFolder, 'temp')):
                os.mkdir(os.path.join(outputFolder, 'temp'))
            f=open(os.path.join(outputFolder, 'temp', 'frameLots{}{}.pkl'.format(options.plate, options.well)), 'w')
            pickle.dump(totalFrameLots, f); f.close()
            
        #ii. From framelots, computing batches of predicted solutions and saving them to be opened in the track builder
            #Need to decide how many batches
            num_batches = len(frameLots[options.plate][options.well])/100+1
            num_frames = np.max([el for el in frameLots[options.plate][options.well]])
            
            for k in range(num_batches):
                currFrameLots=frameLots()
                currFrameLots.lstFrames={options.plate:
                                         {options.well:
                                          {el:totalFrameLots[options.plate][options.well][el] for el in range(k*100, min((k+1)*100, num_frames))}}}
                if training == False:
                    singlets, doublets = currFrameLots.getAllUplets(outputFolder)
                else:
                    singlets, doublets = currFrameLots.getTrainingUplets(outputFolder)
                # print "TIME TIME TIME after getting all uplets", time.clock()
                print "joining uplets now"
                #j est dans tracking.test
                solutions = j(singlets, doublets, FEATURE_NUMBER, training)
                #print "TIME TIME TIME after joining", time.clock()
                print "normalization"
                
        
        
#ICI ON RECUPERE DONC LES SINGLETS ET DOUBLETS AVEC LA VALEUR DU TRAINING DANS CELL.TO SI ILS Y SONT, NONE SINON
    #POUR LE CENTRE ET LES FEATURES C'EST LA MOYENNE DES OBJETS DU SINGLET

# # # 
# # #     solutions = j(singlets, doublets, FEATURE_NUMBER, training)
# # #     #print "TIME TIME TIME after joining", time.clock()
# # #     print "normalization"
# # #     
# # #     fichier = open(os.path.join(loadingFolder,"minMax_data_all.pkl"), "r")  
# # #     minMax = pickle.load(fichier)
# # #     fichier.close()
# # #     try:
# # #         solutions.normalisation(minMax)
# # #     except AttributeError:
# # #         sys.stderr.write('No tracking hypotheses could be computed for this video. It is very likely that there is only one frame.')
# # #         sys.exit()
# # #     #print "TIME TIME TIME after normalization", time.clock()
# # #     
# # #     return solutions




# # #         tSol=gettingSolu(plate, well, loadingFolder, dataFolder, outputFolder, training_only, first, new_cecog_files = new_cecog_files, intensity_qc_dict=intensity_qc_dict,
# # #                      separating_function=separating_function)
# # #         first=False
# # #         new_sol = sousProcessClassify(tSol, loadingFolder)
# # #         print "Building trajectories for predicted data"
# # #         dicTraj, conn, movie_length =trackletBuilder(new_sol, outputFolder, training=False)

#         if d is not None:
#             #saving results
#             
#             f=open(os.path.join(outputFolder, fi), 'w')
#             pickle.dump(dict(zip(['tracklets dictionary', 'connexions between tracklets', 'movie_length'], [d, c, movie_length])), f)
#             f.close()
#     #    pickle.dump(dict(zip(['tracklets dictionary', 'connexions between tracklets', 'tracklets features', 'tracklets coordinates'], 
#     #                         [d, c, features, coordonnees])), f); f.close()
#         else:
#             sys.stderr.write('No output for plate {}, well {}'.format(options.plate, options.well))