import os, shutil, gc
import cPickle as pickle
from optparse import OptionParser
import sys, pdb
import numpy as np

from tracking.importPack.imp import importRawSegFromHDF5, frameLots
from tracking.dataPack import classify, joining
from tracking.PyPack.fHacktrack2 import initXml, ecrireXml, finirXml, sortiesAll
from tracking.trajPack.trajFeatures import trackletBuilder, histogramPreparationFromTracklets
from util.settings import Settings
import psutil

class TrackPrediction(object):
    def __init__(self, plate, w, settings, new_cecog_files=False):
        self.plate=plate
        self.well = w
        
        self.settings=settings
        self.new_cecog_files = new_cecog_files
        
        self.dataFolder=os.path.join(self.settings.dataFolder, plate, 'hdf5')
        
    def _getIntensityQC(self):
        intensity_qc_dict=None
        if self.settings.intensity_qc_file is not None:
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
    
    @staticmethod    
    def frameJoin(singlets, doublets, featSize, training= True):
        solutions= None
        for plate in singlets:
                print plate
                for well in singlets[plate]:
                    print well
                    sys.stderr.write("\n plate {}, well {}\n".format(plate, well))
                    for index in singlets[plate][well]:
                        print '-- ',
                        if index+1 not in singlets[plate][well] or singlets[plate][well][index+1]==[]:
                            continue
                        singletsL = singlets[plate][well][index]
                        nextSinglets = singlets[plate][well][index+1]
                        doubletsL = doublets[plate][well][index]
                        nextDoublets = doublets[plate][well][index+1]
                        if len(nextSinglets)==1:
                            continue
                        solution = joining.Solution(plate, well, index, singletsL, nextSinglets, doubletsL, nextDoublets, featSize, training)
                        if solutions == None:
                            solutions= joining.Solutions(solution, lstSolutions = None)
                        else:
                            solutions.append(solution)    
        return solutions
    
    @staticmethod
    def predict(sol, loadingFolder, loadingFile= None, i =None, n_f =None, n_big_f=None):
        #subprocess.call(["/media/lalil0u/New/software2/downloads/unpacked/svm-python-v204/svm_python_classify", "--m", "test", "-v", "3", "results/data_TEST_fold"+str(n_fold)+".pkl", "results/modelfile_"+c+"_"+str(n_fold)+".pkl"])
        if i==None:
            f = open(os.path.join(loadingFolder,"modelfile_all.pkl"), 'r')
        elif n_f==None:
            c = "{:f}".format(10**i)
            if n_big_f==None:
                f = open(os.path.join(loadingFolder,"modelfile_all"+c+".pkl"), 'r')
            else:
                f = open(os.path.join(loadingFolder,"modelfile_all"+c+"_"+str(n_big_f)+".pkl"), 'r')
        else:
            c = "{:f}".format(10**i)
            if n_big_f is not None:
                f = open(os.path.join(loadingFolder,"modelfile_"+c+"_"+str(n_big_f)+'_'+str(n_f)+".pkl"), 'r')
            else:
                f = open(os.path.join(loadingFolder,"modelfile_"+c+"_"+str(n_f)+".pkl"), 'r')
            
        mesPoids = pickle.load(f)
        f.close()
        new_sol = []
        zz=0
        if loadingFile == None:
            for solu in sol.lstSolutions:
                new_sol.append((solu, solu.truthVec()))
                zz+=1
        else:
            new_sol = classify.read_examples(loadingFile)
    
        for x,_ in new_sol:
            r1=[]
            try:
                ybar = classify.classify_example(x, mesPoids)
            except:
                print "pbl de NaN"
                x.truth = "PROBLEME DE NaN"
                continue
            else:
                for k in range(len(ybar)):
                    r1.extend(ybar[k]) 
                x.truth = r1
         
        return new_sol
    
    def trackletBuilderPrep(self,num_batches, from_batches=True):
        if self.settings.training or not from_batches:
            raise ValueError('You are not in the right place to do what you want to do')
        totalSolutions=[]

        for k in range(num_batches):
            f=open(os.path.join(self.settings.outputFolder, 'temp', 'solutions{}{}{}.pkl'.format(self.plate, self.well,k)))
            totalSolutions.extend(pickle.load(f)); f.close()
            
        return trackletBuilder(totalSolutions, self.settings.outputFolder, training=self.settings.training)
    
    def writeTracks(self):
        '''
        To produce annotated tracklets, either all tracklets (option all_=True) or only those used for feature extraction (option all=False)
        '''
        intensity_qc_dict = self._getIntensityQC()
        
        if 'tracking_annotations' not in os.listdir(self.settings.outputFolder):
            os.mkdir(os.path.join(self.settings.outputFolder, 'tracking_annotations'))
        if 'annotations' not in os.listdir(os.path.join(self.settings.outputFolder, 'tracking_annotations')):
            os.mkdir(os.path.join(os.path.join(self.settings.outputFolder, 'tracking_annotations', 'annotations')))
        
        if self.settings.all_trajectory_xml:
            well_files = filter(lambda x: self.settings.traj_filename.split('{}')[0] in x, os.listdir(os.path.join(self.settings.outputFolder, self.plate)))
        else:
            well_files = filter(lambda x: self.settings.feature_filename.split('{}')[0] in x, os.listdir(os.path.join(self.settings.outputFolder, self.plate)))
        for file_ in well_files:
            well_num=int(file_.split('_')[2][1:])
            well = '{:>05}_01'.format(well_num)
            
            nomFichier = "PL"+self.plate+"___P"+well+"___T00001.xml"
            nomFichier = os.path.join(self.settings.outputFolder,'tracking_annotations', 'annotations', nomFichier)
    
            if self.settings.all_trajectory_xml:
                f=open(os.path.join(self.settings.outputFolder, self.plate, file_), 'r')
                d=pickle.load(f); f.close()
                ensTraj = d['tracklets dictionary'][self.plate][well]
                _ = sortiesAll(nomFichier, ensTraj)
            else:
                f=open(os.path.join(self.settings.outputFolder, self.plate, file_), 'r')
                _,coord,_=pickle.load(f); f.close()
                
                compteur =1
                fichierX = open(nomFichier, "w")
                fichierX.write(initXml())
                if intensity_qc_dict is not None and well_num in intensity_qc_dict:
                    frames_to_skip = intensity_qc_dict[well_num]
                else:
                    frames_to_skip=None
    
                for lstPoints in coord:
                    frameL=lstPoints[0]; X=lstPoints[1]; Y=lstPoints[2]
                    d={(frameL[k]+len(np.where(frames_to_skip<=frameL[k])[0]),0):(X[k], Y[k]) for k in range(len(frameL))}
                    txt, coord = ecrireXml(compteur,d, True)
                    fichierX.write(txt)
                    compteur+=1
                    if compteur>1600:
                        break
            
                fichierX.write(finirXml())
                fichierX.close()
                print 'finally ', compteur, 'trajectories'
            
        return

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
    
    parser.add_option("-f", "--settings_file", dest="settings_file", default='bgeig/settings/settings_trajFeatures.py',
                      help="Settings_file")

    parser.add_option("-p", "--plate", dest="plate",
                      help="The plate which you are interested in")
    
    parser.add_option("-w", "--well", dest="well",
                      help="The well which you are interested in")
    
    parser.add_option("-c", "--choice", dest="choice", default = False, type=int,
                      help="0 to build trajectories, 1 to compute xml files from trajectories, 2 to compute features from existing trajectories")

    parser.add_option("-s", "--simulated", dest="simulated", default = 0, type=int, 
                      help="Use of simulated trajectories or no")
    parser.add_option("--cecog_file", dest="cecog_file", default = False, type=int, 
                       help="True for new type, False for old")    
    
    (options, args) = parser.parse_args()
    
    if options.well==None:
        print "You need to specify which well to treat. Pgm exiting"
        sys.exit()
    if type(options.choice)!=bool: options.choice=int(options.choice)
    
    settings = Settings(options.settings_file, globals())
    if not os.path.isdir(settings.outputFolder):
        os.mkdir(settings.outputFolder)
        
    outputFolder = os.path.join(settings.outputFolder, options.plate)
    if not os.path.isdir(outputFolder):
        os.mkdir(outputFolder)
        
    fi=settings.traj_filename.format(options.well)
    fi_trajfeatures = settings.feature_filename.format(options.well.split('.')[0])
                                                       
    if options.simulated:
        settings.training=True
        outputFolder = os.path.join('../resultData/simulated_traj/simres/plates', options.plate)
        fi='{}--{}.pickle'.format(options.plate, options.well)
    
    if not settings.repeat:
        if options.choice==0 and fi in os.listdir(outputFolder):
            print "Trajectories already generated"
            sys.exit()
        elif options.choice==2 and fi_trajfeatures in os.listdir(outputFolder):
            print "Trajectories features already calculated"
            sys.exit()
    
    print options.cecog_file
    if options.choice==0: 
        print '### \n# \n###\n We are going to predict trajectories for plate {}, well {}'.format(options.plate, options.well)
#FOR PREDICting DATA
        #i. Compute frameLots and save it
        settings.outputFolder=outputFolder
        predictor = TrackPrediction(options.plate, options.well, settings,new_cecog_files=bool(options.cecog_file))
        totalFrameLots, FEATURE_NUMBER = predictor.gettingSolu()
        print psutil.virtual_memory()
        predictor=None; gc.collect()
        print psutil.virtual_memory()
         
        well=options.well.split('.')[0]
         
        if totalFrameLots is None:
            sys.stderr.write('Pbl computing framelots for {}, {}'.format(options.plate, well))
            sys.exit()
             
        if not os.path.isdir(os.path.join(outputFolder, 'temp')):
            os.mkdir(os.path.join(outputFolder, 'temp'))
        f=open(os.path.join(outputFolder, 'temp', 'frameLots{}{}.pkl'.format(options.plate, well)), 'w')
        pickle.dump(totalFrameLots, f); f.close()
         
    #ii. From framelots, computing batches of predicted solutions and saving them to be opened in the track builder
        #Need to decide how many batches
        num_batches = len(totalFrameLots.lstFrames[options.plate][well])/100+1
        num_frames = np.max([el for el in totalFrameLots.lstFrames[options.plate][well]])
         
        #loading normalization
        fichier = open(os.path.join(settings.loadingFolder,"minMax_data_all.pkl"), "r")  
        minMax = pickle.load(fichier)
        fichier.close()
         
#        print psutil.virtual_memory()
        
        for k in range(num_batches):
            if not settings.redo and 'solutions{}{}{}.pkl'.format(options.plate, well,k) in os.listdir(os.path.join(outputFolder, 'temp')):
                continue
            
            currFrameLots=frameLots()
            currFrameLots.lstFrames={options.plate:
                                     {well:#(k+1)*100+1 is important here, it is the link between the consecutive pair frames. num_frames+1 is important as well to take the last frame
                                            #into account
                                      {el:totalFrameLots.lstFrames[options.plate][well][el] for el in range(k*100, min((k+1)*100+1, num_frames+1))}}}
             
#ICI ON RECUPERE DONC LES SINGLETS ET DOUBLETS AVEC LA VALEUR DU TRAINING DANS CELL.TO SI ILS Y SONT, NONE SINON
#POUR LE CENTRE ET LES FEATURES C'EST LA MOYENNE DES OBJETS DU SINGLET
            print "Going from {} to {} in batch mode".format(k*100, min((k+1)*100+1, num_frames))
            if settings.training == False:
                singlets, doublets = currFrameLots.getAllUplets(outputFolder)
            else:
                singlets, doublets = currFrameLots.getTrainingUplets(outputFolder)
            # print "TIME TIME TIME after getting all uplets", time.clock()
            print "Joining uplets now"
            #print psutil.virtual_memory()
            solutions = TrackPrediction.frameJoin(singlets, doublets, FEATURE_NUMBER, settings.training)
            #print psutil.virtual_memory()
            singlets=None ; doublets=None ; gc.collect()
            #print psutil.virtual_memory()
            print "Feature normalization"
            try:
                solutions.normalisation(minMax)
            except AttributeError:
                sys.stderr.write('No tracking hypotheses could be computed for this video. It is very likely that there is only one frame.')
                continue
            else:
                new_sol = TrackPrediction.predict(solutions, settings.loadingFolder)
                solutions=None; gc.collect()
                print "After prediction ", psutil.virtual_memory()
                
                
                #Cleaning the solution before saving it because otherwise it's very heavy, and there's a lot of info used for prediction that
                #we're not going to need again when building tracks
                
                for solution, _ in new_sol:
                    solution.result =filter(lambda x: x[1]==1, zip(solution.hypotheses, solution.truth))
                    
                    solution.constraints=None
                    solution.events=None
                    solution.features=None
                    solution.hypotheses=None
                    
                f=open(os.path.join(outputFolder, 'temp', 'solutions{}{}{}.pkl'.format(options.plate, well,k)), 'w')
                pickle.dump(new_sol, f); f.close()
                new_sol=None; gc.collect()
                #print psutil.virtual_memory()
#
        print "Building trajectories for predicted data"
        predictor = TrackPrediction(options.plate, well, settings,new_cecog_files=bool(options.cecog_file))
        dicTraj, conn, movie_length =predictor.trackletBuilderPrep(num_batches=num_batches, from_batches=True)                

        if dicTraj is not None:
            #saving results
            f=open(os.path.join(outputFolder, fi), 'w')
            pickle.dump(dict(zip(['tracklets dictionary', 'connexions between tracklets', 'movie_length'], [dicTraj, conn, movie_length])), f)
            f.close()

        else:
            sys.stderr.write('No output for plate {}, well {}'.format(options.plate, well))
            
        if settings.removeTempFiles:
            shutil.rmtree(os.path.join(outputFolder, 'temp'))
            
            
    elif options.choice==1:
        print "We are going to write xml files from predicted trajectories"
        predictor = TrackPrediction(options.plate, None, settings,new_cecog_files=bool(options.cecog_file))
        predictor.writeTracks()
        
    elif options.choice==2:
        print "### \n # \n ###\n  We are going to compute features from existing trajectories for plate {}\n Adding density information".format(options.plate)
        try:
            filename = os.path.join(outputFolder, fi)
            f=open(filename, 'r')
            dataDict = pickle.load(f)
            f.close()
        except:
            sys.stderr.write('Folder {} does not contain densities trajectories file.'.format(outputFolder))
            sys.exit()
                    
        else:
            if not options.simulated:
                d,c, movie_length = dataDict['tracklets dictionary'], dataDict['connexions between tracklets'], dataDict['movie_length']
                res = histogramPreparationFromTracklets(d, c, outputFolder, False, verbose, movie_length, name=fi_trajfeatures,
                                                        filtering_fusion=settings.filtering_fusion) 
            else:
                raise AttributeError
#                 d=ensTraj()
#                 for traj in dataDict:
#                     t = trajectoire(1, xi=None, yi=None, frame=None, idC=None, id=1)
#                     t.lstPoints = traj
#                     d.lstTraj.append(t)
#                 res=histogramPreparationFromTracklets({options.plate : {options.well : d}}, None, 
#                                                       outputFolder,training =True, verbose=verbose, movie_length={options.plate : {options.well :99}}, 
#                                                       name=fi_trajfeatures)  #(d,c, outputFolder, False, verbose, tab=True, length=movie_length)
        
        
            