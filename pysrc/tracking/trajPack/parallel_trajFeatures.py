#IMPORT
import pdb, os, time, sys
import numpy as np
import cPickle as pickle
from optparse import OptionParser
import matplotlib
from tracking.PyPack.fHacktrack2 import trajectoire, ensTraj
matplotlib.use('Agg')

from tracking.importPack import FEATURE_NUMBER
from tracking.test import gettingRaw, j
from tracking.dataPack import treatments
from tracking.trajPack import densities
from trajFeatures import trackletBuilder, histogramPreparationFromTracklets
from tracking.trackingF import sousProcessClassify

'''
    To perform trajectory extraction and trajectory feature extraction on a cluster node.
    Relevant information should be passed on the command line, cf infra
'''

    
def gettingSolu(plate, w, loadingFolder, dataFolder, outputFolder, training = False, first=True, new_cecog_files=False):
    global FEATURE_NUMBER
    
    tabF = None
    #tableau qui contiendra toutes les features de tlm pr voir lesquelles contiennent des NaN

    fichier = open(os.path.join(loadingFolder,"featuresToDelete.pkl"), 'r')
    f = pickle.load(fichier)
    fichier.close()

    newFrameLot = None
    listP = os.listdir(dataFolder)
    if w is not None:
        listP = [w] 
    for well in listP:
        well=well[:-5]
        #print well
        if not new_cecog_files:
            filename = os.path.join(dataFolder, well+".hdf5")
        else:
            filename = os.path.join(dataFolder, well+".ch5")
        if training:
            filenameT = '/media/lalil0u/New/workspace2/Tracking/data/trainingset/PL'+plate+"___P"+well+"___T00000.xml"
        else:
            filenameT = None
            #ajout du frameLot et du tabF
        frameLotC, tabFC = gettingRaw(filename, filenameT, plate, well, name_primary_channel='primary__primary')
        if newFrameLot == None:
            newFrameLot = frameLotC 
        else: newFrameLot.addFrameLot(frameLotC)
        tabF = tabFC if tabF == None else np.vstack((tabF, tabFC))
    
    #en ce qui concerne le nettoyage des NaN
    c, f2 = treatments.whichAreNan(tabF)
    #print len(f2.keys()), len(f)
    #if there are features with NaN entries in the predict data but not in the training data
    toZeros = filter(lambda x: x not in f, f2.keys())
    if toZeros !=[]:
        sys.stderr.write("Attention attention, plate {}, some features here have NaN entries, and this was not the case in the training set. They are put to 0".format(plate))
        newFrameLot.zeros(toZeros)

    newFrameLot.clean(f) ##np.delete(X, f, 1)
    if first:
        FEATURE_NUMBER -=len(f)
   
    #ICI ON RECUPERE DONC LES SINGLETS ET DOUBLETS AVEC LA VALEUR DU TRAINING DANS CELL.TO SI ILS Y SONT, NONE SINON
    #POUR LE CENTRE ET LES FEATURES C'EST LA MOYENNE DES OBJETS DU SINGLET
    if training == False:
        singlets, doublets = newFrameLot.getAllUplets(outputFolder)
    else:
        singlets, doublets = newFrameLot.getTrainingUplets(outputFolder)
   # print "TIME TIME TIME after getting all uplets", time.clock()
    print "joining uplets now"

    solutions = j(singlets, doublets, FEATURE_NUMBER, training)
    #print "TIME TIME TIME after joining", time.clock()
    print "normalization"
    
    fichier = open(os.path.join(loadingFolder,"minMax_data_all.pkl"), "r")  
    minMax = pickle.load(fichier)
    fichier.close()
    try:
        solutions.normalisation(minMax)
    except AttributeError:
        sys.stderr.write('No tracking hypotheses could be computed for this video. It is very likely that there is only one frame.')
        sys.exit()
    #print "TIME TIME TIME after normalization", time.clock()
    
    return solutions

def output(plate,  well, allDataFolder, outputFolder, training_only=True):
    #listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
    first = True; 
    dataFolder = os.path.join(allDataFolder, plate, 'hdf5')##################################### 
    loadingFolder = "/cbio/donnees/aschoenauer/workspace2/Xb_screen/prediction"########################################################################################"
    tSol=gettingSolu(plate, well, loadingFolder, dataFolder, outputFolder, training_only, first)
    first=False
    new_sol = sousProcessClassify(tSol, loadingFolder)
    print "Building trajectories for predicted data"
    dicTraj, conn, movie_length =trackletBuilder(new_sol, outputFolder, training=False)
##        del new_sol;del tSol;
    return dicTraj, conn, movie_length

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
'''
    initTime=time.clock()
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    parser.add_option("-p", "--plate", dest="plate",
                      help="The plate which you are interested in")
    
    parser.add_option("-w", "--well", dest="well",
                      help="The well which you are interested in")

    parser.add_option("-d", "--data_folder", dest="dataFolder", default='',
                      help="Give absolute path to data")
    
    parser.add_option("-c", "--choice", dest="choice", default = False, 
                      help="False to build trajectories and true to compute features from existing trajectories")

    parser.add_option("-n", "--name", dest="filename", default = 'hist_tabFeatures_{}.pkl', 
                      help="Filename for trajectory features")
    
    parser.add_option("-s", "--simulated", dest="simulated", default = 0, type=int, 
                      help="Use of simulated trajectories or no")
    
    parser.add_option("-r", "--repeat", dest="repeat", default = False, 
                      help="False to do only videos that haven't been treated yet and true to compute features even if already computed")

    (options, args) = parser.parse_args()
    
    if options.well==None:
        print "You need to specify which well to treat. Pgm exiting"
        sys.exit()
    
    training=False
    predict = True
    TroisD = False
    if type(options.choice)!=bool: options.choice=int(options.choice)
    outputFolder = os.path.join("/share/data20T/mitocheck/tracking_results", options.plate)#######################################################################################
    fi = 'traj_noF_densities_w{}.pkl'.format(options.well)
    
    fi_trajfeatures = options.filename.format(options.well[:-5])
    
    if options.simulated:
        training=True
        outputFolder = os.path.join('../resultData/simulated_traj/simres/plates', options.plate)
        fi='{}--{}.pickle'.format(options.plate, options.well)
    
    try:
        os.mkdir(outputFolder)
    except OSError:
        print "Folder ", outputFolder, 'already created'
    
    if not options.repeat:
        if not options.choice and fi in os.listdir(outputFolder):
            print "Trajectories already generated"
            sys.exit()
        elif options.choice and fi_trajfeatures in os.listdir(outputFolder):
            print "Trajectories features already calculated"
            sys.exit()
        
    if not options.choice: 
        print '### \n # \n ###\n We are going to predict trajectories for plate {}, well {}'.format(options.plate, options.well)
        print 'Densities distance ', densities
#FOR PREDICTED DATA
        d, c, movie_length=output(options.plate,options.well, options.dataFolder,outputFolder, training); 

        if d is not None:
            #saving results
            
            f=open(os.path.join(outputFolder, fi), 'w')
            pickle.dump(dict(zip(['tracklets dictionary', 'connexions between tracklets', 'movie_length'], [d, c, movie_length])), f)
            f.close()
    #    pickle.dump(dict(zip(['tracklets dictionary', 'connexions between tracklets', 'tracklets features', 'tracklets coordinates'], 
    #                         [d, c, features, coordonnees])), f); f.close()
        else:
            sys.stderr.write('No output for plate {}, well {}'.format(options.plate, options.well))
    else:
        print "### \n # \n ###\n  We are going to compute features from existing trajectories for plate {}\n Adding density information".format(options.plate)
        try:
            filename = os.path.join(outputFolder, fi); print filename
            f=open(filename, 'r')
            dataDict = pickle.load(f)
            f.close()
        except:
            sys.stderr.write('Folder {} does not contain densities trajectories file.'.format(outputFolder))
            sys.exit()
                    
        else:
            if not options.simulated:
                d,c, movie_length = dataDict['tracklets dictionary'], dataDict['connexions between tracklets'], dataDict['movie_length']
                res = histogramPreparationFromTracklets(d, c, outputFolder, False, verbose, movie_length, name=options.filename)  #(d,c, outputFolder, False, verbose, tab=True, length=movie_length)
            else:
                d=ensTraj()
                for traj in dataDict:
                    t = trajectoire(1, xi=None, yi=None, frame=None, idC=None, id=1)
                    t.lstPoints = traj
                    d.lstTraj.append(t)
                res=histogramPreparationFromTracklets({options.plate : {options.well : d}}, None, 
                                                      outputFolder,training =True, verbose=verbose, movie_length={options.plate : {options.well :99}}, name=options.filename)  #(d,c, outputFolder, False, verbose, tab=True, length=movie_length)
    
    final = time.clock() - initTime
    print "##################TEMPS FINAL {}".format(final)
    
    
    
    
