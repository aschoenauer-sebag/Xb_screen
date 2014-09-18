import os, pdb, time, random
import cPickle as pickle
import numpy as np
from importPack import imp, FEATURE_NUMBER
from dataPack import treatments, joining, classify
from PyPack import fHacktrack2
from test import gettingRaw, j
import datetime
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as p

from trajPack.trajFeatures import trackletBuilder, ordonnancement
from analyzer.rough_analysis import featureEvolOverTime
#import PyPack

count = 0
monOutput = ""

def gettingSolu(loadingFolder, hdf5Folder = '/media/lalil0u/New/workspace2/Tracking/data/predict/', plate=None, wellL=None, training = False, first=True, xb_screen=False):
    global monOutput
    global FEATURE_NUMBER
    
    tabF = None
    #var pour indiquer si on prend un toy example ou pas
    #tableau qui contiendra toutes les features de tlm pr voir lesquelles contiennent des NaN
#    print "current directory ", os.getcwd()
    fichier = open(loadingFolder+"featuresToDelete.pkl", 'r')
    f = pickle.load(fichier)
    fichier.close()

    newFrameLot = None
        
    listP = os.listdir(hdf5Folder)
    if plate is not None:
        listP=[plate]
        
    for plate in listP:
        print plate
        listW = sorted(os.listdir(os.path.join(hdf5Folder, plate, 'hdf5')))
        if wellL is not None:
            listW = wellL
        for well in listW[:18]:
            well=well.split('_')[0]+'_'+well.split('_')[1][:2]
            if not xb_screen:
                filename = os.path.join(hdf5Folder, plate, 'hdf5', well+".hdf5")
            else:
                filename = os.path.join(hdf5Folder, plate, 'hdf5', well+".ch5")
            print well

                                
            if training:
                filenameT = '/media/lalil0u/New/workspace2/Tracking/data/trainingset/PL'+plate+"___P"+well+"___T00000.xml"
            else:
                filenameT = None
                
            monOutput+="plate = "+plate+",well = "+well+"\n"
            #ajout du frameLot et du tabF
            frameLotC, tabFC = gettingRaw(filename, filenameT, plate, well, name_primary_channel='primary__primary3')
            if newFrameLot == None:
                newFrameLot = frameLotC 
            else: newFrameLot.addFrameLot(frameLotC)
            tabF = tabFC if tabF == None else np.vstack((tabF, tabFC))
    
    #en ce qui concerne le nettoyage des NaN
    c, f2 = treatments.whichAreNan(tabF)
    print len(f2.keys()), len(f)
    #if there are features with NaN entries in the predict data but not in the training data
    toZeros = filter(lambda x: x not in f, f2.keys())
    if toZeros !=[]:
        print "Attention attention, some features here have NaN entries, and this was not the case in the training set. They are put to 0"
        newFrameLot.zeros(toZeros)
    
#    if f!= f2.keys():
#        featuresNames = imp.importFeaturesNames(filename)
#        print filter(lambda x: x not in f, f2.keys())
#        pdb.set_trace()
    newFrameLot.clean(f) ##np.delete(X, f, 1)
    if first:
        FEATURE_NUMBER -=len(f)

    #FEATURE_NUMBER -=len(bou)
    
    print FEATURE_NUMBER
    print "Getting all uplets now"
    print "TIME TIME TIME", time.clock()
    #ICI ON RECUPERE DONC LES SINGLETS ET DOUBLETS AVEC LA VALEUR DU TRAINING DANS CELL.TO SI ILS Y SONT, NONE SINON
    #POUR LE CENTRE ET LES FEATURES C'EST LA MOYENNE DES OBJETS DU SINGLET
    if training == False:
        singlets, doublets = newFrameLot.getAllUplets(loadingFolder)
    else:
        singlets, doublets = newFrameLot.getTrainingUplets()
    print "TIME TIME TIME after getting all uplets", time.clock()
    print "joining uplets now"

    solutions = j(singlets, doublets, FEATURE_NUMBER, training)
    print "TIME TIME TIME after joining", time.clock()
    print "normalization"
    
    fichier = open(loadingFolder+"minMax_data_all.pkl", "r")  
    minMax = pickle.load(fichier)
    fichier.close()
    solutions.normalisation(minMax)
    print "TIME TIME TIME after normalization", time.clock()
    
    return solutions

def sousProcessClassify(sol, loadingFolder, loadingFile= None, i =None, n_f =None, n_big_f=None):
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
    #pdb.set_trace()
    zz=0
    if loadingFile == None:
        for solu in sol.lstSolutions:
            new_sol.append((solu, solu.truthVec()))
            zz+=1
#            if zz>3:
#                break
    else:
        new_sol = classify.read_examples(loadingFile)

    for (x,y) in new_sol:
        r1=[]; r2=[]; truths=[]
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
            truths.append(r1)
     
    return new_sol

def trajBuilder(new_sol):
    stories = ordonnancement(new_sol); res = {}

    for plate in stories:
        res.update({plate : {}})
        for well in stories[plate]:
            print plate, well
            ensembleTraj = fHacktrack2.ensTraj()
            res[plate].update({well:ensembleTraj})
#            print ensembleTraj.numMitoses
#            print ensembleTraj.lstTraj

            num = 0; l=[0]
            beginningFrames=[0,128,129,75]
            for index in stories[plate][well]:
                print "INDEX ", index
                solC = stories[plate][well][index]
                if type(solC.truth)==str:
                    l.append(index+1)
                    continue
                singlets = solC.singlets
                nextSinglets = solC.nextSinglets
                result =filter(lambda x: x[1]==1, zip(solC.hypotheses, solC.truth))
                for hyp in result:
                    s, t = hyp[0]; source = [] ; target = []; 
                    if type(s)==np.int32 and s!=-1:
                        source.append(filter(lambda x: x.label == s, singlets)[0])
                    elif type(s)==tuple:        
                        for c in s:
                            source.append(filter(lambda x: x.label == c, singlets)[0])
                    if type(t)==np.int32 and t!=-1:
                        target.append(filter(lambda x: x.label == t, nextSinglets)[0])
                    elif type(t)==tuple:        
                        for c in t:
                            target.append(filter(lambda x: x.label == c, nextSinglets)[0])  
                    
                    incoming = len(source) ; outcoming= len(target)
                    mit = (outcoming>1) 
#                    if mit: pdb.set_trace()
                    num+=1
                    if incoming == 0:
#                        print "APPEAR"
                        for c in target:
                            newTraj= fHacktrack2.trajectoire(num, c.center[0], c.center[1], index+1, c.label)
                            newTraj.longueurs=[index+1]
                            newTraj.fractions=[1]
                            ensembleTraj.add(newTraj)
                            num+=1
                        continue
                    elif outcoming ==0:
#                        print "DISAPPEAR"
                        if index in l:
                            #pdb.set_trace()#en fait ce sont au moins les trajectoires qui commencent en 0 et finissent en 0 ie les apparitions disparitions en 0
                            for c in source:
                                newTraj = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label)
                                newTraj.longueurs = [index]
                                newTraj.fractions=[1]
                                ensembleTraj.add(newTraj)
                                num+=1
                        continue
                    elif incoming == 1: 
#                        print "FROM ONE"
                        c = source[0]
                        incomingTraj = ensembleTraj.findLabelFrame(c.label, index)
                        if outcoming <= len(incomingTraj):
                            first = True; ll=[]; lf=[]
#                            for traj in incomingTraj:
#                                ll.extend(traj.longueurs)
#                                lf.extend(traj.fractions)
#                            if ll==[] and index>0:
#                                pdb.set_trace()
                            for traj in incomingTraj:
                                try:
                                    t = target.pop()
                                except IndexError:
                                    continue
                                else:
                                    traj.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
#                                    traj.longueurs = ll; traj.fractions = lf
                                    if first: numC = traj.numCellule
                                    if mit: 
#                                        print "ajout mitose"
                                        traj.mitose=1; traj.numCellule = numC; first = False
                                        #traj.ajoutMomentInteressant(c.center[0], c.center[1], index)
                                        ensembleTraj.numMitoses+=1
                        else:
                            i = 0; first = True; ll=[]; lf=[]
#                            for traj in incomingTraj:
#                                ll.extend(traj.longueurs)
#                                lf.extend(traj.fractions)
#                            if ll==[] and index>0:
#                                pdb.set_trace()    
                            for t in target:
                                try:
                                    incomingTraj[i].ajoutPoint(t.center[0], t.center[1], index+1, t.label)
#                                    incomingTraj[i].longueurs = ll; incomingTraj[i].fractions = lf
                                    if first: numC = incomingTraj[i].numCellule                                   
                                    if mit: 
#                                        print "ajout mitose"
                                        incomingTraj[i].mitose=1; incomingTraj[i].numCellule = numC; first = False
                                        #incomingTraj[i].ajoutMomentInteressant(c.center[0], c.center[1], index)
                                        ensembleTraj.numMitoses+=1
                                    i+=1
                                except IndexError:
                                    if not mit: 
                                        #ici on n'y arrive a t=0 quand aucune trajectoire n'est initialise
                                        if index not in beginningFrames: pdb.set_trace()
                                        num+=1
                                        newTraj = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label)
                                        newTraj.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
                                        newTraj.longueurs=[index]; newTraj.fractions = [1]
                                    if first: 
                                        numC = num    
                                    if mit:
                                        newTraj = fHacktrack2.trajectoire(numC, c.center[0], c.center[1], index, c.label)
                                        newTraj.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
                                        newTraj.longueurs = incomingTraj[0].longueurs if incomingTraj !=[] else [index]
                                        newTraj.fractions = incomingTraj[0].fractions if incomingTraj !=[] else [1]
                                        #pdb.set_trace()
#                                        print "ajout mitose"
                                        newTraj.mitose=1
                                        first = False
                                        #newTraj.ajoutMomentInteressant(c.center[0], c.center[1], index)
                                        ensembleTraj.numMitoses+=1
                                    #pdb.set_trace()
                                    ensembleTraj.add(newTraj)
            
                    elif incoming >1:
                        if outcoming>1:
                            raise
#                        print "FROM MORE THAN ONE"
                        t = target[0]; tt=None
                        if index not in beginningFrames:
                            tt=[]; ll=[]; lf=[]
                            for c in source:
                                incomingTraj = ensembleTraj.findLabelFrame(c.label, index)
                                tt.append(incomingTraj[0])
                                k=0
                                for lon in incomingTraj[0].longueurs:
                                    ll.append(lon)
                                    lf.append(incomingTraj[0].fractions[k]/float(incoming))
                                    k+=1
#                        if ll==[] and index>0:
#                            pdb.set_trace()        
                        i=0
                        for c in source:
                            try:
                                trajC = tt[i]
                                trajC.longueurs = ll; trajC.fractions = lf
#                                pdb.set_trace()
                            except:
                                #on tombe ici a t=0 mais est-ce qu'on y tombe apres ?
                                if index not in beginningFrames: pdb.set_trace()
                                num+=1
                                trajC = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label)
                                trajC.longueurs=[index for bb in range(incoming)]; trajC.fractions = [1/float(incoming) for bb in range(incoming)]
                                ensembleTraj.add(trajC)
                            finally:    
                                trajC.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
                            i+=1
                            #incomingTraj[0].ajoutMomentInteressant(c.center[0], c.center[1], index)
    return res

def outputTrajXml(dicTraj, savingFolder):
    compteurMax = 0; first=True
    for plate in dicTraj:
        for well in dicTraj[plate]:
            ensTraj = dicTraj[plate][well]
            nomFichier = "PL"+plate+"___P"+well+"___T00001.xml"
            nomFichier = savingFolder+"/"+nomFichier

            coord = fHacktrack2.sortiesAll(nomFichier, ensTraj)
            #plot3d(coord, plate, well, 10000)
            #c=fHacktrack2.sortiesAllBis(nomFichier, ensTraj, compteurMax)
            #compteurMax = max(compteurMax, c)
    return

if __name__ == '__main__':

    now = datetime.datetime.now()
    
    savingFolder = "/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions"
    hdf5Folder = "/media/lalil0u/New/projects/Xb_screen/plates_new_seg"
    plate = '130814'
    
    well_groups=[['{:>05}_01'.format(k) for k in [1,2,3,5,7,8,9,10,11,12]], #endosulfan
                 ['{:>05}_01'.format(k) for k in [4,15,26]],                  #tgf beta
                 ['{:>05}_01'.format(k) for k in [13,14,16,17,19,20,21,22,23]], #tcdd
                 ['{:>05}_01'.format(k) for k in [6,29,48,52,55]], #dmso 
                 ['{:>05}_01'.format(k) for k in [18,24,41,60,61]], #nonane
                 ['{:>05}_01'.format(k) for k in [25,27,28,30,31,32,33,34,35,36]], #bpa
                 ['{:>05}_01'.format(k) for k in [37,38,39,40,42,43,44,45,46,47]], #pcb
                 ['{:>05}_01'.format(k) for k in [49,50,51,53, 54,56, 57,58,59]], #mehg
                 ]
    well_groups=zip(well_groups, ['alpha-Endosulfan', 'TGF beta 1', 'TCDD', 'DMSO', 'Nonane', 'Bisphenol A', 'PCB 153', 'MeHg'])
    if not os.path.isdir(savingFolder):
        os.mkdir(savingFolder)
    loadingFolder = "../prediction/"
    
    print "TIME TIME TIME", time.clock()
#    try:
    first=True
    for wellL in well_groups[5:]:
        sol = gettingSolu(loadingFolder, hdf5Folder, plate,wellL[0], xb_screen=True, first=first)
        first=False
        
        print "TIME TIME TIME", time.clock()
        print "classifying data from :", loadingFolder
        
        new_sol = sousProcessClassify(sol, loadingFolder)
        
        print "TIME TIME TIME", time.clock()
        print "traj building"
        print "Building trajectories for predicted data"
        if True:
            dicTraj, conn, movie_length =trackletBuilder(new_sol, loadingFolder, training=True)
            featureEvolOverTime(dicTraj, conn, savingFolder, verbose=True, movie_length=movie_length, name = wellL[-1])
        else:
            try:    
                dicTraj =trajBuilder(new_sol)
            except:
                pass
            else:
                f = open(os.path.join(savingFolder, "dicTraj_{}_{}.pkl".format(plate, wellL[-1])), "w")
                pickle.dump(dicTraj, f); f.close()
        #    f = open(loadingFolder+"dicTraj_41_35.pkl", "r")
        #    dicTraj = pickle.load(f); f.close()
        #    
            print "TIME TIME TIME", time.clock()
            print "saving as xml in", loadingFolder
            outputTrajXml(dicTraj, savingFolder)
            print "TIME TIME TIME", time.clock()
