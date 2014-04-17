import os, pdb, time, sys, argparse
import numpy as np
from math import fabs
from PyPack import fHacktrack2
from importPack import imp, nb_folds, EVENTS, FEATURE_NUMBER, c_list
from dataPack import treatments, joining, cplexSolving, classify
from joblib import Parallel, delayed, Memory
import cPickle as pickle
from scipy import spatial

#import PyPack

count = 0
monOutput = ""

def gettingSolu():
    global monOutput
    global FEATURE_NUMBER
    
    newFrameLot = None
    listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
    for plate in listD:
        print plate
        listW = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate)
        for well in listW:
            well=well[:-5]
            print well
            filename = '/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate+"/"+well+".hdf5"
            filenameT = '/media/lalil0u/New/workspace2/Tracking/data/trainingset/PL'+plate+"___P"+well+"___T00000.xml"
            
            
    
            monOutput+="plate = "+plate+",well = "+well+"\n"
            #ajout du frameLot et du tabF
            frameLotC, tabFC = gettingRaw(filename, filenameT, plate, well)
            if newFrameLot == None:
                newFrameLot = frameLotC 
            else: newFrameLot.addFrameLot(frameLotC)
    
    print "final training set content :"
    count, total= newFrameLot.statisticsTraining2()
    print count, total
    print "uplets now"
    #ICI ON RECUPERE DONC LES SINGLETS ET DOUBLETS AVEC LA VALEUR DU TRAINING DANS CELL.TO SI ILS Y SONT, NONE SINON
    #POUR LE CENTRE ET LES FEATURES C'EST LA MOYENNE DES OBJETS DU SINGLET
    singlets, doublets = newFrameLot.getTrainingUplets()
    
#    newFrameLot = None
#    listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/oldRaw')
#    for plate in listD:
#        print plate
#        listW = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/oldRaw/'+plate+'/hdf5')
#        for well in listW:
#            well=well[:-5]
#            print well
#            filename = '/media/lalil0u/New/workspace2/Tracking/data/oldRaw/'+plate+"/hdf5/"+well+".hdf5"
#            filenameT = '/media/lalil0u/New/workspace2/Tracking/data/trainingset/PL'+plate+"___P"+well+"___T00000.xml"
#            
#            
#    
#            monOutput+="plate = "+plate+",well = "+well+"\n"
#            #ajout du frameLot et du tabF
#            frameLotC, tabFC = gettingRaw(filename, filenameT, plate, well, True)
#            if newFrameLot == None:
#                newFrameLot = frameLotC 
#            else: newFrameLot.addFrameLot(frameLotC)
#    
#    print "final training set content :"
#    count, total= newFrameLot.statisticsTraining2()
#    print count, total
#    print "uplets now"
#    #ICI ON RECUPERE DONC LES SINGLETS ET DOUBLETS AVEC LA VALEUR DU TRAINING DANS CELL.TO SI ILS Y SONT, NONE SINON
#    #POUR LE CENTRE ET LES FEATURES C'EST LA MOYENNE DES OBJETS DU SINGLET
#    osinglets, odoublets = newFrameLot.getTrainingUplets()
    
    return singlets#, osinglets, odoublets

def gettingRaw(filename, filenameT, plate, well, old=False):
    global FEATURE_NUMBER
    print "images loading : plate = "+plate+",well = "+well    
    tabF = None
    
    frameLotC = imp.importRawSegFromHDF5(filename, plate, well, old)
    tabF=imp.importFeaturesOnly(filename, plate, well)
    if filenameT is not None:
        print "training set loading : filename ="+filenameT
        frameLotC.addTraining2(filename, filenameT)
        print "current training set content :"
        print frameLotC.statisticsTraining2()
    
    return frameLotC, tabF
    
    


if __name__ == '__main__':
    singlets = gettingSolu()
#    osinglets, singlets = gettingSolu()
#    toCheck = {}
#    for plate in singlets:
#        if plate not in toCheck:
#            toCheck[plate]={}
#        for well in singlets[plate]:
#            if well not in toCheck[plate]:
#                toCheck[plate][well]={}
#            #    compteur = 0
#            for index in singlets[plate][well]:
#                if index not in toCheck[plate][well]:
#                    toCheck[plate][well][index]=[]
#           #     z=0    
#                sing = singlets[plate][well][index]
#                try:
#                    osing = osinglets[plate][well][index]
#                except:
#                    raise
#                centers = []
#                for s in filter(lambda x: x.label !=-1, sing):
#                    centers.append(s.center)
#                if len(centers)<=1:
#                    continue
#                otree = spatial.cKDTree(centers, leafsize = 10)
#                
#                for os in filter(lambda x: x.label !=-1, osing):
#                    dou, i = otree.query(os.center, 1, distance_upper_bound =10)
#                    #pdb.set_trace()
#                    if dou<100:
#                        s = sing[i+1]
#                        roi = s.features[236]; oroi = os.features[236]
#                        if fabs(roi - oroi)>5:
#                            print "old roi ", oroi, "and roi ", roi
#                            toCheck[plate][well][index].append(os.center)
#                            print "PLATE ", plate, "WELL ", well, "INDEX ", index, "CENTER ", os.center, "TO CHECK"
#                            #z+=1
#                    else:
#                        toCheck[plate][well][index].append(os.center)
#                        #z+=1
#                        print "WOUUUUUUUUUUUUUUUUUUUUUUUUU cellule disparue dans la nouvelle version"
#                        print "PLATE ", plate, "WELL ", well, "INDEX ", index, "CENTER ", os.center, "TO CHECK"
#          #compteur = max(compteur, z)
#                        
#            f = open('PL'+plate+"___P"+well+'___T00000.xml', 'w')
#            f.write(fHacktrack2.initXml())
#            for index in toCheck[plate][well]:
#                zz=1
#                for chpt in toCheck[plate][well][index]:
#                    result = "\n      <Marker_Type>\n        <Type>{0}</Type>".format(zz)
#                    result+="\n        <Marker>\n          <MarkerX>{1}</MarkerX>\n          <MarkerY>{2}</MarkerY>\n          <MarkerZ>{0}</MarkerZ>\n        </Marker>".format(index+1, chpt[0],chpt[1])
#        
#                    result+="\n      </Marker_Type>"
#                    f.write(result)
#                    zz+=1
#            f.write (fHacktrack2.finirXml())
#            f.close()
         #   fHacktrack2.ecrireClassDef(compteur)
                        
                    
                #en fait l'ideal ce serait les noyaux douteux de les mettre dans un fichier xml comme des annotations et alors je pourrais avoir
                #en meme temps les tracking annotations et les annotations de bizarrerie    
                #mais que faire une fois ici ? car les labels ont probablement change... il faut donc pour chaque singlet essayer de trouver celui qui est le plus proche
                #en termes de centre puis voir les differences de roisize
    
    #roisize = features 236 idee = comparer les roisize pour voir si la segmentation a change
    
    
    
