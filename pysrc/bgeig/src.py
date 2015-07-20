import os
import numpy as np
import cPickle as pickle

from tracking.PyPack.fHacktrack2 import initXml, ecrireXml, finirXml
from util.settings import Settings

def imageRenaming(folder):
    currFolder=os.getcwd()
    os.chdir(folder)
    plates=os.listdir(os.getcwd())
    for plate in plates:                               
        wells=os.listdir(plate)    
        for well in wells:                                                      
            for image in filter(lambda x: 'TIF' in x, os.listdir(os.path.join(plate, well))):
                frame = int(image.split('.')[1])
                if 'phase' in image:
                    os.rename(os.path.join(plate, well, image), os.path.join(plate, well, '--W{:>04}--P00001_t{:>05}_c1.tif'.format(int(well[1:]), frame)))
                elif 'calibrate' in image:
                    os.rename(os.path.join(plate, well, image), os.path.join(plate, well, '--W{:>04}--P00001_t{:>05}_c0.tif'.format(int(well[1:]), frame)))
                    
    for plate in plates:                               
        wells=os.listdir(plate)    
        for well in wells:                                                               
            os.rename(os.path.join(plate, well), os.path.join(plate, 'W{:>04}'.format(int(well[1:]))))
    os.chdir(currFolder)
    return

#then for predicting trajectories, writing xml files and computing features:
# for w in ['{:>05}_01.ch5'.format(k) for k in range(21, 25)]:
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic3 -w $w -c 0 --cecog_file 1   
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic2 -w $w -c 0 --cecog_file 1
# for w in ['{:>05}_01.ch5'.format(k) for k in range(21, 25)]:
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic3 -w $w -c 1 --cecog_file 1   
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic2 -w $w -c 1 --cecog_file 1
# for w in ['{:>05}_01.ch5'.format(k) for k in range(21, 25)]:
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic3 -w $w -c 2 --cecog_file 1   
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic2 -w $w -c 2 --cecog_file 1

def findingSharpMovementDistribution(setting_file='bgeig/settings/settings_bgeig.py', measure='turningangle'):
    settings=Settings(setting_file, globals())
    
    plates=os.listdir(settings.dataFolder)
    result=[]
    for plate in plates:
        wells=['{:>05}_01'.format(k) for k in range(21, 25)]
        
        for well in wells:
            f=open(os.path.join(settings.outputFolder, plate, settings.feature_filename.format(well)))
            arr, coord, hist=pickle.load(f); f.close()
        
            for el in hist[measure]:
                result.extend(el)
                
    return result

def findingSharpMovements(setting_file='bgeig/settings/settings_bgeig.py', measure='turningangle'):
    settings=Settings(setting_file, globals())
    
    plates=os.listdir(settings.dataFolder)
    for plate in plates:
        wells=['{:>05}_01'.format(k) for k in range(21, 25)]
        
        for well in wells:
            f=open(os.path.join(settings.outputFolder, plate, settings.feature_filename.format(well)))
            arr, coord, hist=pickle.load(f); f.close()
            
            nomFichier = "PL"+plate+"___P"+well+"___T00001.xml"
            nomFichier = os.path.join(settings.outputFolder,'sharp_mov', 'annotations', nomFichier)
            compteur =1
            fichierX = open(nomFichier, "w")
            fichierX.write(initXml())
            
            for i,el in enumerate(hist[measure]):
                if np.any(np.array(el)>0.70):
                    print len(el), len(coord[i][0])
                    wh_=np.where(np.array(el)>0.70)[0]
                    currcoord=np.array(zip(*coord[i]))
                    d={(x[0],0):(x[1], x[2]) for x in currcoord[wh_]}
                    
                    txt, _ = ecrireXml(compteur,d, True)
                    fichierX.write(txt)
                    compteur+=1
            
            fichierX.write(finirXml())
            fichierX.close()
                
    return