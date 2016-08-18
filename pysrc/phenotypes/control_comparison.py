import csv, os
import numpy as np
import cPickle as pickle

from analyzer import interface
from tracking.importPack.imp import importTargetedFromHDF5
from _collections import defaultdict


raw_result_dir_Mitocheck= "/share/data40T/Thomas/mitocheck_full_hdf5/out_data"
primary_channel_name = 'primary__primary3'
focus_classname = 'OutOfFocus'
artefact_classname = 'Artefact'


def extractControlDataFromMitocheck(n=200):
    #i. get n control wells from Mitocheck
    f=open('../data/mapping_2014/qc_export.txt', 'r')
    reader = csv.reader(f, delimiter='\t')
    ctrls = []
    plateList = np.array(os.listdir(raw_result_dir_Mitocheck))
    for el in reader:
        if 'Valid' not in el[0] and int(el[0][2:6])<600 and el[2]=='scrambled' and el[-1]=="ok":
            plate = el[0].split('--')[0]
            well = el[0].split('--')[1]
            truePlate = plateList[np.where([plate in el for el in plateList])[0]][0]
            
            ctrls.append((truePlate, '00{}_01'.format(well)))
    np.random.shuffle(ctrls); chosenCtrls = ctrls[:n]
    
    #ii. load the data under the form {well: data}
    result = loadData(chosenCtrls)

    #iii. save it
    f=open('~/projects/drug_screen/results/plates/mitocheck_ctrls.pkl', 'w')
    pickle.dump(result, f)
    f.close()
    
    return result
    
def loadData(ctrls):
    newFrameLot = None
    dataFolder = raw_result_dir_Mitocheck  
    featureL=['Interphases']
    for plate, well in ctrls:
        filename = os.path.join(dataFolder, plate,"hdf5", '{}.ch5')
        try:
            featureL,classes, frameLotC= importTargetedFromHDF5(filename, plate, well,featureL,primary_channel_name=primary_channel_name,
                                                                secondary=False)
        except ValueError:
            print "Error at loading from hdf5 ", plate, well
            continue
        if newFrameLot == None:
            newFrameLot = frameLotC 
        else: 
            newFrameLot.addFrameLot(frameLotC)
            
            featureL = list(featureL)
    if classes is not None:
        featureL.extend(filter(lambda x: x not in featureL, classes))
        
    print'--------------', featureL
    totalResult=defaultdict(dict)
    
    for plate in frameLotC.lstFrames:
        for well in frameLotC.lstFrames[plate]:
            result={'object_count':[], 'cell_count':[]}
            
            result.update({'{}_ch1'.format(el):[] for el in featureL})
            
            for frame_nb in frameLotC.lstFrames[plate][well]:
                frame = frameLotC.lstFrames[plate][well][frame_nb]
                
            #for each frame, the cell count = object count - [artefacts + out of focus objects]. Only possible if the classification was computed
                if classes is not None:
                    bincount = np.bincount(frame.classification, minlength=len(classes))
                    if np.sum(bincount)==0:
                        print "No info frame {} well {} plate {}".format(frame_nb, well, plate)
                        continue
                    out_of_focus = bincount[np.where(classes==focus_classname)] if focus_classname in classes else 0
#NB: no more smallUnidentified after classifier update on 2/20/15
                    artefacts = bincount[np.where(classes==artefact_classname)] if artefact_classname in classes else 0
                    result["cell_count"].append( frame.centers.shape[0] - (out_of_focus + artefacts))
                    
            #for each frame, the object count = all objects = cells+out of focus objects + artefacts
                result["object_count"].append(frame.centers.shape[0])
                
            #other parameters are specified in the setting file    
                for el in featureL:
                    result['{}_ch1'.format(el)].append(bincount[np.where(classes==el)]/float(np.sum(bincount)))
            
            result["initObjectCount"]=result["object_count"][0]
            result["initCellCount"]=result["cell_count"][0]
            result["proliferation"]=result["cell_count"][-1]/float(result["cell_count"][0])

            for el in classes:
                result['{}_ch1'.format(el)]=np.array(result['{}_ch1'.format(el)])[:,0]
                
            if int(well[-2:])==1:
                totalResult[plate][int(well[:-3])].update(result)
            else:
                totalResult[plate][int(well[:-3])]['pos_02']=result
            
            
    return totalResult