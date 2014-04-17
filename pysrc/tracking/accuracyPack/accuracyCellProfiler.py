#from PIL import Image
import vigra, h5py, random
import numpy as np
import os, pdb, time
import cPickle as pickle
import test
import subprocess
from itertools import product
from importPack import EVENTS, nb_big_folds, nb_folds
from trajPack import d_ctrl
from accuracyPack import accuracy
from accuracy import movement
from scipy import spatial
import warnings

index_file='/media/lalil0u/New/workspace2/Tracking/data/CeproPredictions/INDEXES.h5'

def gettingSolu(plate, well):
    global loadingFolder
    newFrameLot = None

    print "Loading data from ", plate, well
    filename = '/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate+"/"+well+".hdf5"
    filenameT = '/media/lalil0u/New/workspace2/Tracking/data/TSinclMD/PL'+plate+"___P"+well+"___T00000.xml"

    #ajout du frameLot et du tabF
    frameLotC, tabFC = test.gettingRaw(filename, filenameT, plate, well)
    if newFrameLot == None:
        newFrameLot = frameLotC 
    else: newFrameLot.addFrameLot(frameLotC)
#    
    return newFrameLot

def gettingDataIndexes(filename):
    data=h5py.File(filename)
    tabs=[];
    data.visit(tabs.append)
    
    tabsPlate = filter(lambda x: 'Metadata_Plate' in x and 'data' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs);
    tabsWell = filter(lambda x: 'Metadata_Well' in x and 'data' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs);
    plates =np.array(data[tabsPlate[0]])
    plList = np.unique(plates)
    wells = np.array(data[tabsWell[0]])
    weList = np.unique(wells) 
    if not plates.shape==wells.shape:
        raise
    metadata = {}#np.vstack((plates, wells))#ok, size 2, 2351
    realNames = d_ctrl.keys()
    for p, w in product(plList, weList):
        tt= filter(lambda x: x in np.where(wells==w)[0], np.where(plates==p)[0])
        if tt!=[]:
            #+1 for it to be the line number as expected by CellProfiler
            well='00'+str(w)+'_01'
            plate = corresp(realNames,str(p))
            
            metadata.update({(plate, well):(min(tt)+1, max(tt)+1)})
    return metadata

def corresp(realNames, p):
    for name in realNames:
        if name[2:6]==p:
            return name
    return None
    
def gettingCepro(filename):
    data=h5py.File(filename)
    tabs=[]; parents=[]; children=[]
    data.visit(tabs.append)
    if tabs==[]:
        return None, None
    
    tabsCenterXIndx = filter(lambda x: 'Location_Center_X' in x and 'index' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs);
    tabsCenterYIndx = filter(lambda x: 'Location_Center_Y' in x and 'index' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs);
    if not np.all(np.array(data[tabsCenterXIndx[0]])==np.array(data[tabsCenterYIndx[0]])):
        raise
    tabsCenterX = filter(lambda x: 'Location_Center_X' in x and 'data' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs);
    tabsCenterY = filter(lambda x: 'Location_Center_Y' in x and 'data' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs);
    
    tabsC=filter(lambda x: 'Relationship' in x and 'Second' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs); tabsC.sort()
    tabsP=filter(lambda x: 'Relationship' in x and 'First' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs); tabsP.sort()
    f=True
    for name in tabsP:
        if f:
            t=np.array(data[name])
            firstLine =np.min(t)
            t=t-firstLine+1
            parents.append(t)
            f=False
        else:
            parents.append(data[name])
    f=True
    for name in tabsC:
        if f:
            t=np.array(data[name])
            t=t-firstLine+1
            children.append(t)
            f=False
        else:
            children.append(data[name])
    parents = np.vstack((np.array(parents[0]), np.array(parents[1]))); parents=parents.transpose()
    children = np.vstack((np.array(children[0]), np.array(children[1]))); children=children.transpose()
    
    Object_NumberIndx = filter(lambda x: 'Number_Object_Number' in x and 'index' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs)
    Parent_NumberIndx = filter(lambda x: 'TrackObjects_Parent' in x and 'index' in x and isinstance(data[x], h5py._hl.dataset.Dataset), tabs)
    ParentIndx = np.array(data[Parent_NumberIndx[0]])
    
    max = np.max(ParentIndx[np.where(ParentIndx[:,0]<0)])
    ParentIndx2 = ParentIndx[np.where(ParentIndx[:,0]>-1)]
    ParentIndx2[:,1:]=ParentIndx2[:,1:]-max
    if not np.all(np.array(data[Object_NumberIndx[0]])==ParentIndx2):
        raise
    if not np.all(np.array(data[Object_NumberIndx[0]])==np.array(data[tabsCenterYIndx[0]])):
        raise
#    pdb.set_trace()
    Object_Number = np.array(data[filter(lambda x: 'Number_Object_Number' in x and 'data' in x and isinstance(data[x], 
                            h5py._hl.dataset.Dataset), tabs)[0]], dtype=int)
    Parent_Number = np.array(data[filter(lambda x: 'TrackObjects_ParentObjectNumber' in x and 'data' in x and isinstance(data[x], 
                            h5py._hl.dataset.Dataset), tabs)[0]], dtype=int)
    
    trees = {}
    output=[[],[]]
    frames=[]
    FrameIndex = np.array(data[Object_NumberIndx[0]])
    FrameIndex[:,0]=FrameIndex[:,0]-np.min(FrameIndex[:,0])+1
    centersX = data[tabsCenterX[0]]; centersY= data[tabsCenterY[0]]
#    pdb.set_trace()
    print "Building trees"
    for k in range(FrameIndex.shape[0]):
        
        frames.extend([k+1 for ttt in range(FrameIndex[k][1],FrameIndex[k][2])])
        #if k==92: 
        centers = zip(centersX[FrameIndex[k][1]:FrameIndex[k][2]], centersY[FrameIndex[k][1]:FrameIndex[k][2]])
        tree = spatial.cKDTree(centers, leafsize = 10)
        trees.update({k+1:tree})
        output[0].extend(Parent_Number[FrameIndex[k][1]:FrameIndex[k][2]])
        output[1].extend(Object_Number[FrameIndex[k][1]:FrameIndex[k][2]])
    output=np.vstack((frames, output[0], output[1]))
 #   pdb.set_trace()
    
    #pour trouver les appear
    appP = [[],[]]; appC=[[],[]]
    for k in range(output.shape[1]):
        if output[1][k]==0 and output[0][k]>1:
            if output[2][k] not in children[np.where(children[:,0]==output[0][k])][:,1]:
                appP[0].append(output[0][k]-1)
                appP[1].append(0)
                appC[0].append(output[0][k])
                appC[1].append(output[2][k])

    appP = np.vstack((np.array(appP[0]), np.array(appP[1]))); appP=appP.transpose()
    appC = np.vstack((np.array(appC[0]), np.array(appC[1]))); appC=appC.transpose()
    parents = np.vstack((parents, appP))
    children = np.vstack((children, appC))
    #pdb.set_trace()
    print "Merges and split"
    out={}; lastChildrenC=output[:,np.where(output[0]==1)][2][0]
    for k in range(1,np.max(parents[:,0]+1)):
        out.update({k:[[],[]]})
        parentsC=parents[np.where(parents[:,0]==k)][:,1]
        childrenC=children[np.where(children[:,0]==k+1)][:,1]
        #pdb.set_trace()
        merge=list(np.bincount(childrenC)); split=list(np.bincount(parentsC))
        for cell in range(len(split)):
#            if cell==0:
#                pdb.set_trace()
            if split[cell]==0:
                continue
            if split[cell]>1:
                out[k][0].append(changeApp(cell))
                out[k][1].append(list(childrenC[np.where(parentsC==cell)]))
            elif split[cell]==1:
                #try:
                child = list(childrenC[np.where(parentsC==cell)])
                if merge[child[0]]==1:
                    out[k][0].append(changeApp(cell))
                    out[k][1].append(changeApp(child))
                elif merge[child[0]]>1:
                    #if the merge is already recorded
                    if child in out[k][1]:
                        continue
                    out[k][0].append(list(parentsC[np.where(childrenC==child)]))
                    out[k][1].append(changeApp(child))
#                except:
#                    pdb.set_trace()
        #if k>1:
        disp = filter(lambda x: x not in parentsC, lastChildrenC)
        for c in disp:
            out[k][0].append(changeApp(c))
            out[k][1].append([-1])
        lastChildrenC=childrenC
        #if k==1: pdb.set_trace()
    return out, trees

def CeproLabelsTraduction(frameC, source, treeC):
    rSource = []; pasTrouves=[]
    for uplet in source:
        rUplet = []
        if uplet==[-1]:
            rUplet=uplet
            rSource.append(rUplet)
            continue
        for u in uplet:
            center = list(frameC.getCenter(u))
            dist, i = treeC.query(center, k=2, distance_upper_bound = 40)
            try:
                tradCP=filter(lambda x: x[1]<100, zip(i,dist))[0][0]+1
            except IndexError:
                pasTrouves.append(u)
                continue
            #print "traduction de",u, ":", tradCP
            rUplet.append(tradCP)
        rSource.append(rUplet)
        
    return rSource, pasTrouves

def changeApp(label):
    if label==0:
        return [-1]
    elif type(label)==list:
        return label
    else: return [label]

def calcAccuracy(newFrameLot, output, trees):
    visu={}
    acc = np.zeros(shape=(len(EVENTS), 2), dtype=np.int)
    acc2 = np.zeros(shape=(len(EVENTS), len(EVENTS)), dtype=np.int)
    compteur = [[0,0] for x in range(len(EVENTS))]
    for plate in newFrameLot.lstFrames:
        outPasTrouves = '\n plate {}'.format(plate)
        if plate not in visu:
            visu[plate]={}
        for well in newFrameLot.lstFrames[plate]:
            print "Accuracy for", plate, well
            outPasTrouves += 'well {}'.format(well)
            if well not in visu[plate]:
                visu[plate][well]={}
            for index in newFrameLot.lstFrames[plate][well]:
                #print index,
                if index not in visu[plate][well]:
                    visu[plate][well][index]=[]
                if index+1 not in visu[plate][well]:
                    visu[plate][well][index+1]=[]

                
            #DATA CELL PROFIleR
                try:
                    outputC =output[index+1]
                except KeyError:
                    continue
                else:
                    treeC = trees[index+1]
                    nextTree = trees[index+2]
            #DATA TRAINING SET
                frameC= newFrameLot.lstFrames[plate][well][index]
                nextFrame = newFrameLot.lstFrames[plate][well][index+1]
                sourceC, targetC =accuracy.filtre(frameC)
                
                sourceC, pasTrouvesS = CeproLabelsTraduction(frameC, sourceC , treeC)
                targetC, pasTrouvesT = CeproLabelsTraduction(nextFrame, targetC, nextTree)
                if pasTrouvesS!=[]:
                    print '########################################################################PAS TROUVES DANS LA SOURCE', pasTrouvesS
                    visu[plate][well][index].append(pasTrouvesS)
                if pasTrouvesT!=[]:
                    print '#######################################################################PAS TROUVES DANS LA CIBLE', pasTrouvesT
                    visu[plate][well][index+1].append(pasTrouvesS)
                if pasTrouvesT!=[] or pasTrouvesS!=[]:
                    outPasTrouves += '\n index {}'.format(index)
                    outPasTrouves+=str(pasTrouvesS)+' target '+str(pasTrouvesT)
                i=0
                
#                for k in range(len(outputC[0])):
#                    print 's', outputC[0][k], outputC[1][k],
                if [] in sourceC:
                    for k, s in enumerate(sourceC):
                        if s==[]:
                            del sourceC[k]
                            del targetC[k]
                for uplet in sourceC:
#                    if plate =='LT0024_41' and index==76 and uplet==[236]:
#                        pdb.set_trace()
                    t, ttarg = movement(uplet, s = None, tC=targetC[i])
                    if t==None:
                        if targetC[i]==[]:
                            if pasTrouvesT==[]:
                                pdb.set_trace()
                            else:
                                continue
                    OUT = "truth {},{},{}".format(uplet, t, ttarg)
                    #compteur[0]=nb total de liens, compteur[1]=nb total de matchings 
                    
                    compteur[t][0]+=max(len(uplet), len(ttarg)); compteur[t][1]+=1 
                    first = True
                    for u in uplet:
                        try:
                            p, ptarg = movement(u, s=outputC, tC = None)
                        except IndexError:
                            #here it means that there is no appear in the predicted sequence
                            if u==-1:
                                for tt in ttarg:
                                    z, ztarg = movement (u, s=outputC, tC = tt, incomplet = True)
                                    acc2[t][z]+=1
                                    continue
                            else:
                                pdb.set_trace()

                        OUT+="predit {},{}".format(p, ptarg)
                        
                        if t==p:
                            test = (ttarg==ptarg) if t!=1 else (ttarg in ptarg)
                            if test:
                                #le matching est le bon et la cible egalement : j'enregistre acc[0]=nb de bons liens, acc[1]=nb de bons matchings
                                acc[t][0]+=max(1, len(ttarg))
                                if first: acc[t][1]+=1
                                first = False
                                #acc2[t][p]+=1
                            else:
#                                    donc la j'enregistre les erreurs, ici il y en a 
                                if t==1:
                                    #print OUT, '------------------------FALSE MOVE'
                                    for tt in ttarg:
                                        z, ztarg = movement (u, s=outputC, tC = tt, incomplet = True)
                                        acc2[t][z]+=1
                                else:
                                    #print OUT,"------------------------FALSE TARGET"
                                    for tt in ttarg:
                                        if tt in ptarg:
                                            acc[t][0]+=1
                                        else:
                                            z, ztarg = movement(u, s=outputC, tC=tt, incomplet=True)
                                            acc2[t][z]+=1
                        else:
                            #print OUT,'------------------------FALSE MOVE'
                            if len(ttarg)>1:
                                if first:
                                    for tt in ttarg:
                                        z, ztarg = movement (u, s=outputC, tC = tt, incomplet = True)
                                        acc2[t][z]+=1
                                first = False
                            else:
                                acc2[t][p]+=1
                            
                    i+=1
    compteur =np.array(compteur)
    print acc, '\n', acc2, '\n', compteur
    return acc, acc2, np.array(compteur), outPasTrouves
            
                
                
                

def makeTif(pathSearch, pathSave, filename):  
    frameNum = os.path.splitext(filename)[0][-5:]; begName = os.path.splitext(filename)[0][:-5]; ext = os.path.splitext(filename)[1]
    newFrameNum = '{:0>5}'.format(int(frameNum)*30)
    try:
        im=vigra.readImage(os.path.join(pathSearch, filename))
    except IOError:
        print "File pbl with file ", filename
        return 0
    else:
        im2=vigra.Image(im, dtype=np.uint8)
        im2[im2>0]=1
        im2.writeImage(os.path.join(pathSave, 'Mask_'+begName+newFrameNum+ext))
        return 1
    
if __name__=='__main__':
    path = "/media/lalil0u/New/projects/groundtruth/results/all2"
    #pour enregistrer la segmentation
    pathSave = "/media/lalil0u/New/data/ground_truth_tracking/CecogSegmentation"
    pathData = '/media/lalil0u/New/workspace2/Tracking/data/CeproPredictions'
    nom='accuracy_Cepro_TestSet_2'; param =['split_alternative_cost', 'max_split_score', 'merge_alternative_cost', 'max_merge_score']
    
    first_pass=[(40,80),(60,100), (100, 140)]
    second_pass = [(100,300),(150,300), (200, 300)]
    third_pass = [(40, 140),(60, 140), (80, 140)]
    
    split_p=[]; split_p.extend(first_pass); split_p.extend(third_pass); split_p.extend(second_pass)
    first_pass_m = [(40,100),(40,120), (20, 80), (60,80), (70,80)]
    merge_p = [(40,80)]; merge_p.extend(first_pass_m)
    split_param_list=[] ; merge_param_list=[]; param_list=[] 
    metadata=gettingDataIndexes(index_file)
#    for j in range(len(split_p)):
#        parameter = 'm{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format(first_pass[0], split_p[j])
#        split_param_list.append(parameter)
#        print "PARAMETER S", parameter
#    for j in range(len(merge_p)):
#        parameter = 'm{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format(merge_p[j], first_pass[2])
#        merge_param_list.append(parameter)
#        print "PARAMETER M", parameter
#    merge_param_list.extend([ 'm{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (120,140)), 
#                             'm{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (150,200))])
    merge_param_list.extend([ #'m{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (120,140)), 
                             'm{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (100,120)),
                             'm{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (80,100)),
                             'm{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (100,140)),
                             'm{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (120,140))
                             ])

    #merge_param_list=['m{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (100,140))]
    print merge_param_list

#CROSS VALIDATION
#    N=25; CVfolds = []
#    r=list(random.sample(range(25), 25))
#    for k in range(nb_big_folds):
#        CVfolds.append(list(r[k*5:(k+1)*5]))                     
#    outputCV=open("CVfolds_2.pkl", 'w'); pickle.dump(CVfolds, outputCV); outputCV.close()
    outputCV=open("CVfolds_2.pkl", 'r'); 
    CVfolds=pickle.load(outputCV); outputCV.close()
    print CVfolds
    fold=True; testSet=True
    output2=open(nom+'TOTAL.pkl', 'w')
    output3=open(nom+'TOTAL.txt', 'w')
    listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
    if not fold:
        r = np.zeros(shape=(len(param_list),len(EVENTS), 2), dtype=np.int)
        r2 = np.zeros(shape=(len(param_list),len(EVENTS), len(EVENTS)), dtype=np.int)
        compteur = np.zeros(shape=(len(EVENTS),2))
        first=True
        for plate in listD:
            listW = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate)
            for well in listW:
                well=well[:-5]
                newFrameLot = gettingSolu(plate, well)
                first=True
                for i, param in enumerate(param_list):
                    print "param ", param
                    CeproFilename = os.path.join('/media/lalil0u/New/workspace2/Tracking/data/CeproPredictions/', param, plate, well+".h5")
                    #try:
                    output, trees = gettingCepro(CeproFilename)
                    if output==None:
                        warnings.warn("Le fichier est vide")
                        print "##################### FICHIER VIDE"
                        continue
#                    except:
#                        print "FILE NOT FOUND"
#                        continue
        
                    acc, acc2, compt, pasTrouves = calcAccuracy(newFrameLot, output, trees)
                
                    output3.write(pasTrouves)
                    #pdb.set_trace()
                    r[i]+=acc; r2[i]+=acc2; 
                    if first: compteur+=compt
                    first=False
        result = []
        for i, param in enumerate(param_list):
            #recall
            recall = np.sum(r[i], 0)/np.array(np.sum(compteur, 0), dtype=float)
                
            #pourcentages de reussites dans les differents mouvements, en nb de mouvements
            p = r[i]/np.array(compteur, dtype=float)
                
            #enfin la precision par mouvement et en nb de liens
            precision =np.array(np.sum(r[i], 0)[0], dtype=float)/(np.sum(r[i], 0)[0]+np.sum(r2[i], 0))
            #pdb.set
            result.append((recall, p, precision))
        pickle.dump([param_list, result, r, r2, compteur], output2)
    else:
        r = np.zeros(shape=(nb_big_folds, len(merge_param_list)+len(split_param_list),len(EVENTS), 2), dtype=np.int)
        r2 = np.zeros(shape=(nb_big_folds,len(merge_param_list)+len(split_param_list),len(EVENTS), len(EVENTS)), dtype=np.int)
        compteur = np.zeros(shape=(nb_big_folds,len(EVENTS),2))
        first=True
        
        for k in range(nb_big_folds):
            print 'FOLD ', k
            l=0
            for plate in listD:
                listW = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate)
                for well in listW:
                    well=well[:-5]
                    print l
                    if l not in CVfolds[k]:
                        l+=1
                        print "PAS PRIS"
                        continue
                    else:
                        l+=1 
                    newFrameLot = gettingSolu(plate, well)
                    first=True
                    for i, param in enumerate(merge_param_list):
                        print "param ", param
                        CeproFilename = os.path.join('/media/lalil0u/New/workspace2/Tracking/data/CeproPredictions/', param, plate, well+".h5")
                        output, trees = gettingCepro(CeproFilename)
                        
                        if output==None:
                            warnings.warn("Le fichier est vide")
                            continue
#                        pdb.set_trace()
                        acc, acc2, compt, pasTrouves = calcAccuracy(newFrameLot, output, trees)
                    
                        output3.write(pasTrouves)
                        #pdb.set_trace()
                        if not testSet:
                            r[k][i]+=acc; r2[k][i]+=acc2;
                        else:
                            r[0][i]+=acc; r2[0][i]+=acc2;
                        if first: 
                            if not testSet:
                                compteur[k]+=compt
                            else:
                                compteur[0]+=compt
                        first=False
                    for i, param in enumerate(split_param_list):
                        print "param ", param
                        CeproFilename = os.path.join('/media/lalil0u/New/workspace2/Tracking/data/CeproPredictions/', param, plate, well+".h5")
                        #try:
                        output, trees = gettingCepro(CeproFilename)
                        if output==None:
                            warnings.warn("Le fichier est vide")
                            continue
#                        except:
#                            print "FILE NOT FOUND"
#                            continue
            
                        acc, acc2, compt, pasTrouves = calcAccuracy(newFrameLot, output, trees)
                    
                        r[k][i+len(merge_param_list)]+=acc; r2[k][i+len(merge_param_list)]+=acc2; 
        outputRawResults = open(nom+'RawResults.pkl', 'w')
        pickle.dump([merge_param_list, split_param_list, r, r2, compteur], outputRawResults);  outputRawResults.close()
        splitR=[]; mergeR=[]
        if not test:
            for k in range(nb_big_folds):
            #MERGE
                mergeR.append([])
                for i, param in enumerate(merge_param_list):
                    #recall
                    recall = np.sum(r[k][i], 0)/np.array(np.sum(compteur[k], 0), dtype=float)
                        
                    #pourcentages de reussites dans les differents mouvements, en nb de mouvements
                    p = r[k][i]/np.array(compteur[k], dtype=float)
                        
                    #enfin la precision par mouvement et en nb de liens
                    pp=np.sum(np.array(r[k][i][:,0], dtype=float))/np.sum((r[k][i][:,0]+np.sum(r2[k][i], 0)))
                    precision =np.array(r[k][i][:,0], dtype=float)/(r[k][i][:,0]+np.sum(r2[k][i], 0))
                    #pdb.set
                    #mergeR[k].append((recall[1], p[3][1], precision[3]))
                    mergeR[k].append((recall[1], p[3][1], precision[3]))
            #SPLIT
                splitR.append([])
                for i, param in enumerate(split_param_list):
                    #recall
                    recall = np.sum(r[k][i+len(merge_param_list)], 0)/np.array(np.sum(compteur[k], 0), dtype=float)
                        
                    #pourcentages de reussites dans les differents mouvements, en nb de mouvements
                    p = r[k][i+len(merge_param_list)]/np.array(compteur[k], dtype=float)
                        
                    #enfin la precision par mouvement et en nb de liens
                    precision =np.array(r[k][i+len(merge_param_list)][:,0], dtype=float)/(r[k][i+len(merge_param_list)][:,0]
                                                                                                +np.sum(r2[k][i+len(merge_param_list)], 0))
                    #pdb.set
                    splitR[k].append((recall[1], p[4][1], precision[4]))
        else:
            k=0
            for i, param in enumerate(merge_param_list):
                #recall
                recall = np.sum(r[k][i], 0)/np.array(np.sum(compteur[k], 0), dtype=float)
                    
                #pourcentages de reussites dans les differents mouvements, en nb de mouvements
                p = r[k][i]/np.array(compteur[k], dtype=float)
                    
                #enfin la precision par mouvement et en nb de liens
                pp=np.sum(np.array(r[k][i][:,0], dtype=float))/np.sum((r[k][i][:,0]+np.sum(r2[k][i], 0)))
                precision =np.array(r[k][i][:,0], dtype=float)/(r[k][i][:,0]+np.sum(r2[k][i], 0))
                #pdb.set
                #mergeR[k].append((recall[1], p[3][1], precision[3]))
                mergeR.append((recall, pp, p, precision))
        #pickle.dump([(merge_param_list, mergeR),(split_param_list, splitR)], output2)
        pickle.dump({'merge':[merge_param_list, mergeR], 'split':[split_param_list, splitR]}, output2)
    output2.close()    
    output3.close()
        
        
        
        
        
        
        