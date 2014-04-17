import os, pdb, sys
import numpy as np
from importPack import imp, EVENTS, FEATURE_NUMBER, nb_folds#, c_list
from dataPack import treatments, joining, classify
import cPickle as pickle
import test, tracking
from PyPack import fHacktrack2
import gc
from trajPack import d_ctrl

count = 0
monOutput = ""
loadingFolder = "../results/"
'''
Use function count() at the end of this file to count the content of the training set (how many for each move type).

Parameter: phenoOnly. Default value: False. Put it to True if you want only the content of the training set from
pheno experiments (as opposed to control experiments).

Output: a event number x 2 array. First column is number of links, second column number of events of each type.
So a merge will add +2 in the first column, +1 in the second one. 
'''

def gettingSolu(plate, phenoOnly=False):
    global monOutput
    global loadingFolder
    
    print "current directory ", os.getcwd()
    
    #print os.getcwd()
    newFrameLot = None

    print plate
    listW = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate)
    for well in listW:
        well=well[:-5]
        print well
        if phenoOnly:
            if well not in d_ctrl[plate]:
                print "------------------------PAS PRIS------------------------"
                continue
                
            else:
                print "CONTROL"
                
        filename = '/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate+"/"+well+".hdf5"
        filenameT = '/media/lalil0u/New/workspace2/Tracking/data/TSinclMD/PL'+plate+"___P"+well+"___T00000.xml"

        monOutput+="plate = "+plate+",well = "+well+"\n"
        #ajout du frameLot et du tabF
        frameLotC, tabFC = test.gettingRaw(filename, filenameT, plate, well)
        if newFrameLot == None:
            newFrameLot = frameLotC 
        else: newFrameLot.addFrameLot(frameLotC)
        
    predictions = {}
    
    return newFrameLot, predictions



def movement(uplet, s=None, tC=None, incomplet = False):
    if tC != None:
        if incomplet:
            up= filter(lambda x: tC in x, s[1])[0]
            ind = s[1].index(up)
            up = s[0][ind]
            return movement(up, None,s[1][ind])

        t = tC
        t.sort()
        r=None
        #print uplet, t
        if uplet==[-1]:
            r=1
        elif len(uplet)>1:
            if t==[-1]:
                r=2
            else:
                r= 3
        elif len(uplet) ==1:
            l = len(t)
            if l>1:
                r= 4
            if l ==1:
                if t==[-1]:
                    r= 2
                else:
                    r= 0
        else:
            raise
        #print "resultat :", r
        return r, t
    else:
        z = s
        tu = filter(lambda x: uplet in x, s[0])[0]
        return movement(tu, s= None, tC= z[1][z[0].index(tu)])

def formulation2(pred):
    #j'ai une liste de trajectoires sous la forme [(frame, label)]
    result = {}
    for traj in pred:
        t = traj.lstPoints.keys(); t.sort()
        if t==[]:
            continue
        k=0; continuer = True
        while continuer:
            (frame, label)=t[k]
            
            if frame not in result:
                result[frame]=[[],[]]
            try:
                targ = t[k+1][1]
            except IndexError:
                targ = -1
            k+=1
            if k==len(t):
                continuer = False 
            if label not in result[frame][0] or targ!=-1:
                result[frame][0].append(label)
                result[frame][1].append(targ)
            
                
    return result

def formulation(hypotheses):
    res = [[], []]
    app = []
    for hyp, t in hypotheses:
        if hyp[0]==-1:
            app.append(hyp[1])
        else:
            try:
                h0 = list(hyp[0])
            except TypeError:
                h0 = [hyp[0]]
            try:
                h1 = list(hyp[1])
            except TypeError:
                h1 = [hyp[1]]
                
            res[0].append(h0)
            res[1].append(h1)
    res[0].append([-1])
    res[1].append(app)
    return res

def filtre(frame):
    source= frame.source ; target = frame.target
    s = [] ; t=[]
    for i in range(len(source)):
        if len(source[i])>4 or len(target[i])>4:
            continue
        else:
            s.append(source[i]); t.append(target[i])
    return s, t

            
def calcAccuracy(newFrameLot):
    compteur = [[0,0] for x in range(len(EVENTS))]

    for plate in newFrameLot.lstFrames:
        for well in newFrameLot.lstFrames[plate]:
            for index in newFrameLot.lstFrames[plate][well]:
                
                sourceC, targetC = filtre(newFrameLot.lstFrames[plate][well][index])
                i=0
                for uplet in sourceC:
                    t, ttarg = movement(uplet, s = None, tC=targetC[i])

                    #compteur[0]=nb total de liens, compteur[1]=nb total de matchings 
                    compteur[t][0]+=max(len(uplet), len(ttarg)); compteur[t][1]+=1 
                    i+=1
    print compteur
    return np.array(compteur)

def filtre2(sol, tdict, plate):
    l = sol.lstSolutions
    r = filter(lambda x: x.well in tdict[plate] and x.index in tdict[plate][x.well], l)
    sol2 = joining.Solutions(solution = None, lstSolutions = r)
    return sol2
            
def count(phenoOnly = False): 
    compteurT = [[0,0] for x in range(len(EVENTS))]
    
    for plate in os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw/'):
        nF, predictions = gettingSolu(plate, phenoOnly)
        if nF==None:
            continue
        compteur = calcAccuracy(nF)
        
        compteurT+=compteur
        print compteurT
        del nF ; del predictions ; gc.collect()
    
    print compteurT
        
    return compteurT
    #print>>outputTotal,"TIME TIME TIME TIME TIME", time.clock()