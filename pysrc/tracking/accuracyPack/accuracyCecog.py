import os, pdb, time, sys
import numpy as np
from importPack import imp, EVENTS, FEATURE_NUMBER, nb_folds#, c_list
from dataPack import treatments, joining, classify
import cPickle as pickle
import test, tracking
from joblib import Parallel, delayed
from PyPack import fHacktrack2
import gc

count = 0
monOutput = ""
loadingFolder = "../results/"

#donc cette fonction va recuperer le training set VRAI ie incl merge divides, elle va charger un modele et mettre qqpart la prediction du modele,
# et elle va retourner le tout pour faire un calcul d'accuracy

def gettingSolu(plate):
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
        filename = '/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate+"/"+well+".hdf5"
        filenameT = '/media/lalil0u/New/workspace2/Tracking/data/TSinclMD/PL'+plate+"___P"+well+"___T00000.xml"

        monOutput+="plate = "+plate+",well = "+well+"\n"
        #ajout du frameLot et du tabF
        frameLotC, tabFC = test.gettingRaw(filename, filenameT, plate, well)
        if newFrameLot == None:
            newFrameLot = frameLotC 
        else: newFrameLot.addFrameLot(frameLotC)
        
    predictions = {}

    for well in listW:
        well=well[:-5]
        print well
        filenameD = '/media/lalil0u/New/workspace2/Tracking/data/CecogPredictions/'+plate+'/tracking_graph___P'+well+'.dot'
        predictions.update({well: lireDot(filenameD)})
    
    return newFrameLot, predictions

def lireDot(nomDot):
    fDot = open(nomDot, "r")

    mesLignes=fDot.readlines()
    fDot.close()
    del mesLignes[0:3]
    mesLignes.reverse()

    continuer = True
    ensembleTraj = fHacktrack2.ensTraj()
    print ensembleTraj.numMitoses
    print ensembleTraj.lstTraj

    trajCourante = fHacktrack2.trajectoire(0,0,0,0,0)
    frameSuivante =-1
    idSuivante =-1

    while continuer:
        ligneCourante=mesLignes.pop()
        if ("#AAAAAA" in ligneCourante):
            continuer = False
            continue

        decomp = ligneCourante.split(" ")
        decomp[0]=decomp[0][1:-1]

        frameCourante = int(decomp[0].split("_")[0])
        idCourante = int(decomp[0].split("_")[1])
        decomp[2]=decomp[2][1:-1]
        
#        if idCourante==105 and frameCourante ==1:
#            pdb.set_trace()

        if (int(frameCourante) != int(frameSuivante) or int(idCourante) != int(idSuivante)):
            if frameSuivante!=-1 and idSuivante!=-1:
                trajCourante.ajoutPoint(0,0, frameSuivante, idSuivante)
              
            trajectoires = ensembleTraj.trouver(trajCourante.numCellule)
            lstPointsTotale = []
            for trajec in trajectoires.lstTraj:
                for element in trajec.lstPoints.keys():
                    #print element
                    lstPointsTotale.append(element)

            #if lstPointsTotale <> []: print lstPointsTotale
            if (frameCourante, idCourante) in lstPointsTotale or (frameCourante, idCourante) in trajCourante.lstPoints.keys():
                #print "mitose", frameCourante, idCourante
                trajCourante.mitose = 1
                ensembleTraj.numMitoses +=1
                if (frameCourante, idCourante) in trajCourante.lstPoints.keys():
                    trajCourante.ajoutPoint(0, 0, frameCourante, idCourante)
                ensembleTraj.ajouter(trajCourante)
                num = trajCourante.numCellule
  #              print num
                trajCourante.supp()
                trajCourante = fHacktrack2.trajectoire(num, 0,0, frameCourante, idCourante)
                trajCourante.mitose = 1

            else:
   #             print "nouvelle trajectoire", frameCourante, idCourante
                ensembleTraj.ajouter(trajCourante)
                num = trajCourante.numCellule+1

                trajCourante.supp()

                trajCourante = fHacktrack2.trajectoire(num, 0,0, frameCourante, idCourante)             
#                print trajCourante.lstPoints

        else:
 #           print "ajout point"
            trajCourante.ajoutPoint(0, 0, frameCourante, idCourante)
  #          print trajCourante.lstPoints
     

        frameSuivante = int(decomp[2].split("_")[0])
        idSuivante = int(decomp[2].split("_")[1])

    ensembleTraj.ajouter(trajCourante)

    return ensembleTraj

def findNext(pred, label):
    result = []
    k=0
    for l in pred[0]:
        if l==label and pred[1][k] not in result:
            result.append(pred[1][k])
        k+=1
#    if len(result)>1:
#        pdb.set_trace()
    return result

def findPrevious(pred, label):
    result = []
    k=0
    for l in pred[1]:
        if l in label and pred[0][k] not in result:
            result.append(pred[0][k])
        k+=1
#    if len(result)>1:
#        pdb.set_trace()
 #   pdb.set_trace()
    return result


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

def postProcessing(r, r2, compteur):
    #recall
    recall = [0 for z in range(len(r))]
    for z in range(len(r)):
        recall[z]=np.sum(r[z], 0)/np.array(np.sum(compteur, 0), dtype=float)
        
    #pourcentages d'erreurs dans les differents mouvements
    p = [0 for y in range(len(r))]
    for z in range(len(r)):
        p[z]=np.sum(r2[z], 1)/compteur.transpose()[1]*100
        
    #enfin la precision
    precision = [0 for y in range(len(r))]
    pdb.set_trace()
    for y in range(len(r)):
        precision[y]=np.array(r[y].transpose()[1], dtype=float)/(r[y].transpose()[1]+np.sum(r2[y], 1))

    return [r, recall, r2, p, precision, compteur]
            
def calcAccuracy(newFrameLot, predictions):
    acc = np.zeros(shape=(len(EVENTS), 2), dtype=np.int)
    acc2 = np.zeros(shape=(len(EVENTS), len(EVENTS)), dtype=np.int)
    compteur = [[0,0] for x in range(len(EVENTS))]

    for plate in newFrameLot.lstFrames:
        for well in newFrameLot.lstFrames[plate]:
            predGlobal = formulation2(predictions[well].lstTraj)
            for index in newFrameLot.lstFrames[plate][well]:
                
                print plate, well, index
                sourceC, targetC = filtre(newFrameLot.lstFrames[plate][well][index])
                pred = predGlobal[index]
                i=0
                for uplet in sourceC:
                    t, ttarg = movement(uplet, s = None, tC=targetC[i])

                    #compteur[0]=nb total de liens, compteur[1]=nb total de matchings 
                    compteur[t][0]+=max(len(uplet), len(ttarg)); compteur[t][1]+=1 
                    first = True ; ptarg = [] ; previous = {}; predit ={} ; next=[]
                    
                    if uplet==[-1]:
                        previous= findPrevious(pred, ttarg)
                        if previous==[]:
                            previous=[-1]
                        p, pp = movement(previous, s=None, tC =ttarg)
                        if p==1:
                            acc[p][1]+=1; acc[p][0]+=1
                            print "vrai appear"
                        else:
                            acc2[1][p]+=1
                        #pdb.set_trace()
                    elif len(uplet)>1:
                        print "plus d'un", uplet
                        for u in uplet:
                            ptarg= findNext(pred, u)
                            predit.update({u:ptarg})
                            if ptarg not in next:
                                next.append(ptarg)
                        print 'predit :', predit
                    else:
                        for u in uplet:
                            ptarg.extend(findNext(pred, u))
                        previous= findPrevious(pred, ptarg)
                        ptarg.sort()
                        print "uplet ", uplet
                        print 'predit :', previous, ptarg
                        if previous ==[]:
                            #pdb.set_trace()
                            previous=[-1]
                            ptarg = ttarg
#                    if index==26:
#                       
                    if len(uplet)>1:
                        print "BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH", next, ttarg
                        if next==ttarg:
                            print uplet, "predit Cecog ", next
                            acc[t][0]+=max(len(uplet), len(ttarg))
                            if first: acc[t][1]+=1
                            first = False
                            continue
                        else:
                            for u in uplet:
                                print u, predit[u]
                                p, pp=movement([u], s=None, tC =predit[u])
                                acc2[t][p]+=1   
                                print "mouvement predit ", u, predit[u], p
                                
                    elif len(ttarg)>1:
                        if ttarg==ptarg:
                            print "+1 vrai"
                            #le matching est le bon et la cible egalement : j'enregistre acc[0]=nb de bons liens, acc[1]=nb de bons matchings
                            acc[t][0]+=max(len(uplet), len(ttarg))
                            if first: acc[t][1]+=1
                            first = False
                        else:
                            for label in ttarg:
                                previous = findPrevious(pred, [label])
                                if previous==[]:
                                    previous = [-1]
                                p, pp = movement(previous, s=None, tC=[label])
                                print 'BIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII mouvement predit', p
                                acc2[t][p]+=1
                    elif uplet!=[-1]:
                        p, pp = movement(previous, s=None, tC =ptarg)
               #     pdb.set_trace()
                        if t==p:
                        #    print "+1 vrai"
                            #le matching est le bon et la cible egalement : j'enregistre acc[0]=nb de bons liens, acc[1]=nb de bons matchings
                            acc[t][0]+=max(1, len(ttarg))
                            if first: acc[t][1]+=1
                            first = False
                        else:
                            print 'mouvement predit : ', p
                            acc2[t][p]+=1

                    i+=1
    return acc, acc2, np.array(compteur)

def filtre2(sol, tdict, plate):
    l = sol.lstSolutions
    r = filter(lambda x: x.well in tdict[plate] and x.index in tdict[plate][x.well], l)
    sol2 = joining.Solutions(solution = None, lstSolutions = r)
    return sol2
            
def accu(): 
    accT = np.zeros(shape=(len(EVENTS), 2), dtype=np.int)
    acc2T = np.zeros(shape=(len(EVENTS), len(EVENTS)), dtype=np.int)
    compteurT = [[0,0] for x in range(len(EVENTS))]
    
    for plate in os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw/'):
        nF, predictions = gettingSolu(plate)
        acc, acc2, compteur = calcAccuracy(nF, predictions)
        print acc, acc2, compteur
        accT+=acc ; acc2T+=acc2; compteurT+=compteur
        del nF ; del predictions ; gc.collect()
    
    print accT, acc2T, compteurT
        
    return accT, acc2T, compteurT

if __name__ == '__main__':

    nom="accuracy_Cecog2"
    acc, acc2, compteur= accu()
    fi = open(nom+'.pkl', 'w')
    pickle.dump(postProcessing([acc], [acc2], compteur), fi)
    fi.close()
    #print>>outputTotal,"TIME TIME TIME TIME TIME", time.clock()
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
