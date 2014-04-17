import os, pdb, time, sys, argparse
import numpy as np
from PyPack import fHacktrack2
from importPack import imp, EVENTS, FEATURE_NUMBER, nb_folds, nb_big_folds, c_list
from trajPack import d_ctrl
from dataPack import treatments, joining, classify
import cPickle as pickle
import test, tracking
from joblib import Parallel, delayed
import gc, warnings
from accuracyPack.accuracyCecog import loadingFolder

count = 0
monOutput = ""
#loadingFolder = "../results/halfTS/"
ctrl=7; pheno = 6

#donc cette fonction va recuperer le training set VRAI ie incl merge divides, elle va charger un modele et mettre qqpart la prediction du modele,
# et elle va retourner le tout pour faire un calcul d'accuracy

def gettingSolu(plate, pheno_only=False):
    global monOutput
    global FEATURE_NUMBER
    global loadingFolder
    global first
    global ctrl; global pheno
    #tableau qui contiendra toutes les features de tlm pr voir lesquelles contiennent des NaN
    tabF= None
    
    print "current directory ", os.getcwd()
    fichier = open(loadingFolder+"featuresToDelete.pkl", 'r')
    f = pickle.load(fichier)
    fichier.close()
    
    #print os.getcwd()
    newFrameLot = None

    print plate
    listW = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate)
    for well in listW:
        well=well[:-5]
        print well
        if pheno_only:
            if well not in d_ctrl[plate]:
                if pheno>0:
                    print "------------------------PAS PRIS------------------------"
                    pheno-=1
                    continue
                else:
                    print "PHENO"
            else:
                if ctrl>0:
                    print "------------------------PAS PRIS------------------------"
                    ctrl-=1
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
        tabF = tabFC if tabF == None else np.vstack((tabF, tabFC))
    
    #    print "final training set content :"
    #    count, total= newFrameLot.statisticsTraining2()
    #    print count, total
    if newFrameLot is None:
        return None, None
        
    c, f2 = treatments.whichAreNan(tabF)
    #print len(f2.keys()), len(f)
#    if len(f)<len(f2.keys()):
#        pdb.set_trace()
#        
    #if there are features with NaN entries in the predict data but not in the training data
    toZeros = filter(lambda x: x not in f, f2.keys())
    #pdb.set_trace()
    if toZeros !=[]:
        msg="WARNING WARNING WARNING, some features here have NaN entries, and this was not the case in the training set. They are put to 0"
        warnings.warn(msg)
        newFrameLot.zeros(toZeros)
    newFrameLot.clean(f)
    if first:
        FEATURE_NUMBER -=len(f)
    
    print "Feature number", FEATURE_NUMBER
    print "Getting all uplets now"
    print "TIME TIME TIME", time.clock()
    #ICI ON RECUPERE DONC LES SINGLETS ET DOUBLETS AVEC LA VALEUR DU TRAINING DANS CELL.TO SI ILS Y SONT, NONE SINON
    #POUR LE CENTRE ET LES FEATURES C'EST LA MOYENNE DES OBJETS DU SINGLET
    singlets, doublets = newFrameLot.getAllUplets()
    #pour l'instant je ne garde que le passage de l'image 0 a 1
    print "TIME TIME TIME after getting all uplets", time.clock()
    print "Joining uplets now"
    
    solutions = joining.j(singlets, doublets, FEATURE_NUMBER, training = False)
    print "TIME TIME TIME after joining", time.clock()
    
    return solutions, newFrameLot


def movement(uplet, s=None, tC=None, incomplet = False):
    if tC != None:
        if incomplet:
            #pdb.set_trace()
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

def finalPostProcessing(r, r2, compteur):
    #recall
    recall = np.sum(np.sum(r, 0), 0)/np.array(np.sum(np.sum(compteur, 0), 0), dtype=float)
        
    #pourcentages de reussites dans les differents mouvements, en nb de mouvements
    p = np.sum(r, 0)/np.array(np.sum(compteur, 0), dtype=float)
        
    #enfin la precision par mouvement et en nb de liens
    precision =np.array(np.sum(r, 0)[:,0], dtype=float)/(np.sum(r, 0)[:,0]+np.sum(np.sum(r2, 0),0))
    
    return [r, r2, compteur], [recall, p, precision]

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
    for y in range(len(r)):
        precision[y]=np.array(r[y].transpose()[1], dtype=float)/(r[y].transpose()[1]+np.sum(r2[y], 1))

    return [r, recall, r2, p, precision, compteur]
            
def calcAccuracy(new_sol, newFrameLot, tdict=None):
    stories = tracking.ordonnancement(new_sol) 
    acc = np.zeros(shape=(len(EVENTS), 2), dtype=np.int)
    acc2 = np.zeros(shape=(len(EVENTS), len(EVENTS)), dtype=np.int)
    compteur = [[0,0] for x in range(len(EVENTS))]
    
    visu = {}
    
    
    if tdict == None:
        for plate in newFrameLot.lstFrames:
            if plate not in visu:
                visu[plate]={}
            for well in newFrameLot.lstFrames[plate]:
                print plate, well, 'calculating accuracy on the whole video'
                if well not in visu[plate]:
                    visu[plate][well]={}
                for index in newFrameLot.lstFrames[plate][well]:
                    if index not in visu[plate][well]:
                        visu[plate][well][index]=[]
                    if index+1 not in visu[plate][well]:
                        visu[plate][well][index+1]=[]
                        
                    visuC=visu[plate][well][index] ; nextVisu=visu[plate][well][index+1]
                    #print index,
                    sourceC, targetC = filtre(newFrameLot.lstFrames[plate][well][index])
    #                for b in range(len(sourceC)):
    #                    print "training set :", sourceC[b], targetC[b]
                    if index not in stories[plate][well]:
                        continue
                    solC = stories[plate][well][index]
                    if type(solC.truth)==str:
                        print "attention, il n'y a pas de prediction pour cet index"
                        continue
                   
                    predicted =formulation(filter(lambda x: x[1]==1, zip(solC.hypotheses, solC.truth)))
    #                for k in range(len(predicted[0])):
    #                    print "predictions :", predicted[0][k], predicted[1][k]
                    i=0
                    for uplet in sourceC:
                        t, ttarg = movement(uplet, s = None, tC=targetC[i])
                        #print "truth", uplet, t, ttarg
                        
                        #compteur[0]=nb total de liens, compteur[1]=nb total de matchings 
                        compteur[t][0]+=max(len(uplet), len(ttarg)); compteur[t][1]+=1 
                        first = True
                        for u in uplet:
                            p, ptarg = movement(u, s=predicted, tC = None)
                        #   print "predit pour", u, " : ", p, ptarg
                            #pdb.set_trace()
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
                                    if u!=-1:
                                        visuC.append(filter(lambda x: x.label==u, solC.singlets)[0].center)
                                    else:
                                        nextVisu.append(filter(lambda x: x.label==ttarg, solC.nextSinglets)[0].center)
                                    if t==1:
                                        for tt in ttarg:
                                            z, ztarg = movement (u, s=predicted, tC = tt, incomplet = True)
                                            acc2[t][z]+=1
                                    else:
                                        for tt in ttarg:
                                            if tt in ptarg:
                                                acc[t][0]+=1
                                            else:
                                                z, ztarg = movement(u, s=predicted, tC=tt, incomplet=True)
                                                acc2[t][z]+=1
                            else:
                                #print '-------------------------------------------------------pas pareil...'
                                if u!=-1:
                                    visuC.append(filter(lambda x: x.label==u, solC.singlets)[0].center)
                                else:
                                    nextVisu.append(filter(lambda x: x.label==ttarg, solC.nextSinglets)[0].center)
                                if len(ttarg)>1:
                                    if first:
                                        for tt in ttarg:
                                            z, ztarg = movement (u, s=predicted, tC = tt, incomplet = True)
                                            acc2[t][z]+=1
                                    first = False
                                else:
                                    acc2[t][p]+=1
                        i+=1
    else:
        for plate in newFrameLot.lstFrames:
            for well in tdict[plate]:
                for index in tdict[plate][well]:
                    print plate, well, index
                    sourceC, targetC = filtre(newFrameLot.lstFrames[plate][well][index])

    #                for b in range(len(sourceC)):
    #                    print "training set :", sourceC[b], targetC[b]
                    if index not in stories[plate][well]:
                        continue
                    solC = stories[plate][well][index]
                    if type(solC.truth)==str:
                        print "attention, il n'y a pas de prediction pour cet index"
                        continue
                   
                    predicted =formulation(filter(lambda x: x[1]==1, zip(solC.hypotheses, solC.truth)))
                    
    #                for k in range(len(predicted[0])):
    #                    print "predictions :", predicted[0][k], predicted[1][k]
                    i=0
                    for uplet in sourceC:
                        t, ttarg = movement(uplet, s = None, tC=targetC[i])
                        #print "truth", uplet, t, ttarg
                        
                        #compteur[0]=nb total de liens, compteur[1]=nb total de matchings 
                        compteur[t][0]+=max(len(uplet), len(ttarg)); compteur[t][1]+=1 
                        first = True
                        for u in uplet:
                            p, ptarg = movement(u, s=predicted, tC = None)
                         #   print "predit", p, ptarg
                            #pdb.set_trace()
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
                                        for tt in ttarg:
                                            z, ztarg = movement (u, s=predicted, tC = tt, incomplet = True)
                                            acc2[t][z]+=1
                                    else:
                                        for tt in ttarg:
                                            if tt in ptarg:
                                                acc[t][0]+=1
                                            else:
                                                z, ztarg = movement(u, s=predicted, tC=tt, incomplet=True)
                                                acc2[t][z]+=1
                            else:
                                #print '-------------------------------------------------------pas pareil...'
                                if len(ttarg)>1:
                                    if first:
                                        for tt in ttarg:
                                            z, ztarg = movement (u, s=predicted, tC = tt, incomplet = True)
                                            acc2[t][z]+=1
                                    first = False
                                else:
                                    acc2[t][p]+=1
                                
                        i+=1
    return acc, acc2, np.array(compteur), visu

def ecrireXml(visu, p):
    for plate in visu:
        for well in visu[plate]:
            if p==None:
                fichier = open("PL"+plate+"___P"+well+"___T00000.xml", 'w')
            else:
                fichier = open(str(p)+"PL"+plate+"___P"+well+"___T00000.xml", 'w')
                
            fichier.write(fHacktrack2.initXml())
            compteur=1
            for index in visu[plate][well]:
                for c in visu[plate][well][index]:
                    result = "\n      <Marker_Type>\n        <Type>{0}</Type>".format(compteur)
                    numFramePourXml = index+1
                    result+="\n        <Marker>\n          <MarkerX>{1}</MarkerX>\n          <MarkerY>{2}</MarkerY>\n          <MarkerZ>{0}</MarkerZ>\n        </Marker>".format(numFramePourXml, c[0], c[1])
                    result+="\n      </Marker_Type>"
                    fichier.write(result)
            fichier.write(fHacktrack2.finirXml())
    return

def accuBig(nom, c_list, nb_fold, n_big_fold, bestD = None):
    global nb_folds
    global first
    
    #ACCURACY POUR TOUT LE MONDE avec cross-validation
    if nb_fold == -1 and bestD is not None:
        output2=open(nom+'TOTAL.pkl', 'w')
        d_e={}; d_c={}
        listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
        r = np.zeros(shape=(n_big_fold,len(EVENTS), 2), dtype=np.int)
        r2 = np.zeros(shape=(n_big_fold, len(EVENTS), len(EVENTS)), dtype=np.int)
        compteur = np.zeros(shape=(n_big_fold,len(EVENTS),2))
        for plate in listD:
            sol, newFrameLot = gettingSolu(plate)
            first = False; firstFold=True

            for n_b_f in range(n_big_fold):
                print "big fold ", n_b_f, "and c = 10**", bestD[n_b_f]
                fichier = open(os.path.join(loadingFolder, "data_TEST_fold_BIG"+str(n_b_f)+".pkl"), 'r')
                tdict = pickle.load(fichier); fichier.close()
                if plate not in tdict:
                    continue
    
                if not firstFold:
                    sol.denormalisation(minMax)
    
                print "normalization"
    
                fichier = open(os.path.join(loadingFolder,"minMax"+str(n_b_f)+".pkl"), "r")  
                minMax = pickle.load(fichier)
                fichier.close()
                sol.normalisation(minMax)
                print "TIME TIME TIME after normalization", time.clock()
                firstFold=False
                
                acc, acc2, compt = accuPerPlate(nom, [bestD[n_b_f]], n_b_f, None, tdict, plate, sol, newFrameLot)
                r[n_b_f]+=acc[0]; r2[n_b_f]+=acc2[0]; compteur[n_b_f]+=compt[0]
                    
            outputPlate = open(nom+str(plate)+'.pkl', 'w')
            pickle.dump([r, r2, compteur], outputPlate)
            outputPlate.close()

        pickle.dump(finalPostProcessing(r, r2, compteur), output2)
        output2.close()
        return            
    elif nb_fold == -1 and bestD is None:
#        first = True
        print "______________________________First", first
        listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
        d_e={}; d_c={}
        output2=open(nom+'TOTAL.pkl', 'w')
        for plate in listD:
            sol, newFrameLot = gettingSolu(plate)
            first = False; firstFold=True
            
            for n_b_fold in range(n_big_fold):    
                r = np.zeros(shape=(len(c_list),len(EVENTS), 2), dtype=np.int)
                r2 = np.zeros(shape=(len(c_list), len(EVENTS), len(EVENTS)), dtype=np.int)
                compteur = np.zeros(shape=(len(c_list),len(EVENTS),2))
                
                print "big fold ", n_b_fold, "and c_list", c_list
                fichier = open(os.path.join(loadingFolder, "data_TEST_fold_BIG"+str(n_b_fold)+".pkl"), 'r')
                tdict = pickle.load(fichier); fichier.close()
                if plate not in tdict:
                    continue
    
                if not firstFold:
                    sol.denormalisation(minMax)
    
                print "normalization"
    
                fichier = open(os.path.join(loadingFolder,"minMax"+str(n_b_fold)+".pkl"), "r")  
                minMax = pickle.load(fichier)
                fichier.close()
                sol.normalisation(minMax)
                print "TIME TIME TIME after normalization", time.clock()
                firstFold=False
                
                acc, acc2, compt = accuPerPlate(nom, c_list, n_b_fold, None, tdict, plate, sol, newFrameLot)
                r+=acc; r2+=acc2; compteur+=compt
                if n_b_fold not in d_e: 
                    d_e[n_b_fold]=[r, r2]
                else:
                    d_e[n_b_fold][0]+=r
                    d_e[n_b_fold][1]+=r2 
                if n_b_fold not in d_c: 
                    d_c[n_b_fold]=compteur[0]
                else:
                    d_c[n_b_fold]+=compteur[0]
            
        pickle.dump([d_e, d_c], output2)
        output2.close()
        return d_e, d_c
        #ACCURACY POUR pheno seulement, sans cross-validation
#    elif nb_fold == -1 and bestD==None:
#        output2=open(nom+'TOTAL.pkl', 'w')
#
#        listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
#        r = np.zeros(shape=(1, len(EVENTS), 2), dtype=np.int)
#        r2 = np.zeros(shape=(1, len(EVENTS), len(EVENTS)), dtype=np.int)
#        compteur = np.zeros(shape=(1,len(EVENTS),2))
#        global ctrl; global pheno
#        for plate in listD:
#            sol, newFrameLot = gettingSolu(plate, pheno_only=True)
#            if sol==None:
#                continue
#            else:
#                first = False
#    
#            print "normalization"
#    
#            fichier = open(loadingFolder+"minMax_data_all.pkl", "r")  
#            minMax = pickle.load(fichier)
#            fichier.close()
#            sol.normalisation(minMax)
#            print "TIME TIME TIME after normalization", time.clock()
#                
#            acc, acc2, compt = accuPerPlate(nom, c_list, None, None, None, plate, sol, newFrameLot)
#            r[0]+=acc[0]; r2[0]+=acc2[0]; compteur[0]+=compt[0]
#                    
#            outputPlate = open(nom+str(plate)+'.pkl', 'w')
#            pickle.dump([r, r2, compteur], outputPlate)
#            outputPlate.close()
#
#        pickle.dump(finalPostProcessing(r, r2, compteur), output2)
#        output2.close()
#        return 1
    #ACCURACY for small folds
    elif nb_fold == nb_folds:
#        first = True
        print "______________________________First", first
        listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
        d_e={}; d_c={}
        output2=open(nom+'.pkl', 'w')
        for plate in listD:
            sol, newFrameLot = gettingSolu(plate)
            first = False; firstFold=True
            
            for n_b_fold in range(n_big_fold):    
                print "big number", n_b_fold
                
                r = np.zeros(shape=(len(c_list),len(EVENTS), 2), dtype=np.int)
                r2 = np.zeros(shape=(len(c_list), len(EVENTS), len(EVENTS)), dtype=np.int)
                compteur = np.zeros(shape=(len(c_list),len(EVENTS),2))
                
                for n_fold in range(nb_fold):
                    print "small number", n_fold
                    fichier = open("../results/data_TEST_fold_SMALL"+str(n_b_fold)+'_'+str(n_fold)+".pkl", 'r')
                    tdict = pickle.load(fichier); fichier.close()
                    if plate not in tdict:
                        continue
    
                    if not firstFold:
                        sol.denormalisation(minMax)
    
                    print "normalization"
        
                    fichier = open(loadingFolder+"minMax"+str(n_b_fold)+'_'+str(n_fold)+".pkl", "r")  
                    minMax = pickle.load(fichier)
                    fichier.close()
                    sol.normalisation(minMax)
    
                    firstFold=False
    
                    print "TIME TIME TIME after normalization", time.clock()
                    
                    
                    acc, acc2, compt = accuPerPlate(nom, c_list, n_b_fold, n_fold, tdict, plate, sol, newFrameLot)
                    r+=acc; r2+=acc2; compteur+=compt
                if n_b_fold not in d_e: 
                    d_e[n_b_fold]=[r, r2]
                else:
                    d_e[n_b_fold][0]+=r
                    d_e[n_b_fold][1]+=r2 
                if n_b_fold not in d_c: 
                    d_c[n_b_fold]=compteur[0]
                else:
                    d_c[n_b_fold]+=compteur[0]
            
        pickle.dump([d_e, d_c], output2)
        output2.close()
        return d_e, d_c


def accu(nom, c_list, nb_fold = -1):
    global nb_folds
    #i: charger les data
    r = [np.zeros(shape=(len(EVENTS), 2), dtype=np.int) for x in range(len(c_list))]
    r2 = [np.zeros(shape=(len(EVENTS), len(EVENTS)), dtype=np.int) for x in range(len(c_list))]
    compteur = [np.zeros(shape=(len(EVENTS),2)) for x in range(len(c_list))]
    #ACCURACY POUR TOUT LE MONDE
    if nb_fold == -1:
        listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
        first = True
        for plate in listD:
            output2=open(nom+plate+'.pkl', 'w')
            sol, newFrameLot = gettingSolu(first, plate)
            
            fichier = open(loadingFolder+"minMax_data_all.pkl", "r")  
            minMax = pickle.load(fichier)
            fichier.close()
            sol.normalisation(minMax)
            
            first = False
            j=0
            for p in c_list:
                #ii. predire et changer la forme de new_sol pour qu'elles soient exploitables
                print loadingFolder
                new_sol = tracking.sousProcessClassify(sol, loadingFolder, integers=True, loadingFile= None, i=p)
                dicTraj =tracking.trajBuilder(new_sol)
                tracking.outputTrajXml(dicTraj, loadingFolder+'/predictions/')

                #iii. construire le doublet source, target pour chaque frame, avec uniquement les hypotheses concernant les individus du training set
                acc, acc2, c, visu = calcAccuracy(new_sol, newFrameLot)
                #pickle.dump(visu, open('visu'+plate+'_'+str(p)+'.pkl', 'w'))
                
                #ecrireXml(visu, p)
                
                print acc, acc2, c
                r[j]+=acc; r2[j]+=acc2; compteur[j] +=c
                j+=1
            sol = None ;newFrameLot = None ; new_sol = None
                
            pickle.dump(postProcessing(r, r2, compteur[0]), output2)
            output2.close()
            
    elif nb_fold == nb_folds:
        first = True
        outputTotal = open(nom+"TOTALE.txt", 'w')
        output2=open(nom+'TOTALE.pkl', 'w')
        listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
        for plate in listD:
            sol, newFrameLot = gettingSolu(first, plate)
            first = False; firstFold=True
            for n_fold in range(nb_fold):
                fichier = open("../results/data_TEST_fold"+str(n_fold)+".pkl", 'r')
                tdict = pickle.load(fichier); fichier.close()
                if plate not in tdict:
                    continue

                if not firstFold:
                    sol.denormalisation(minMax)

                print "normalization"
    
                fichier = open(loadingFolder+"minMax"+str(n_fold)+".pkl", "r")  
                minMax = pickle.load(fichier)
                fichier.close()
                sol.normalisation(minMax)

                firstFold=False

                print "TIME TIME TIME after normalization", time.clock()
                
                
                acc, acc2, compt = accuPerPlate(nom, c_list, n_fold, tdict, plate, sol, newFrameLot)
                for pou in range(len(r)):
                    r[pou]+=acc[pou]; r2[pou]+=acc2[pou]
                    compteur[pou]+=compt[pou]
                    
            outputPlate = open(nom+str(plate)+'.pkl', 'w')
            pickle.dump(postProcessing(r, r2, compteur[0]), outputPlate)
            outputPlate.close()
        print>>outputTotal, "c list :", c_list
        for k in range(len(c_list)):
            print>>outputTotal, r[k], compteur[k]
            #print>>outputTotal, d
        outputTotal.close()
        pickle.dump(postProcessing(r, r2, compteur[0]), output2)
        output2.close()
    return r, r2, compteur

def filtre2(sol, tdict, plate):
    l = sol.lstSolutions
    r = filter(lambda x: x.well in tdict[plate] and x.index in tdict[plate][x.well], l)
    sol2 = joining.Solutions(solution = None, lstSolutions = r)
    return sol2
            
def accuPerPlate(nom, c_list, n_big_fold, n_fold, tdict, plate, sol, newFrameLot): 
    r = np.zeros(shape=(len(c_list),len(EVENTS), 2), dtype=np.int)
    r2 = np.zeros(shape=(len(c_list), len(EVENTS), len(EVENTS)), dtype=np.int)
    compteur = np.zeros(shape=(len(c_list),len(EVENTS),2))
    j=0

    for p in c_list:
        #ii. predire et changer la forme de new_sol pour qu'elles soient exploitables
        if tdict is not None:
            sol = filtre2(sol, tdict, plate)
        new_sol = tracking.sousProcessClassify(sol, loadingFolder, loadingFile= None, i=p, n_f=n_fold, n_big_f=n_big_fold)

        #iii. construire le doublet source, target pour chaque frame, avec uniquement les hypotheses concernant les individus du training set
        acc, acc2, compt, visu = calcAccuracy(new_sol, newFrameLot, tdict)
        r[j]+=acc ; r2[j]+=acc2; compteur[j] +=compt
        j+=1
    del sol ;del newFrameLot ; del new_sol ; del tdict; gc.collect()
            
#    print>>outputTotal, "n fold", n_fold, "plate", plate, "c list :", c_list
#    for k in range(len(c_list)):
#        print>>outputTotal, r[k] 
#    #print>>outputTotal, d
#    outputTotal.close()
        
    return r, r2, compteur

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Accuracy")
    parser.add_argument('-s', '--step', default=0, type=int, dest ='step')#possible values : 0 for internal accuracy 1 for

    parser.add_argument('-n', '--name', default='accuracy', type = str, dest='nom')
    parser.add_argument('-lf', '--loading-folder', default='../results/', type = str, dest='loadingFolder')
     
    args = parser.parse_args()

    step = args.step
    nom = args.nom
    loadingFolder = args.loadingFolder
    #nom="Naccuracy_CVinterne_nomorpho"
    global first
    first= True
    print 'Here we measure accuracy'
    print 'c_list', c_list
    #print>>outputTotal,"TIME TIME TIME TIME TIME", time.clock()
    output = []
#CROSS-VALIDATION calculating accuracy on small folds to choose model
    if step==0:
        f=open(nom+'bestModel.pkl', 'w')
        d_e, d_c =accuBig(nom, c_list, nb_folds, nb_big_folds)
        for n_b_f in d_e:
            if n_b_f not in d_c:
                raise
            
            compteur = d_c[n_b_f]
            rights = d_e[n_b_f][0]
            recall = [0 for z in range(len(rights))]
            for z in range(len(recall)):
                recall[z]=np.sum(rights[z], 0)/np.array(np.sum(compteur, 0), dtype=float)
            output.append(recall)
        pickle.dump(np.array(output), f)
#def meana(fold, c):
#    return np.mean(d_e[fold][0][c][:,1]/d_c[fold][:,1])
#def meanp(fold, c):
#    return np.mean(d_e[fold][0][c][:,0]/np.array(d_e[fold][0][c][:,0]+np.sum(d_e[fold][1][c], 0), dtype=float))


#CROSS VALIDATION once best model was chosen on this big fold, we train it
    elif step==1:
        b={0:None, 1:1.5, 2:1.5, 3:1.5, 4: 1.5}
        f=open(nom+'bestModel.pkl', 'r')
        d= pickle.load(f)
        
        for n_b_f in range(d.shape[0]):
            b[n_b_f]=c_list[np.argmax(d[n_b_f].transpose()[1])]
        print b; pdb.set_trace()
        Parallel(n_jobs=3, verbose=5)(delayed(test.wholeProcess)(b[n_big_fold], n_big_fold, -1, True, loadingFolder) for n_big_fold in b)
#CROSS-VALIDATION final step, calculating accuracy on each new model for each big fold
    elif step==2:
        b={0:2, 1:2, 2:2, 3:2, 4:2}
        print "loading folder", loadingFolder
        print 'step', step
        f=open(nom+'bestModel.pkl', 'r')
        d= pickle.load(f); f.close()
        
        for n_b_f in range(d.shape[0]):
            b[n_b_f]=c_list[np.argmax(d[n_b_f].transpose()[1])]
        print b
#        
#        f=open(nom+'bestModel_2.pkl', 'w')
        d_e, d_c = accuBig(nom, c_list, -1, nb_big_folds, b)
        #pickle.dump([d_e, d_c], f)
#ONLY FOR MODEL CHOICE here we do simple cross-validation to choose a model on the whole training set
    elif step==3:
        print "loading folder", loadingFolder
        print 'step', step
        accuBig(nom, c_list, -1, nb_big_folds, None)
    else:
        print "je n'ai pas compris"
    try:
        f.close()
    except NameError:
        pass
    #r = [np.zeros(shape=(len(EVENTS), 2), dtype=np.int) for x in range(len(c_list))]
    #r2 = [np.zeros(shape=(len(EVENTS), len(EVENTS)), dtype=np.int) for x in range(len(c_list))]
#    print n_big_fold
#    recall = [0 for z in range(len(acc))]
#    for z in range(len(acc)):
#        recall[z]=np.sum(acc[z], 0)/np.array(np.sum(compteur, 0), dtype=float)
#        
#    r=[]
#    for tt in recall:r.append(tt[1])
#    index_max = r.index(max(r)); output[index_max]+=1
#    bestC=c_list[index_max]
#    print "best model for this fold: 10**", bestC 
#    nomC=nomC+str(bestC)

#    for n_big_fold in range(nb_big_folds):
#        test.wholeProcess(1.5, n_big_fold, -1, True, loadingFolder)
    
        
#    for n_big_fold in range(nb_big_folds):
#        n_big_fold+=4
#        nomC=nom+str(n_big_fold)
#        acc, acc2, compteur =accuBig(nomC, c_list, nb_folds, n_big_fold)
#        #r = [np.zeros(shape=(len(EVENTS), 2), dtype=np.int) for x in range(len(c_list))]
#        #r2 = [np.zeros(shape=(len(EVENTS), len(EVENTS)), dtype=np.int) for x in range(len(c_list))]
#        print n_big_fold
#        recall = [0 for z in range(len(acc))]
#        for z in range(len(acc)):
#            recall[z]=np.sum(acc[z], 0)/np.array(np.sum(compteur, 0), dtype=float)
#            
#        r=[]
#        for tt in recall:r.append(tt[1])
#        index_max = r.index(max(r)); output[index_max]+=1
#        bestC=c_list[index_max]
#        print "best model for this fold: 10**", bestC 
#        nomC=nomC+str(bestC)
#        
#        test.wholeProcess(bestC, n_big_fold, -1, True, loadingFolder)
#       # pdb.set_trace()
#        acc,acc2,compteur = accuBig(nomC, [bestC], -1, n_big_fold)        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
