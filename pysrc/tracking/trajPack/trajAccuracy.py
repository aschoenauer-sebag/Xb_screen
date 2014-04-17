import os, pdb, time, sys
from math import sqrt
import numpy as np
np.seterr(all='warn')
from PyPack import fHacktrack2
from importPack import imp, EVENTS, FEATURE_NUMBER, nb_folds, c_list
from trajPack import moments, featuresList, basic_colors, d_ctrl
from dataPack import treatments, joining, classify
import cPickle as pickle
import test, tracking
from trajPack import trajFeatures
from scipy.stats import linregress
import matplotlib.pyplot as p
from mpl_toolkits.mplot3d import Axes3D

count = 0
monOutput = ""
loadingFolder = "../results/"
movie_length={}
#donc cette fonction va recuperer le training set VRAI ie incl merge divides, elle va charger un modele et mettre qqpart la prediction du modele,
# et elle va retourner le tout pour faire un calcul d'accuracy
# trajAccuracy.gettingSolu


#accuracy.movement, formulation, filtre, postProcessing, filtre2, accuPerPlate

def trajBuilder(new_sol, training=False):
    
    stories = tracking.ordonnancement(new_sol); res = {}
    global movie_length
    
    for plate in stories:
        res.update({plate : {}})
        movie_length.update({plate : {}})
        for well in stories[plate]:
            print plate, well
            ensembleTraj = fHacktrack2.ensTraj()
            res[plate].update({well:ensembleTraj})
            movie_length[plate].update({well : max(stories[plate][well].keys())+1})
#            print ensembleTraj.numMitoses
#            print ensembleTraj.lstTraj
            num = 0; l=[0]; mergeList=-1; trajC=None
            for index in stories[plate][well]:
                print " INDEX ", index,
                solC = stories[plate][well][index]
                if type(solC.truth)==str:
                    l.append(index+1)
                    continue
                singlets = solC.singlets
                nextSinglets = solC.nextSinglets
                result =filter(lambda x: x[1]==1, zip(solC.hypotheses, solC.truth))
                
                for hyp in result:
                    #print hyp
                    #if index==86: pdb.set_trace()
                    s, t = hyp[0]; source = [] ; target = []; 
                    if type(s) in [np.int32, np.uint16] and s!=-1:
                        source.append(filter(lambda x: x.label == s, singlets)[0])
                    elif type(s)==tuple:        
                        for c in s:
                            source.append(filter(lambda x: x.label == c, singlets)[0])
                    if type(t) in [np.int32, np.uint16] and t!=-1:
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
                                        traj.split(index); traj.numCellule = numC; first = False
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
                                    if first: 
                                        numC = incomingTraj[i].numCellule
                                        trajC = incomingTraj[i]  
                                        mergeList = incomingTraj[i].fusion                   
                                    if mit: 
#                                        print "ajout mitose"
                                        incomingTraj[i].split(index); incomingTraj[i].numCellule = numC; first = False
                                        #incomingTraj[i].ajoutMomentInteressant(c.center[0], c.center[1], index)
                                        ensembleTraj.numMitoses+=1
                                    i+=1
                                except IndexError:
                                    if not mit: 
                                        #ici on y arrive a t=0 quand aucune trajectoire n'est initialise ou cas de split a quatre sorti du TS
                                        if index>0 and training==False: pdb.set_trace()
                                        num+=1
                                        newTraj = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label)
                                        newTraj.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
                                        newTraj.longueurs=[index]; newTraj.fractions = [1]
                                    if first: 
                                        numC = num    
                                        if index>0 and training==False:
                                            print hyp
                                            pdb.set_trace()
                                    if mit:
                                        #ici, quand j'ai une mitose je copie toute la trajectoire precedente dans la nouvelle OU PAS ?
                                        newTraj = fHacktrack2.trajectoire(numC, c.center[0], c.center[1], index, c.label)
                                        newTraj.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
                                        newTraj.longueurs = incomingTraj[0].longueurs if incomingTraj !=[] else [index]
                                        newTraj.fractions = incomingTraj[0].fractions if incomingTraj !=[] else [1]
                                        #print id(trajC), numC, hyp, mergeList
                                        #pdb.set_trace()
                                        newTraj.fusion=list(mergeList) if type(mergeList)==list else mergeList
                                        #pdb.set_trace()
#                                        print "ajout mitose"
                             
                                        newTraj.split(index)
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
                        if index>0:
                            tt=[]; ll=[]; lf=[]
                            for c in source:
                                incomingTraj = ensembleTraj.findLabelFrame(c.label, index)
                                try:
                                    tt.append(incomingTraj[0])
                                except IndexError:
                                    pass
                                else:
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
                                #on tombe ici a t=0 mais est-ce qu'on y tombe apres ? oui si on calcule les trajectoires du training
                                #set pcq il y a des evenements elimines par defaut (mvt a quatre)
                                if index>0 and training==False: pdb.set_trace()
                                num+=1
                                trajC = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label)
                                trajC.longueurs=[index for bb in range(incoming)]; trajC.fractions = [1/float(incoming) for bb in range(incoming)]
                                ensembleTraj.add(trajC)
                            finally:    
                                trajC.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
                                #print id(trajC), trajC.numCellule, trajC.fusion
                                #pdb.set_trace()
                                trajC.merge(index)
                            i+=1
                            #incomingTraj[0].ajoutMomentInteressant(c.center[0], c.center[1], index)
    return res

def trajPreparation(dicT, training):
    r={}; coord={}; nuCount=[0 for x in range(len(moments))]
    for plate in dicT:
        if plate not in r:
            r[plate]={}
            coord[plate]={}
        for well in dicT[plate]:
            if well not in r[plate]:
                r[plate][well]=[]
                coord[plate][well]=[]
            print plate, well
            dicC = dicT[plate][well]
            fichier=open("PL{0}___P{1}___T00000.xml{2}".format(plate, well, training), 'w'); fichier.write(fHacktrack2.initXml())

            tracklets=[]
            compteur = 1
            #signature de featuredTraj.__init__ : (self, ide, labelsSequence, beginning, end, num, features, to=None, fr=None):
            #going into all trajectories to get small tracklets
            #AND keeping in memory which tracklet they're related to
            
            #once we have all tracklets we can do the findCorresp
            
            #FINALLY we can recoller les morceaux and get model accuracy/trajectory instead of /movement from t to t+1
            print 'Pass of all tracks'
            for traj in dicC.lstTraj:                
                t=traj.lstPoints.keys(); t.sort(); 
                
#                if traj.numCellule ==36:
#                    fichier.write(fHacktrack2.ecrireXml(36, traj))
                    
                print "une trajectoire", traj.numCellule
                print "merge :", traj.fusion
                print "split :", traj.mitose
                print 'debut de la trajectoire :', min(t), "fin ", max(t)
#CHECK qu'est-ce que min(t), max(t)
#PLUS compute to and from for all tracklets
                #Here we don't want to delete all that happens after a merge
                #we want to cut the trajectory into pieces between all split and merge
                events = list(traj.fusion); events.extend(traj.mitose); events.sort(reverse=True)
                events = filter(lambda x: x>min(t)[0] and x<max(t)[0], events)[0]
#                if traj.fusion !=-1:
#                    if min(traj.fusion)<min(t):
#                        continue
#                    traj=traj.copyBefore(min(traj.fusion))
#                    t=sorted(traj.lstPoints.keys())
#                        
#                print 'apres le check merge, plus de merge normalement : ', traj.fusion , 'et que des splits avant : ', traj.mitose
                print "evenements de la trajectoire ", events
                if events !=[]:
                    end = max(t)[0]
                    for event in events:
                    #    if end-mit>9:
                        print "event ", event
                        track = traj.copyBetween(event, end)
                        r[plate][well].append(trajFeatures.featuredTraj(compteur, labelsSeq(track), track.numCellule, None, to, fr))
                        compteur +=1
                        #tracklets.append()
                        end = event
                    
                    track= traj.copyBefore(end)
                    r[plate][well].append(trajFeatures.featuredTraj(compteur, labelsSeq(track), track.numCellule, None, to, fr))
                    compteur +=1
                    #tracklets.append()
                else:#if end-t[0][0]>9:
                    track = traj
                    r[plate][well].append(trajFeatures.featuredTraj(compteur, labelsSeq(track), track.numCellule, None, to, fr))
                    compteur +=1
                    #tracklets.append(traj)    
                        

#ATTENTION here we throw away tracklets whose length is smaller than 10 (otherwise we don't have diffusion info because
#I use delta t in range (1, len(track)/3 +1) ie (1,2,3) if len(track)>9, and then do regression on it,
#so if we have less than three points to do linear regression it seems a bit dubious to me)    
#            if not len(tracklets)==len(filter(lambda x: len(x.lstPoints)>9 ,tracklets)):
#                raise             
            print 'Second pass of all tracks to extract features of tracks longer than 9 frames'
            #euh no here we're not doing this
#            pbl=dict(zip(moments, [[] for x in range(len(moments))]))
#            for track in filter(lambda x: len(x.lstPoints)>9 ,tracklets):
#                t=track.lstPoints.keys(); t.sort(); 
#                labelsSequence = [k[1] for k in t]
#                x= [track.lstPoints[k][0] for k in t]
#                y =[track.lstPoints[k][1] for k in t]
#                #pdb.set_trace()
#                coord[plate][well].append([x, y])
#                print "debut de l'extrait de trajectoire", min(t)[0], "fin ", max(t)[0]
#                
#    #Pour l'instant on ne normalise pas par rapport au deplacement moyen sur toute l'image car on ne considere qu'une partie de l'image
#                f=computingFeatures(track, movie_length[plate][well], np.mean(average, 0))
#                if 'pbl' in f:
#                    nu=f['pbl']
#                    #nuCount[nu].append(len(track.lstPoints))
#                    pbl[nu].append(track)
#                    fichier.write(fHacktrack2.ecrireXml(nu, track))
#                r[plate][well].append(featuredTraj(compteur, labelsSequence, min(t)[0], max(t)[0], track.numCellule, f))
#                if track.numCellule ==12:
#                    fichier.write(fHacktrack2.ecrireXml(12, traj))    
#                compteur +=1
##            
#            fichier.write(fHacktrack2.finirXml()); fichier.close()  
#    f=open('pblCoeffCorr{0}.pkl'.format(training), 'w')
#    pickle.dump(pbl, f); f.close() 
    return r, coord

def labelsSeq (track):
    t=track.lstPoints.keys(); t.sort(); 
    labelsSequence = [k[1] for k in t]
    return labelsSequence, min(t)[0], max(t)[0]

if __name__ == '__main__':
    training=True
    predict = True
    TroisD = False
    d, td=trajFeatures.output(training, predict); 
    print movie_length
    #FOR TRAINING SET DATA
    tFeatures, tCoordonnees = trajPreparation(td, True)
    trajFeatures.plot(tFeatures, tCoordonnees, TroisD)

    #FOR PREDICTED DATA
    features, coordonnees = trajPreparation(d, False)
    
#    plot(features, coordonnees, TroisD)
    
    #here I consider every predicted tracklet
    correspTF={}; correspF={}; N=0
    for plate in features:
        if plate not in correspF:
            correspF[plate]={}
            correspTF[plate]={}
            
        for well in features[plate]:
            if well not in correspF[plate]:
                correspF[plate][well]=[]
                correspTF[plate][well]=[]
                
            fC = features[plate][well]; 
            tfC = tFeatures[plate][well]

#from the point of view of predicted data
##first i see which predicted tracklet is present when
#            r=WhoIsWhen(tFeatures, plate, well)
#            corresp = {}
#            toPlotF = []; toPlotCoord=[]
#            k=0
#            for tracklet in fC:
#                rC = r[tracklet.beginning: tracklet.end]; zz=[]
#                for l in rC: zz.extend(l)
#                h=np.bincount(zz)
#                
#                correspondance, nbCommonFrames, feat=tracklet.findCorresp(filter(lambda x: x.id in np.where(h==max(h))[0], tfC)) #filter renvoie un tuple avec une liste comme seul element
#                #print correspondance, nbCommonFrames, tracklet.end-tracklet.beginning+1
#                
#                if correspondance is not None:
#                    corresp.update({tracklet.id:correspondance})
#                    correspTF[plate][well].append(feat)
#                    toPlotF.append(fC[k]); toPlotCoord.append(coordonnees[plate][well][k])
#                    #pdb.set_trace()
#                k+=1
#            correspF[plate][well]=list(toPlotF)
#            N+=1

#from the point of view of training data
#first i see which predicted tracklet is present when
            r=trajFeatures.WhoIsWhen(features, plate, well)
            corresp = {}
            toPlotF = []; toPlotCoord=[]
            k=0
            for tracklet in tfC:
                rC = r[tracklet.beginning: tracklet.end]; zz=[]
                for l in rC: zz.extend(l)
                h=np.bincount(zz)
                
                correspondance, nbCommonFrames, feat=tracklet.findCorresp(filter(lambda x: x.id in np.where(h==max(h))[0], fC)) #filter renvoie un tuple avec une liste comme seul element
                #print correspondance, nbCommonFrames, tracklet.end-tracklet.beginning+1
                
                if correspondance is not None:
                    corresp.update({tracklet.id:correspondance})
                    correspTF[plate][well].append(feat)
                    toPlotF.append(fC[k]); toPlotCoord.append(coordonnees[plate][well][k])
                else:
                    pdb.set_trace()
                    #pdb.set_trace()
                k+=1
            correspF[plate][well]=list(toPlotF)
            N+=1

            #pdb.set_trace()
            plot(toPlotF, toPlotCoord, TroisD, training=False, plate=plate, well=well)
    calculus = plotComparison(correspF,correspTF, N)
#    ouput=open("plotComparOutput.pkl", 'w')
#    pickle.dump(calculus, output)
#    output.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
