import sys, os, os.path, random, csv, pdb
import vigra as v
import vigra.impex as vi
import numpy as n
from tracking.PyPack import dossierSauvegarde as ds
from fHacktrack2 import ensTraj
from fHacktrack2 import trajectoire
from fHacktrack2 import coord
import fHacktrack2

LENGTH=None
EVENTS = []
EVENTS.append("move")
EVENTS.append("appear")
EVENTS.append("disappear")
EVENTS.append("merge")
EVENTS.append("split")

#POUR L'INSTANT ON TRAVAILLE EN DISANT QUE LES SPLIT ET MERGE SONT MAX A CINQ CELLULES

EVENTSSIZE = {}
EVENTSSIZE["move"]=[1,1]
EVENTSSIZE["appear"]=[1,1]
EVENTSSIZE["disappear"]=[1,1]
EVENTSSIZE["merge"]=[5,1]
EVENTSSIZE["split"]=[1,5]

#This function uses segmentation information from the hdf5 files to retrieve cell labels from coordinates in the xml file. They therefore
#translate the annotations files into trajectories thereafter used to produce annotations data files usable by the algo.
#a new trajectory at time T has a point at time T-1 with label -1, and at time T+end of trajectory+1 a point with label -2. This enables the
#following function, ecrireTraining, to find appearing and disappearing cells

#f=open('d.pkl', 'w')

class SameLabelException(Exception):
    def __init__(self):
        pass
    
    def __str__(self):
        return "attention same labels in merge or split event"

def lireTrajectoires(nomFichier, segmentationTab):
    f=open(nomFichier[-36:], 'w'); toCheck={}
#    (X, Y, upX, upY, dX, dY) = coord(nomTxt)
#    frameLabels = X.keys()
    try:
        fichierX =  open(nomFichier, "r")
    except IOError:
        print nomFichier, " n'existe pas"
        return None
    listK = []
    ensembleTraj = ensTraj()
    traj = trajectoire(0,0,0,0,0)
    k=0
    continuer = True
    mesLignes = fichierX.readlines()
    mesLignes.reverse()

    while continuer:
        ligneCourante = mesLignes.pop()
        if "<Type>" in ligneCourante:
            if not ensembleTraj.lstTraj==[]:
                traj.ajoutPoint(0, 0, frame+1, -2)

            ensembleTraj.ajouter(traj)
            k = int(ligneCourante[14:-8])
            listK.append(k)
            nouveau = True

        elif "<MarkerX>" in ligneCourante:
            xi = int(ligneCourante[19:-11])
            ligne2 = mesLignes.pop()
            yi = int(ligne2[19:-11])
            ligne3 = mesLignes.pop()
            frame = int(ligne3[19:-11])-1
            #print xi, yi, frame
            
            #ici, plutot charger le fichier hdf5 avec la segmentation et regarder a quel label la coordonnee correspond dans le fichier hdf5
            #donc on va utiliser le tableau arraySegmentation qui contient toute la segmentation de l'image

            segmentationF = segmentationTab[0][frame][0][:][:]
            l=segmentationF[yi][xi]
            if int(l)==0:
                print "ATTENTION CELL IDENTIFICATION PROBLEM ON FRAME ", frame, xi, yi
                if frame not in toCheck:
                    toCheck.update({frame:[]})
                toCheck[frame].append((xi, yi))
                #raise
#            if frame==34:
#                pdb.set_trace()
#            for f in frameLabels:
#                if int(f[0])==int(frame):
#                    upxi =int(upX[f])
#                    dxi=int(dX[f])
#                    upyi = int(upY[f])
#                    dyi=int(dY[f])
#                    if int(dxi)>=int(xi) and int(xi)>=int(upxi) and int(dyi)>=int(yi) and int(yi)>=int(upyi):
#                        l2 = f[1]
#                        #if frame == 91 or frame==92: print "TROUVE "+ str(dxi), upxi, xi, dyi, upyi, yi, l, frame
                        
#            if (l!=l2): print "les deux labels ne sont pas identiques attention"
            
            if nouveau:
              #  print "ajout de la trajectoire "+str(k)
                if int(frame)>0:
                    traj = trajectoire(k, 0, 0, frame-1, -1)
                    traj.ajoutPoint(xi, yi, frame, l)
                else: 
                    traj = trajectoire(k, xi, yi, frame, l)
               # ensembleTraj.ajouter(traj)

            else:
                #traj= ensembleTraj.trouver(k).lstTraj[0]
                traj.ajoutPoint(xi, yi, frame, l)

            nouveau = False

        elif "</Marker_Data>" in ligneCourante:
            traj.ajoutPoint(0, 0, frame+1, -2)
            ensembleTraj.ajouter(traj)
            continuer = False

 #   for traj in ensembleTraj.lstTraj: print traj.lstPoints
    f.write(fHacktrack2.initXml())
    for index in toCheck:
        zz=1
        for chpt in toCheck[index]:
            result = "\n      <Marker_Type>\n        <Type>{0}</Type>".format(zz)
            result+="\n        <Marker>\n          <MarkerX>{1}</MarkerX>\n          <MarkerY>{2}</MarkerY>\n          <MarkerZ>{0}</MarkerZ>\n        </Marker>".format(index+1, chpt[0],chpt[1])
        
            result+="\n      </Marker_Type>"
            f.write(result)
            zz+=1
    f.write (fHacktrack2.finirXml())
            
    fichierX.close(); f.close()
    return ensembleTraj

#This function uses the trajectories which are produced by the above function to write data files, used as input for the tracking algo
def ecrireTraining(ensembleTraj, fileToSave, LENGTH):
    solutions = {}
    fileS=open(fileToSave, "w")
    fileSHDF5 = fileToSave[:-3]+"hdf5"

    going = {}
    coming = {}
    appearCandidates = {}
    disappearCandidates = {}

    for trajectoire in ensembleTraj.lstTraj:
        frameLabels = sorted(trajectoire.lstPoints.iterkeys())
        lastLabels = []
        lastFrame = -1
        lastMoveSplit = [False for x in range(100)]

        for f in frameLabels:
            frame = int(f[0])
            nextFrame = frame+1

            if nextFrame==LENGTH:
                break

            if lastLabels == []:
                lastLabels=trajectoire.findFrame(frame)

            if frame == lastFrame:
                continue

            nextLabels = trajectoire.findFrame(nextFrame)

            if not frame in solutions:
                solutions[frame]={}

            l = len(lastLabels)
            nextL = len(nextLabels)
            #print "FRAME :"+str(frame),l, nextL, lastLabels, nextLabels, trajectoire.numCellule


            if int(l)==int(nextL):
                if int(l)==1 and int(lastLabels[0])==-1:
                    appearing = nextLabels[0]
                    if not nextFrame in appearCandidates:
                        appearCandidates[nextFrame]=[]
                    appearCandidates[nextFrame].append(appearing)

                elif int(l)==1 and int(nextLabels[0])==-2:
                    disappearing = lastLabels[0]
                    if not frame in disappearCandidates:
                        disappearCandidates[frame]=[]
                    disappearCandidates[frame].append(disappearing)
                    #continue

                elif int(l)==1:
                    #print "MOVE"
                                    
                    if not "move" in solutions[frame]:
                        solutions[frame]["move"]=[[], []]

                    solutions[frame]["move"][0].append(lastLabels)
                    solutions[frame]["move"][1].append(nextLabels)
                    if not nextFrame in coming:
                        coming[nextFrame]=[]
                    for label in nextLabels: 
                        coming[nextFrame].append(label)
                        
                    if not frame in going:
                        going[frame]=[]
                    for label in lastLabels:
                        going[frame].append(label)

                else:
                    print "problem", lastLabels, nextLabels, frame
                    raise

            else:
                print "problem at frame "+str(frame), lastLabels, nextLabels
                raise

            lastFrame = frame
            if not lastMoveSplit[frame]:
                 lastLabels = nextLabels
            if int(lastLabels[0])==-2:
                break
            
#looking at merges : if two cells have the same target, it means it is not a move but a merge

    for f in solutions:
        if "move" not in solutions[f]:
            continue
        
        listeLabelsTargetMove = solutions[f]["move"][1]
        listeLabelsSourceMove = solutions[f]["move"][0]
        nextFrame = int(f+1)
        #print listeLabelsSourceMove, listeLabelsTargetMove
        
        newSourceMove = []
        newTargetMove = []
        
        merges = []
        moves = []
        for t_el in listeLabelsTargetMove:
            if listeLabelsTargetMove.count(t_el) > 1:
                merges.append(t_el)
            else:
                moves.append(t_el)

        indexToDel = []
        for label1 in merges:
          #  print label1
            nextLabels=label1
            lastLabels=[]
            for i in range(len(listeLabelsTargetMove)):
                if listeLabelsTargetMove[i]==label1:
                    lastLabels.append(listeLabelsSourceMove[i][0])
                    indexToDel.append(i)
#                    

            if not "merge" in solutions[f]:
                solutions[f]["merge"]=[[], []]
            if nextLabels in solutions[f]["merge"][1]:
                continue
            solutions[f]["merge"][0].append(lastLabels)
            solutions[f]["merge"][1].append(nextLabels)

            if not f in going:
                going[f]=[]
            if not nextFrame in coming:
                coming[nextFrame]=[]
            for label in nextLabels:
                coming[nextFrame].append(label)
            for label in lastLabels:
                going[f].append(label)
        
#looking at splits : if two cells have the same source, it means it is not a move but a split
        
        splits = []
        moves = []
        for s_el in listeLabelsSourceMove:
            if listeLabelsSourceMove.count(s_el) > 1:
                splits.append(s_el)
            else:
                moves.append(s_el)
        
#        indexToDel = []
        for label1 in splits:
            #print label1
            lastLabels=label1
            nextLabels=[]
            for i in range(len(listeLabelsSourceMove)):
                if listeLabelsSourceMove[i]==label1:
                    indexToDel.append(i)
                    nextLabels.append(listeLabelsTargetMove[i][0])
#                    
            if not "split" in solutions[f]:
                solutions[f]["split"]=[[], []]
            if lastLabels in solutions[f]["split"][0]:
                continue
            solutions[f]["split"][0].append(lastLabels)
            solutions[f]["split"][1].append(nextLabels)

            if not f in going:
                going[f]=[]
            if not nextFrame in coming:
                coming[nextFrame]=[]
            for label in nextLabels:
                coming[nextFrame].append(label)
            for label in lastLabels:
                going[f].append(label)
                
        for i in range(len(listeLabelsSourceMove)):
            if i not in indexToDel:
                newSourceMove.append(listeLabelsSourceMove[i])
                newTargetMove.append(listeLabelsTargetMove[i])
 #       print listeLabelsSourceMove, listeLabelsTargetMove
        #print newSourceMove, newTargetMove
        solutions[f]["move"][1] = newTargetMove
        solutions[f]["move"][0] = newSourceMove
        
                                           

#appear candidates 
    for frame in appearCandidates.keys():
        appearlist = appearCandidates[frame]
        for label in appearlist:
            if frame not in coming.keys() or (frame in coming.keys() and label not in coming[frame]):
                print "APPEAR "+str(label)+" on frame "+str(frame)
                if not "appear" in solutions[frame-1]:
                    solutions[frame-1]["appear"]=[[], []]

                solutions[frame-1]["appear"][0].append([-1])
                solutions[frame-1]["appear"][1].append(label)

#disappear candidates
    for frame in disappearCandidates.keys():
        disappearlist = disappearCandidates[frame]
        for label in disappearlist:
            if frame not in going.keys() or (frame in going.keys() and label not in going[frame]):
                print "DISAPPEAR "+str(label)+" on frame "+str(frame)

                if not "disappear" in solutions[frame]:
                    solutions[frame]["disappear"]=[[], []]

                solutions[frame]["disappear"][0].append(label)
                solutions[frame]["disappear"][1].append([])

    count={}
    out_merge = "\n MERGE \n"
    out_split = "\n SPLIT \n"
    for e in EVENTS:
        count[e]=0

    for f in solutions:
        fileS.write("\n --------------------------FRAME "+str(f)+"------------------------------------")
        for e in solutions[f]:
            s = len(solutions[f][e][0])
#attention a l'ordre des axes : ici cx=nombre d'evenements, "taille" de l'evenement marche bien (de toutes facons lorsqu'on appelle la fonction writeHDF5 elle retablit l'ordre numpy)
            shapeS = (s,EVENTSSIZE[e][0],)
            shapeT = (s,EVENTSSIZE[e][1],)

            tabSource = v.VigraArray(shapeS,n.int32, axistags = v.VigraArray.defaultAxistags('cx'), init=True)
            tabTarget = v.VigraArray(shapeT,n.int32, axistags = v.VigraArray.defaultAxistags('cx'), init=True)
            tabSource = tabSource-1
            tabTarget = tabTarget-1
            path = "Training/{:0>6}/{}/".format(f, e)
            pathSource = path+"Source"
            pathTarget = path+"Target"
            
            fileS.write("\n EVENT "+str(e)+"*******\n SOURCES :\n")
            j=0
            for label in solutions[f][e][0]:
                if isinstance(label, int) or isinstance(label, n.uint16) :
                    length = 1
                else:
                    length=int(len(label))
                fileS.write("\n label"+str(label))

                if label == []:
                    continue

                diff = EVENTSSIZE[e][0] - length
                if diff==0:
                    tabSource[j]=label
                else:
                    while diff>0:
                        try:
                            label.append(-1)
                        except AttributeError:
                            print "evenement "+str(e)+" probleme de taille entre evenement"+str(EVENTSSIZE[e][0])+" et le training set "+str(length)
                        else:
                            diff-=1
                            
                    tabSource[j]=label

#                if e == "split":print e, f, label, j, tabSource
                j+=1
                count[e]+=1
            i=0
            fileS.write("\n TARGETS :\n")
            for label in solutions[f][e][1]:
                if isinstance(label, int) or isinstance(label, n.uint16) :
                    length = 1
                else:
                    length=int(len(label))
                fileS.write("\n label"+str(label))
                
                diff = EVENTSSIZE[e][1] - length
                if diff==0:
                    tabTarget[i]=label
                else:
                    while diff>0:
                        try:
                            label.append(-1)
                        except AttributeError:
                            print "evenement "+str(e)+" probleme de taille entre evenement"+str(EVENTSSIZE[e][1])+" et le training set "+str(length)
                        else:
                            diff-=1
                    tabTarget[i]=label

                i+=1

            if i<>j: print "probleme : difference de longueurs entre Source et Target at fr "+str(f)+" pour l'evenement "+e

            vi.writeHDF5(tabSource, fileSHDF5, pathSource)
            vi.writeHDF5(tabTarget, fileSHDF5, pathTarget)
            
            if e =="merge":
                out_merge +="\n"+str(solutions[f][e][0])+" on frame "+str(f)+" to "+str(solutions[f][e][1])
            if e=="split":
                out_split +="\n"+str(solutions[f][e][0])+" on frame "+str(f)+" to "+str(solutions[f][e][1])
    print count, out_merge, out_split
    fileS.close()
    return solutions

def ecrireTrainingSansHDF5(ensembleTraj, fileToSave, LENGTH):
    solutions = {}
    #fileS=open(fileToSave, "w")

    going = {}
    coming = {}
    appearCandidates = {}
    disappearCandidates = {}

    for trajectoire in ensembleTraj.lstTraj:
        frameLabels = sorted(trajectoire.lstPoints.iterkeys())
        lastLabels = []
        lastFrame = -1
        lastMoveSplit = [False for x in range(100)]

        for f in frameLabels:
            frame = int(f[0])
            nextFrame = frame+1

            if nextFrame==LENGTH:
                break

            if lastLabels == []:
                lastLabels=trajectoire.findFrame(frame)

            if frame == lastFrame:
                continue

            nextLabels = trajectoire.findFrame(nextFrame)

            if not frame in solutions:
                solutions[frame]={}

            l = len(lastLabels)
            nextL = len(nextLabels)
            #print "FRAME :"+str(frame),l, nextL, lastLabels, nextLabels, trajectoire.numCellule


            if int(l)==int(nextL):
                if int(l)==1 and int(lastLabels[0])==-1:
                    appearing = nextLabels[0]
                    if not nextFrame in appearCandidates:
                        appearCandidates[nextFrame]=[]
                    appearCandidates[nextFrame].append(appearing)

                elif int(l)==1 and int(nextLabels[0])==-2:
                    disappearing = lastLabels[0]
                    if not frame in disappearCandidates:
                        disappearCandidates[frame]=[]
                    disappearCandidates[frame].append(disappearing)
                    #continue

                elif int(l)==1:
                    #print "MOVE"
                                    
                    if not "move" in solutions[frame]:
                        solutions[frame]["move"]=[[], []]

                    solutions[frame]["move"][0].append(lastLabels)
                    solutions[frame]["move"][1].append(nextLabels)
                    if not nextFrame in coming:
                        coming[nextFrame]=[]
                    for label in nextLabels: 
                        coming[nextFrame].append(label)
                        
                    if not frame in going:
                        going[frame]=[]
                    for label in lastLabels:
                        going[frame].append(label)

                else:
                    print "problem", lastLabels, nextLabels, frame
                    raise

            else:
                print "problem at frame "+str(frame), lastLabels, nextLabels
                raise

            lastFrame = frame
            if not lastMoveSplit[frame]:
                 lastLabels = nextLabels
            if int(lastLabels[0])==-2:
                break
            
#looking at merges : if two cells have the same target, it means it is not a move but a merge

    for f in solutions:
        if "move" not in solutions[f]:
            continue
        
        listeLabelsTargetMove = solutions[f]["move"][1]
        listeLabelsSourceMove = solutions[f]["move"][0]
        nextFrame = int(f+1)
         #print listeLabelsSourceMove, listeLabelsTargetMove
        
        newSourceMove = []
        newTargetMove = []
        
        merges = []
        moves = []
        for t_el in listeLabelsTargetMove:
            if listeLabelsTargetMove.count(t_el) > 1:
                merges.append(t_el)
            else:
                moves.append(t_el)

        indexToDel = []
        for label1 in merges:
            #print label1
            nextLabels=label1
            lastLabels=[]
            for i in range(len(listeLabelsTargetMove)):
                if listeLabelsTargetMove[i]==label1:
                    lastLabels.append(listeLabelsSourceMove[i][0])
                    indexToDel.append(i)
#                    
            for l in lastLabels:
                if lastLabels.count(l)>1:
                    print f, l
                    raise SameLabelException

            if not "merge" in solutions[f]:
                solutions[f]["merge"]=[[], []]
            if nextLabels in solutions[f]["merge"][1]:
                continue
            solutions[f]["merge"][0].append(lastLabels)
            solutions[f]["merge"][1].append(nextLabels)

            if not f in going:
                going[f]=[]
            if not nextFrame in coming:
                coming[nextFrame]=[]
            for label in nextLabels:
                coming[nextFrame].append(label)
            for label in lastLabels:
                going[f].append(label)
        
#looking at splits : if two cells have the same source, it means it is not a move but a split
        
        splits = []
        moves = []
        for s_el in listeLabelsSourceMove:
            if listeLabelsSourceMove.count(s_el) > 1:
                splits.append(s_el)
            else:
                moves.append(s_el)
        
#        indexToDel = []
        for label1 in splits:
            #print label1
            lastLabels=label1
            nextLabels=[]
            for i in range(len(listeLabelsSourceMove)):
                if listeLabelsSourceMove[i]==label1:
                    indexToDel.append(i)
                    nextLabels.append(listeLabelsTargetMove[i][0])
                    
            for l in nextLabels:
                if nextLabels.count(l)>1:
                    print f, l
                    raise SameLabelException
#                    
            if not "split" in solutions[f]:
                solutions[f]["split"]=[[], []]
            if lastLabels in solutions[f]["split"][0]:
                continue
#            if maxTwos and len(nextLabels)>2:
#                continue
            solutions[f]["split"][0].append(lastLabels)
            solutions[f]["split"][1].append(nextLabels)

            if not f in going:
                going[f]=[]
            if not nextFrame in coming:
                coming[nextFrame]=[]
            for label in nextLabels:
                coming[nextFrame].append(label)
            for label in lastLabels:
                going[f].append(label)
                
        for i in range(len(listeLabelsSourceMove)):
            if i not in indexToDel:
                newSourceMove.append(listeLabelsSourceMove[i])
                newTargetMove.append(listeLabelsTargetMove[i])
 #       print listeLabelsSourceMove, listeLabelsTargetMove
        #print newSourceMove, newTargetMove
        solutions[f]["move"][1] = newTargetMove
        solutions[f]["move"][0] = newSourceMove
        
 
#appear candidates 
    for frame in appearCandidates.keys():
        appearlist = appearCandidates[frame]
        for label in appearlist:
            if frame not in coming.keys() or (frame in coming.keys() and label not in coming[frame]):
                if not "appear" in solutions[frame-1]:
                    solutions[frame-1]["appear"]=[[], []]
                if [label] in solutions[frame-1]["appear"][1]:
                    continue
                #print "APPEAR "+str(label)+" on frame "+str(frame)
                solutions[frame-1]["appear"][0].append([-1])
                solutions[frame-1]["appear"][1].append([label])

#disappear candidates
    for frame in disappearCandidates.keys():
        disappearlist = disappearCandidates[frame]
        for label in disappearlist:
            if frame not in going.keys() or (frame in going.keys() and label not in going[frame]):
                if not "disappear" in solutions[frame]:
                    solutions[frame]["disappear"]=[[], []]
                if [label] in solutions[frame]["disappear"][0]:
                    continue
                #print "DISAPPEAR "+str(label)+" on frame "+str(frame)
                solutions[frame]["disappear"][0].append([label])
                solutions[frame]["disappear"][1].append([-1])

    count={}
    out_merge = "\n MERGE \n"
    out_split = "\n SPLIT \n"
    for e in EVENTS:
        count[e]=0

    for f in solutions:
        #fileS.write("\n --------------------------FRAME "+str(f)+"------------------------------------")
        for e in solutions[f]:
         #   fileS.write("\n EVENT "+str(e)+"*******\n SOURCES :\n")
            j=0
            for label in solutions[f][e][0]:
          #      fileS.write("\n label"+str(label))
                if label == []:
                    continue
                j+=1
                count[e]+=1
            i=0
           # fileS.write("\n TARGETS :\n")
            for label in solutions[f][e][1]:
            #    fileS.write("\n label"+str(label))
                i+=1

            if i<>j: print "probleme : difference de longueurs entre Source et Target at fr "+str(f)+" pour l'evenement "+e

            if e =="merge":
                out_merge +="\n"+str(solutions[f][e][0])+" on frame "+str(f)+" to "+str(solutions[f][e][1])
            if e=="split":
                out_split +="\n"+str(solutions[f][e][0])+" on frame "+str(f)+" to "+str(solutions[f][e][1])
    #fileS.write(str(count))
    print count#, out_merge, out_split
    #fileS.close()
    return solutions

#def ecrireTrainingMergesSplitsOnly(ensembleTraj, fileToSave, LENGTH):
#    solutions = {}
#    fileS=open(fileToSave, "w")
#
#    going = {}
#    coming = {}
#    appearCandidates = {}
#    disappearCandidates = {}
#
#    for trajectoire in ensembleTraj.lstTraj:
#        frameLabels = sorted(trajectoire.lstPoints.iterkeys())
#        lastLabels = []
#        lastFrame = -1
#        lastMoveSplit = [False for x in range(100)]
#
#        for f in frameLabels:
#            frame = int(f[0])
#            nextFrame = frame+1
#
#            if nextFrame==LENGTH:
#                break
#
#            if lastLabels == []:
#                lastLabels=trajectoire.findFrame(frame)
#
#            if frame == lastFrame:
#                continue
#
#            nextLabels = trajectoire.findFrame(nextFrame)
#
#            if not frame in solutions:
#                solutions[frame]={}
#
#            l = len(lastLabels)
#            nextL = len(nextLabels)
#            #print "FRAME :"+str(frame),l, nextL, lastLabels, nextLabels, trajectoire.numCellule
#
#
#            if int(l)==int(nextL):
#                if int(l)==1 and int(lastLabels[0])==-1:
#                    appearing = nextLabels[0]
#                    if not nextFrame in appearCandidates:
#                        appearCandidates[nextFrame]=[]
#                    appearCandidates[nextFrame].append(appearing)
#
#                elif int(l)==1 and int(nextLabels[0])==-2:
#                    disappearing = lastLabels[0]
#                    if not frame in disappearCandidates:
#                        disappearCandidates[frame]=[]
#                    disappearCandidates[frame].append(disappearing)
#                    #continue
#
#                elif int(l)==1:
#                    #print "MOVE"
#                                    
#                    if not "move" in solutions[frame]:
#                        solutions[frame]["move"]=[[], []]
#
#                    solutions[frame]["move"][0].append(lastLabels)
#                    solutions[frame]["move"][1].append(nextLabels)
#                    if not nextFrame in coming:
#                        coming[nextFrame]=[]
#                    for label in nextLabels: 
#                        coming[nextFrame].append(label)
#                        
#                    if not frame in going:
#                        going[frame]=[]
#                    for label in lastLabels:
#                        going[frame].append(label)
#
#                else:
#                    print "problem", lastLabels, nextLabels, frame
#                    raise
#
#            else:
#                print "problem at frame "+str(frame), lastLabels, nextLabels
#                raise
#
#            lastFrame = frame
#            if not lastMoveSplit[frame]:
#                 lastLabels = nextLabels
#            if int(lastLabels[0])==-2:
#                break
#            
##looking at merges : if two cells have the same target, it means it is not a move but a merge
#
#    for f in solutions:
#        if "move" not in solutions[f]:
#            continue
#        
#        listeLabelsTargetMove = solutions[f]["move"][1]
#        listeLabelsSourceMove = solutions[f]["move"][0]
#        nextFrame = int(f+1)
#         #print listeLabelsSourceMove, listeLabelsTargetMove
#        
#        newSourceMove = []
#        newTargetMove = []
#        
#        merges = []
#        moves = []
#        for t_el in listeLabelsTargetMove:
#            if listeLabelsTargetMove.count(t_el) > 1:
#                merges.append(t_el)
#            else:
#                moves.append(t_el)
#
#        indexToDel = []
#        for label1 in merges:
#            #print label1
#            nextLabels=label1
#            lastLabels=[]
#            for i in range(len(listeLabelsTargetMove)):
#                if listeLabelsTargetMove[i]==label1:
#                    lastLabels.append(listeLabelsSourceMove[i][0])
#                    indexToDel.append(i)
##                    
#
#            if not "merge" in solutions[f]:
#                solutions[f]["merge"]=[[], []]
#            if nextLabels in solutions[f]["merge"][1]:
#                continue
##            if maxTwos and len(lastLabels)>2:
##                continue
#            solutions[f]["merge"][0].append(lastLabels)
#            solutions[f]["merge"][1].append(nextLabels)
#
#            if not f in going:
#                going[f]=[]
#            if not nextFrame in coming:
#                coming[nextFrame]=[]
#            for label in nextLabels:
#                coming[nextFrame].append(label)
#            for label in lastLabels:
#                going[f].append(label)
#        
##looking at splits : if two cells have the same source, it means it is not a move but a split
#        
#        splits = []
#        moves = []
#        for s_el in listeLabelsSourceMove:
#            if listeLabelsSourceMove.count(s_el) > 1:
#                splits.append(s_el)
#            else:
#                moves.append(s_el)
#        
##        indexToDel = []
#        for label1 in splits:
#            #print label1
#            lastLabels=label1
#            nextLabels=[]
#            for i in range(len(listeLabelsSourceMove)):
#                if listeLabelsSourceMove[i]==label1:
#                    indexToDel.append(i)
#                    nextLabels.append(listeLabelsTargetMove[i][0])
##                    
#            if not "split" in solutions[f]:
#                solutions[f]["split"]=[[], []]
#            if lastLabels in solutions[f]["split"][0]:
#                continue
##            if maxTwos and len(nextLabels)>2:
##                continue
#            solutions[f]["split"][0].append(lastLabels)
#            solutions[f]["split"][1].append(nextLabels)
#
#            if not f in going:
#                going[f]=[]
#            if not nextFrame in coming:
#                coming[nextFrame]=[]
#            for label in nextLabels:
#                coming[nextFrame].append(label)
#            for label in lastLabels:
#                going[f].append(label)
#                
#        for i in range(len(listeLabelsSourceMove)):
#            if i not in indexToDel:
#                newSourceMove.append(listeLabelsSourceMove[i])
#                newTargetMove.append(listeLabelsTargetMove[i])
# #       print listeLabelsSourceMove, listeLabelsTargetMove
#        #print newSourceMove, newTargetMove
#        solutions[f]["move"][1] = newTargetMove
#        solutions[f]["move"][0] = newSourceMove
#        
## 
###appear candidates 
##    for frame in appearCandidates.keys():
##        appearlist = appearCandidates[frame]
##        for label in appearlist:
##            if frame not in coming.keys() or (frame in coming.keys() and label not in coming[frame]):
##                if not "appear" in solutions[frame-1]:
##                    solutions[frame-1]["appear"]=[[], []]
##                if [label] in solutions[frame-1]["appear"][1]:
##                    continue
##                print "APPEAR "+str(label)+" on frame "+str(frame)
##                solutions[frame-1]["appear"][0].append([-1])
##                solutions[frame-1]["appear"][1].append([label])
##
###disappear candidates
##    for frame in disappearCandidates.keys():
##        disappearlist = disappearCandidates[frame]
##        for label in disappearlist:
##            if frame not in going.keys() or (frame in going.keys() and label not in going[frame]):
##                print "DISAPPEAR "+str(label)+" on frame "+str(frame)
##
##                if not "disappear" in solutions[frame]:
##                    solutions[frame]["disappear"]=[[], []]
##                if [label] in solutions[frame]["disappear"][0]:
##                    continue
##
##                solutions[frame]["disappear"][0].append([label])
##                solutions[frame]["disappear"][1].append([-1])
#
#    count={}
#    out_merge = "\n MERGE \n"
#    out_split = "\n SPLIT \n"
#    for e in EVENTS:
#        count[e]=0
#
#    for f in solutions:
#        fileS.write("\n --------------------------FRAME "+str(f)+"------------------------------------")
#        for e in solutions[f]:
##            pdb.set_trace()
#            if e!="merge" and e!='split':
#                solutions[f][e]=None
#                continue
#            
#            fileS.write("\n EVENT "+str(e)+"*******\n SOURCES :\n")
#            j=0
#            for label in solutions[f][e][0]:
#                fileS.write("\n label"+str(label))
#
#                if label == []:
#                    continue
#                j+=1
#                count[e]+=1
#            i=0
#            fileS.write("\n TARGETS :\n")
#            for label in solutions[f][e][1]:
#                fileS.write("\n label"+str(label))
#
#                i+=1
#
#            if i<>j: print "probleme : difference de longueurs entre Source et Target at fr "+str(f)+" pour l'evenement "+e
#
#            if e =="merge":
#                out_merge +="\n"+str(solutions[f][e][0])+" on frame "+str(f)+" to "+str(solutions[f][e][1])
#            if e=="split":
#                out_split +="\n"+str(solutions[f][e][0])+" on frame "+str(f)+" to "+str(solutions[f][e][1])
#    fileS.write(str(count))
#    print count, out_merge, out_split
#    fileS.close()
#    return solutions


def initXml():
    return "<?xml version=\"1.0\" encoding=\"utf8\"?>\n  <CellCounter_Marker_File>\n    <Marker_Data>"

def finirXml():
    return "\n    </Marker_Data>\n  </CellCounter_Marker_File> "

def ecrireXml(ligne, trajectoire):
    result = "\n      <Marker_Type>\n        <Type>{0}</Type>".format(ligne)
    
    for point in trajectoire.lstPoints:
        numFramePourXml = point[0]+1
        
        result+="\n        <Marker>\n          <MarkerX>{1}</MarkerX>\n          <MarkerY>{2}</MarkerY>\n          <MarkerZ>{0}</MarkerZ>\n        </Marker>".format(numFramePourXml, trajectoire.lstPoints[point][0], trajectoire.lstPoints[point][1])
        
    result+="\n      </Marker_Type>"
    return result


def arraySegmentation(plaque, puits, dir):
    fileHDF5 = dir[:-8]+"hdf5/"+puits+".hdf5"
    path = "/sample/0/plate/"+plaque+"/experiment/"+puits[:-3]+"/position/"+puits[-1]+"/image/region"
    tab = vi.readHDF5(fileHDF5, path)

    return tab


def HacktrackReverse(plaque, puits, dir):
#    global LENGTH
    dossier = dir+"/"+puits+"/statistics"
    boolean = False
    tab = arraySegmentation(plaque, puits, dir)
    print str(tab.shape)
    LENGTH = tab.shape[1]-1
    
    if "pheno" in dir:
        boolean = True
        sauvegarde = "/input data/pheno/"
    else:
        sauvegarde = "/input data/controls/"
    print "Pheno?"+str(boolean), "Numero de plaque :"+plaque+" numero de puits : "+ puits

#    nomDot = dossier + ("/tracking_graph___P%s.dot" % puits)    
    nomTxt =  dossier + ("/tracking_graph___P%s_features.txt" % puits)

    nomFichier = "PL"+plaque+"___P"+puits+"___T00000.xml"
    nomFichier = ds+"/"+nomFichier
    print os.getcwd()

    print nomTxt, nomFichier

    ensembleTraj = lireTrajectoires(nomFichier, tab)
    if ensembleTraj == None:
        return
    fileToSave =  ("%s_P%s.txt" % (plaque, puits))
    #fileToSave=sauvegarde+fileToSave
    solutions = ecrireTrainingSansHDF5(ensembleTraj, fileToSave, LENGTH)

  #      ensemTraj = lireDot(nomDot, X, Y)
   #     print ensemTraj.numMitoses
    #    sortiesCarre(nomFichier, ensemTraj)
   # except IOError as e :
    #    print "probleme fichier"+e.message


    return solutions
