import sys, os, os.path, random, csv, pdb
import numpy as np
import cPickle as pickle

couleurs = []
couleurs.append("A50026")
couleurs.append("D73027")
couleurs.append("F46D43")
couleurs.append("FDAE61")
couleurs.append("FEE090")
couleurs.append("FFFFBF")
couleurs.append("E0F3F8")
couleurs.append("ABD9E9")
couleurs.append("74ADD1")
couleurs.append("4575B4")
couleurs.append("ABDDA4")
couleurs.append("A50026")
couleurs.append("D73027")
couleurs.append("F46D43")
couleurs.append("FDAE61")
couleurs.append("FEE090")
couleurs.append("FFFFBF")
couleurs.append("E0F3F8")
couleurs.append("ABD9E9")
couleurs.append("74ADD1")
couleurs.append("4575B4")
couleurs.append("ABDDA4")
couleurs.append("A50026")
couleurs.append("D73027")
couleurs.append("F46D43")
couleurs.append("FDAE61")
couleurs.append("FEE090")
couleurs.append("FFFFBF")
couleurs.append("E0F3F8")
couleurs.append("ABD9E9")
couleurs.append("74ADD1")
couleurs.append("4575B4")
couleurs.append("ABDDA4")
couleurs.append("A50026")
couleurs.append("D73027")
couleurs.append("F46D43")
couleurs.append("FDAE61")
couleurs.append("FEE090")
couleurs.append("FFFFBF")
couleurs.append("E0F3F8")
couleurs.append("ABD9E9")
couleurs.append("74ADD1")
couleurs.append("4575B4")
couleurs.append("ABDDA4")
couleurs.append("A50026")
couleurs.append("D73027")
couleurs.append("F46D43")
couleurs.append("FDAE61")
couleurs.append("FEE090")
couleurs.append("FFFFBF")
couleurs.append("E0F3F8")
couleurs.append("ABD9E9")
couleurs.append("74ADD1")
couleurs.append("4575B4")
couleurs.append("ABDDA4")
couleurs.append("A50026")
couleurs.append("D73027")
couleurs.append("F46D43")
couleurs.append("FDAE61")
couleurs.append("FEE090")
couleurs.append("FFFFBF")
couleurs.append("E0F3F8")
couleurs.append("ABD9E9")
couleurs.append("74ADD1")
couleurs.append("4575B4")
couleurs.append("ABDDA4")

dossierSauvegarde = "annotations"

def coord(nomTxt):
    #print repr(fTxt)
    fTxt = open(nomTxt, "r")
    fTxt.readline()
    X = {}
    Y = {}
    upX = {}
    upY = {} 
    dX = {}
    dY = {}
    entier = 0
    for line in fTxt.readlines():
    #print line
        maLigne = line.split("\t")
       # print maLigne
#ATTENTION CES PARAM DEPENDENT DES DONNEES QUI ONT ETE EXPORTEES
        X[(int(maLigne[0]), int(maLigne[1]))]=int(maLigne[-6])
        Y[(int(maLigne[0]), int(maLigne[1]))]=int(maLigne[-5])
        upX[(int(maLigne[0]), int(maLigne[1]))]=int(maLigne[-4])
        upY[(int(maLigne[0]), int(maLigne[1]))]=int(maLigne[-3])
        dX[(int(maLigne[0]), int(maLigne[1]))]=int(maLigne[-2])
        dY[(int(maLigne[0]), int(maLigne[1]))]=int(maLigne[-1])
 #       print Y
 #       entier +=1
 #       if entier > 9: break
    fTxt.close()
    return (X,Y, upX, upY, dX, dY)

def lireDot(nomDot, X, Y):
    fDot = open(nomDot, "r")

    mesLignes=fDot.readlines()
    fDot.close()
    del mesLignes[0:3]
    mesLignes.reverse()

    continuer = True
    ensembleTraj = ensTraj()
    print ensembleTraj.numMitoses
    print ensembleTraj.lstTraj

    trajCourante = trajectoire(0,0,0,0,0)
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
        
        if (int(frameCourante) != int(frameSuivante) or int(idCourante) != int(idSuivante)):
            if frameSuivante!=-1 and idSuivante!=-1:
                trajCourante.ajoutPoint(X[(frameSuivante, idSuivante)], Y[(frameSuivante, idSuivante)], frameSuivante, idSuivante)
              
            trajectoires = ensembleTraj.trouver(trajCourante.numCellule)
            lstPointsTotale = []
            for trajec in trajectoires.lstTraj:
                for element in trajec.lstPoints.keys():
                    lstPointsTotale.append(element)

            if (frameCourante, idCourante) in lstPointsTotale or (frameCourante, idCourante) in trajCourante.lstPoints.keys():
                trajCourante.mitose = 1
                ensembleTraj.numMitoses +=1
                if (frameCourante, idCourante) in trajCourante.lstPoints.keys():
                    trajCourante.ajoutPoint(X[(frameCourante, idCourante)], Y[(frameCourante, idCourante)], frameCourante, idCourante)
                ensembleTraj.ajouter(trajCourante)
                num = trajCourante.numCellule
                trajCourante.supp()
                trajCourante = trajectoire(num, X[(frameCourante, idCourante)], Y[(frameCourante, idCourante)], frameCourante, idCourante)
                trajCourante.mitose = 1

            else:
                ensembleTraj.ajouter(trajCourante)
                num = trajCourante.numCellule+1

                trajCourante.supp()

                trajCourante = trajectoire(num, X[(frameCourante, idCourante)], Y[(frameCourante, idCourante)], frameCourante, idCourante)             

        else:
            trajCourante.ajoutPoint(X[(frameCourante, idCourante)], Y[(frameCourante, idCourante)], frameCourante, idCourante)
     

        frameSuivante = int(decomp[2].split("_")[0])
        idSuivante = int(decomp[2].split("_")[1])

    ensembleTraj.ajouter(trajCourante)

    return ensembleTraj

def sortiesCarre(nomFichier, ensembleTraj):

    compteur = len(couleurs)
    MIN = 200
    absmax=300
    ordmax = 300
    MAX = 1000
    PAS = 150
    iterationX = 0
    iterationY=0

    fichierX = open(nomFichier, "w")
    fichierX.write(initXml())
 
    nbDansBoiteMax = 0
    absFinale = 0
    ordFinale = 0

    while absmax <= MAX:
        absmin = MIN + iterationX*PAS
        absmax = MIN + (iterationX+1)*PAS
        while ordmax<=MAX:
            ordmin = MIN + iterationY*PAS
            ordmax = MIN + (iterationY+1)*PAS
            nbDansBoite = 0
            listeCellulesVues = []
                
            for trajectoire in ensembleTraj.lstTraj:
                keyList = trajectoire.lstPoints.keys()
        #try:
                debutFrame =100
                (xi, yi)=(0,0)
                for i in keyList:
                    if int(i[0])<=debutFrame:
                        debutFrame = int(i[0])
                        debutId = i[1]
                        (xi, yi)=trajectoire.lstPoints[(debutFrame, debutId)]

                    if xi>absmin and xi<absmax and yi>ordmin and yi<ordmax:
                        #print "Cellule numero "+str(trajectoire.numCellule)+" dans la boite ("+str(absmin)+", "+str(ordmin)+")"
                        nbDansBoite +=1
                        if trajectoire.numCellule not in listeCellulesVues:
                            listeCellulesVues.append(trajectoire.numCellule)
           # print listeCellulesVues
               
            if nbDansBoite > nbDansBoiteMax:
                absFinale = absmin
                ordFinale = ordmin
                #print listeCellulesVues
                listeCellulesVuesMax = []
                for num in listeCellulesVues : listeCellulesVuesMax.append(num)
 #               listeCelluleVuesMax = listeCellulesVues
                #print listeCellulesVuesMax
                nbDansBoiteMax = nbDansBoite
            iterationY+=1
        iterationY=0
        ordmax = 300
        iterationX+=1

    print "Finalement : "+str(nbDansBoiteMax)+" dans la boite ("+str(absFinale)+", "+str(ordFinale)+")"
    print listeCellulesVuesMax

   
    for numCellule in listeCellulesVuesMax:
        print numCellule
        trajectoires = ensembleTraj.trouver(numCellule)
        try:
            for trajectoire in trajectoires.lstTraj:
                if compteur >0:
                    fichierX.write(ecrireXml(compteur,trajectoire))
                    compteur-=1
                else: 
                    break
        except:
            print "il y a un probleme entre les numeros des cellules trouvees et les cellules presentes dans la liste de trajectoires."

        #except:
         #   print keyList
    fichierX.write(finirXml())
    fichierX.close()
    return

def sortiesAll(nomFichier, ensembleTraj):
#pour sortir toutes les predictions Cecog
    compteur =1
    ecrites = []
    fichierX = open(nomFichier, "w")
    fichierX.write(initXml())
    
    for trajectoire in ensembleTraj.lstTraj:
        txt, coord = ecrireXml(compteur,trajectoire.lstPoints, True)
        fichierX.write(txt);ecrites.append(coord)
        compteur+=1
        if compteur>1600:
            break

    fichierX.write(finirXml())
    fichierX.close()
#    fichierCoordonnees3D=open(nomFichier[:-3]+'pkl', 'w')
#    pickle.dump(ecrites, fichierCoordonnees3D)
#    fichierCoordonnees3D.close()
    print 'finally ', compteur, 'trajectories'
    return ecrites

def sortiesAllBis(nomFichier, ensembleTraj, compteurMax):

    compteur = 5
    listeCellulesVues = []
    first = True
    ecrites = []
#
#    for trajectoire in ensembleTraj.lstTraj:
#        print "trajectoire :", trajectoire.numCellule, ' et ses points'
#        for point in sorted(trajectoire.lstPoints): print point, trajectoire.lstPoints[point]

    fichierX = open(nomFichier, "w")
    fichierX.write(initXml())

#    for trajectoire in filter(lambda x: x.mtsInteressants !=[] and x.numCellule not in listeCellulesVues, ensembleTraj.lstTraj):
#        if compteur >0:
#            fichierX.write(ecrireXml(compteur, trajectoire))
#            #fichierX.write(ecrireXml(compteurM,trajectoire))
#            listeCellulesVues.append(trajectoire.numCellule)
#            compteur -=1
#            for index in trajectoire.mtsInteressants:
#                x, y = trajectoire.mtsInteressants[index]
#                for traj in filter(lambda x: x.numCellule not in listeCellulesVues, ensembleTraj.lstTraj):
#                    if traj.findFrame(index, x, y) and compteur >0:
#                        fichierX.write(ecrireXml(compteur, traj))
#                #fichierX.write(ecrireXml(compteurM,trajectoire))
#                        listeCellulesVues.append(traj.numCellule)
#                        #compteur -=1
#
#            compteur -=1
            
    for trajectoire in ensembleTraj.lstTraj: 
#        if trajectoire.numCellule in [14,18,30,38,22]:
#            continue
        keyList = trajectoire.lstPoints.keys()
        debutFrame =100
        (xi, yi)=(0,0)
        for i in keyList:
            if int(i[0])<=debutFrame:
                debutFrame = int(i[0])
                debutId = i[1]
                (xi, yi)=trajectoire.lstPoints[(debutFrame, debutId)]
        #print debutFrame, debutId, xi, yi
        
        if xi>800 and yi<400 and trajectoire.numCellule not in listeCellulesVues:
#        if xi>500 and xi<800 and yi>200 and yi<600 and trajectoire.numCellule not in listeCellulesVues:    
            #print "une trajectoire : compteur ", compteur, " numCellule", trajectoire.numCellule
        #if (trajectoire.mitose == 1 or len(trajectoire.lstPoints)<3) and trajectoire.numCellule not in listeCellulesVues and compteur>0:
            txt, coordonnees = ecrireXml(compteur, trajectoire, sortie=True)
            fichierX.write(txt); ecrites.append(coordonnees)
            listeCellulesVues.append(trajectoire.numCellule) 
            compteur +=1
            
            for traj in filter(lambda x: x.numCellule == trajectoire.numCellule and x!=trajectoire, ensembleTraj.lstTraj):
                #pdb.set_trace()
                #print "ajout de ses soeurs"
                txt, coordonnees = ecrireXml(compteur, traj, sortie=True)
                fichierX.write(txt);ecrites.append(coordonnees)
                compteur +=1
    
#            for traj in filter(lambda x: x.numCellule not in listeCellulesVues, ensembleTraj.lstTraj):
#                frame = min(trajectoire.mtsInteressants.keys()); x, y = trajectoire.mtsInteressants[frame]
#                if traj.findFrame(frame, x, y):
#                    fichierX.write(ecrireXml(compteur, trajectoire))
#            #fichierX.write(ecrireXml(compteurM,trajectoire))
#
#                    listeCellulesVues.append(trajectoire.numCellule)
#            compteur -=1

    fichierX.write(finirXml())
    fichierX.close()
    fichierCoordonnees3D=open(nomFichier[:-3]+'pkl', 'w')
    pickle.dump(ecrites, fichierCoordonnees3D)
    fichierCoordonnees3D.close()
    
    if compteur>compteurMax:
        print "compteur final", compteur, "donc je cree classdef en consequence"
        ecrireClassDef(compteur)
    return max(compteur, compteurMax)

def ecrireClassDef(longueur = None):
    fichierT = open("class_definition.txt", "w")
    ecrivain = csv.writer(fichierT, delimiter='\t', quoting=csv.QUOTE_NONE)
    if longueur == None:
        longueur = len(couleurs)

    for i in range(longueur):
        bou = i/len(couleurs)
        ecrivain.writerow([i+1, "Cellule"+str(i+1), "#"+couleurs[i-len(couleurs)*bou]])

    fichierT.close()
    return

def sortiesRandom(nomFichier, ensembleTraj):

    compteurM = 5
    compteurN = 6
    listeCellulesMitosesVues = []

    fichierX = open(nomFichier, "w")
    fichierX.write(initXml())

    for trajectoire in ensembleTraj.lstTraj:

        if trajectoire.mitose == 1 and trajectoire.numCellule not in listeCellulesMitosesVues and compteurM>0:
            randM = random.randint(1,10)
            if randM == 1:

                fichierX.write(ecrireXmlM(trajectoire.numCellule, compteurM, ensembleTraj))
                #fichierX.write(ecrireXml(compteurM,trajectoire))

                listeCellulesMitosesVues.append(trajectoire.numCellule)
                compteurM -=1

        elif trajectoire.numCellule not in listeCellulesMitosesVues and compteurN>0:
            randN = random.randint(1,10)
            if randN == 1:
                fichierX.write(ecrireXml(compteurN+5,trajectoire))

                compteurN -=1

    fichierX.write(finirXml())
    fichierX.close()
    return

   


def initXml():
    return "<?xml version=\"1.0\" encoding=\"utf8\"?>\n  <CellCounter_Marker_File>\n    <Marker_Data>"

def finirXml():
    return "\n    </Marker_Data>\n  </CellCounter_Marker_File> "

def ecrireXml(ligne, lstPoints, sortie=False):
#typiquement un truc qui devrait etre une methode de la class trajectoire haha
    coordonnees = {}
    result = "\n      <Marker_Type>\n        <Type>{0}</Type>".format(ligne)

    for point in lstPoints:
        coordonnees.update({point[0] :lstPoints[point]});
        numFramePourXml = point[0]+1
        
        result+="\n        <Marker>\n          <MarkerX>{1}</MarkerX>\n          <MarkerY>{2}</MarkerY>\n          <MarkerZ>{0}</MarkerZ>\n        </Marker>".format(numFramePourXml, lstPoints[point][0], lstPoints[point][1])
        
    result+="\n      </Marker_Type>"
    if sortie:
        x=[]; y = []; z=coordonnees.keys()
        for frame in z:
            x.append(coordonnees[frame][0])
            y.append(coordonnees[frame][1])
        return result, [x, y, z]
    else:
        return result

class trajectoire():

    def __init__(self, num, xi=None, yi=None, frame=None, idC=None, id=None):
   #     print "creation d'une trajectoire"
        self.id=id
        self.numCellule = num
        self.fusion = False
        self.mitose = -1
        self.lstPoints = {}
        self.density = None
#        self.fractions=[]
        if xi is not None and yi is not None and frame is not None and idC is not None:
            self.lstPoints[(frame, idC)]=(xi, yi)
            
    def addDensity(self, center, tree, nb_voisins_max, densities):
        #print 'id ', self.id
        #print 'longueur traj ', len(self.lstPoints.keys())
        r=[]
        m=max(densities)
        dou, ind=tree.query(center, k=nb_voisins_max, distance_upper_bound=m)
        for i, dens in enumerate(densities):
            li = filter(lambda k: k in np.where(dou>0)[0] and k in np.where(dou<dens)[0], range(10))
            r.append(len(li))
        r=np.array(r)
        self.density = np.vstack((self.density, r)) if self.density is not None else r
        #print 'density', self.density

    def ajoutPoint(self, x, y, frame, idC):
        self.lstPoints[(frame, idC)]=(x,y)
        
    def firstFrame(self):
        return sorted(self.lstPoints.keys(), key=lambda tup:tup[0])[0][0]
        
    def findLast(self):
        keyList = sorted(self.lstPoints.keys(), key=lambda tup:tup[0])
#        lastFrame =0
#        (xi, yi)=(0,0)
#         for i in keyList:
#             if int(i[0])>=lastFrame:
#                 lastFrame = int(i[0])
# #                debutId = i[1]
# #                (xi, yi)=self.lstPoints[(debutFrame, debutId)]
#         if lastFrame !=keyList[-1][0]:
#             raise
#         return lastFrame
        return keyList[-1][0]

    def copyBefore(self, index):
        newKeys = filter(lambda x: x[0]<index+1, self.lstPoints)
        if newKeys==[]:
            pdb.set_trace()
            return None
        lstPoints = {k:self.lstPoints[k] for k in newKeys}
        newTraj = trajectoire(self.numCellule)
        newTraj.lstPoints = lstPoints
        newTraj.mitose=filter(lambda x: x<index, self.mitose) if self.mitose !=-1 else -1
        newTraj.fusion=filter(lambda x: x<index, self.fusion) if self.fusion !=-1 else -1
        newTraj.density=self.density[:index-min(self.lstPoints)[0]]
        pdb.set_trace()
        return newTraj
    
    def copyBetween(self, beginning, end):
        '''
        Function to copy a trajectory between beginning (included) to end (included)
        '''
        newKeys = filter(lambda x: x[0]>=beginning and x[0]<end+1, self.lstPoints)
        lstPoints = {k:self.lstPoints[k] for k in newKeys}
        newTraj = trajectoire(self.numCellule)
        newTraj.lstPoints = lstPoints
        newTraj.mitose=self.mitose#filter(lambda x: x>beginning and x<end, self.mitose) if self.mitose !=-1 else -1
        newTraj.fusion=self.fusion#filter(lambda x: x>beginning and x<end, self.fusion) if self.fusion!=-1 else -1
        newTraj.density=self.density[max(0, beginning-min(self.lstPoints)[0]):end-min(self.lstPoints)[0]+1]
        if min(self.lstPoints)[0]>=beginning and max(self.lstPoints)[0]<end:
            pdb.set_trace()
        return newTraj

    def findFrame(self, frame, x=None, y = None):
        result = []
        if x==None:
            for point in self.lstPoints:
                if int(point[0])==int(frame):
                    result.append(point[1])
            return result
        else:
           for point in self.lstPoints:
                if int(point[0])==int(frame) and (self.lstPoints[point][0]-x)**2+(self.lstPoints[point][1]-y)**2<50**2:
                    return True 
    def __get__(self):
        return self.lstPoints
    
    def split(self, frame):
        if self.mitose<0:
            self.mitose=[frame]
        else:
            self.mitose.append(frame)
        #print self.mitose
    
    def merge(self, frame):
        if self.fusion<0:
            self.fusion=[frame]
        else:
            self.fusion.append(frame)

    def supp(self):
    #    print "suppression de trajectoire"
        self.lstPoints = {}
        self.mitose = -1        

class ensTraj(object):
    def __init__(self):
        self.numMitoses = 0
        self.lstTraj = []
        
    def __get__(self):
        return self.lstTraj

    def ajouter(self, traj):
     #   print "ajout d'une trajectoire"
        mit = traj.mitose
        num = traj.numCellule
        points = traj.lstPoints
        copieTraj = trajectoire(num, 0,0,0,0)
        copieTraj.mitose = mit
        for point in points:
            copieTraj.lstPoints[point]=points[point]
        del copieTraj.lstPoints[(0,0)]
       # liste
        self.lstTraj.append(copieTraj)
    
    def add(self, traj):
     #   print "ajout d'une trajectoire"
        self.lstTraj.append(traj)

    def trouver(self, numC):
        result = ensTraj()
        for traj in self.lstTraj:
            if traj.numCellule ==numC:
                #print "c'est fait de trouver la cellule "+str(numC)
                result.ajouter(traj)
        return result
    
    def findLabelFrame(self, label, index):
        result = []
        if index>=0:
            for traj in self.lstTraj:
                point = traj.findFrame(index)
                if label in point:
                    result.append(traj)
        return result
                    
def Hacktrack(plaque, puits, dir):
   #dir = "/media/lalil0u/New/projects/groundtruth/results"+/pheno ou /controls
   
    dossier = dir+plaque+"/analyzed/"+puits+"/statistics"
    print "Numero de plaque :"+plaque+" numero de puits : "+ puits
    print "Dossier de sauvegarde :"+dir+"tracks_annotations/"+dossierSauvegarde

    nomDot = dossier + ("/tracking_graph___P%s.dot" % puits)    
    nomTxt =  dossier + ("/tracking_graph___P%s_features.txt" % puits)

    nomFichier = "PL"+plaque+"___P"+puits+"___T00000.xml"
    nomFichier = dossierSauvegarde+"/"+nomFichier

    try:
        (X,Y) = (coord(nomTxt)[0],coord(nomTxt)[1])
        ensemTraj = lireDot(nomDot, X, Y)
        print ensemTraj.numMitoses
        sortiesAll(nomFichier, ensemTraj)
    except IOError as e :
        print "probleme fichier"+e.message


    ecrireClassDef()
    return
