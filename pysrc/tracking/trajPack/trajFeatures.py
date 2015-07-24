import pdb, os, sys, time
import numpy as np
import cPickle as pickle
#import matplotlib.pyplot as plt

from math import sqrt, log, fabs
from optparse import OptionParser
from warnings import warn
from scipy import odr
from scipy.stats import linregress
from scipy.spatial import cKDTree, ConvexHull

from tracking.plots import plot, plotComparison, plotBoth, plotFeaturesTraj
from tracking.PyPack import fHacktrack2
from tracking.trajPack import moments, rayon, densities, windows,\
    TIMEDEPTH, featuresHisto, featuresSaved
    
from util.plots import couleurs, makeColorRamp

count = 0
monOutput = ""
loadingFolder = "../results/"
movie_length={}

#appC :  0 for cells present at the beginning of the movie
#        1 for appearing cells
#        2 for cells coming from division

#dispC : 0 for cells present at the end of the movie
        #1 for disappearing cells
        #2 for dividing cells
        #3 for merging cells
        
def ordonnancement(new_sol):
    stories = {}
    for sol, d in new_sol:#attention new_sol est une liste de solutions pas forcement dans l'ordre de leurs plates wells etc
        if sol.plate not in stories:
            stories[sol.plate]= {}
        if sol.well not in stories[sol.plate]:
            stories[sol.plate][sol.well]= {}
        if sol.index not in stories[sol.plate][sol.well]:
            stories[sol.plate][sol.well].update({sol.index:sol});
    return stories

def trackletBuilder(new_sol, centresFolder, training=False):
    '''
    J'ajoute un compteur HoleError dans le cas ou Cplex ne respecte pas les contraintes : des cellules ne sont pas atteintes
    et par consequent elles ont l'air de surgir de nulle part sans etre apparues via l'algorithme. Je considere que dix erreurs de ce type sont acceptables
    par movie
    
    '''
    nb_voisins_max=10
    #distance_densi=40
    stories = ordonnancement(new_sol); res = {}
    global movie_length
    connexions = {}
    for plate in stories:
        res.update({plate : {}})
        connexions.update({plate : {}})
        movie_length.update({plate : {}})
        
        for well in stories[plate]:
            
            count_HoleError=0
            
            print plate, well
            ensembleTraj = fHacktrack2.ensTraj()
            res[plate].update({well:ensembleTraj})
            connexions[plate].update({well:{}})
            num = 0; l=[0]; trajC=None
            idd=0; firstIndex=min(stories[plate][well].keys()); 

            movie_length[plate].update({well : max(stories[plate][well].keys())-firstIndex+2})
            #+2 because there are no tracking in the last frame but it exists
            f=open(os.path.join(centresFolder, 'centers_P{}_w{}.pkl'.format(plate, well)), 'r')
            lCentersDict=pickle.load(f); f.close()
            
            for index in stories[plate][well]:
                print " -----INDEX ", index,
                dC={}
                lCenters =lCentersDict[index]
                nlCenters =lCentersDict[index+1]
                treeC = cKDTree(lCenters, leafsize = 10)      
                nTree = cKDTree(nlCenters, leafsize = 10)          
                
                connexions[plate][well].update({index:dC})
                solC = stories[plate][well][index]
                if type(solC.truth)==str:
                    l.append(index+1)
                    continue
                singlets = solC.singlets
                nextSinglets = solC.nextSinglets
                try:
                    result =solC.result
                except AttributeError:
                    result = filter(lambda x: x[1]==1, zip(solC.hypotheses, solC.truth))
                
                for hyp in result:
                    s, t = hyp[0]; source = [] ; target = []; 
                    #print hyp
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
                    mer = (incoming>1)
#                    if mit or mer:
#                        print'#####', hyp
                    if incoming == 0:
                        #print "APPEAR"
                        for c in target:
                            newTraj= fHacktrack2.trajectoire(num, c.center[0], c.center[1], index+1, c.label, idd)
                            newTraj.addDensity(c.center, nTree, nb_voisins_max, densities)
                            #print num, idd
                            ensembleTraj.add(newTraj)
                            num+=1; idd+=1
                            
                        continue
                    elif outcoming ==0:
#                        print "DISAPPEAR"
                        if index in l:
                            #pdb.set_trace()#en fait ce sont au moins les trajectoires qui commencent en 0 et finissent en 0 ie les apparitions disparitions en 0
                            for c in source:
                                newTraj = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label, idd)
                                #ici c'est une disparititon de trajectoire qui dure une image donc inutile d'ajouter la densite
                                ensembleTraj.add(newTraj)
                                num+=1; idd+=1
                        continue
                    elif incoming == 1: 
                #FROM ONE
                        c = source[0]
                        incomingTraj = ensembleTraj.findLabelFrame(c.label, index)
                        if incomingTraj==[]:
                            #It can happen only at the first frame : trajectories have not been initiated yet
                            #otherwise I raise HoleError if it has been observed more than 10 times
                            if index>firstIndex and training==False: 
                                count_HoleError+=1
                                if count_HoleError>=10:
                                    raise HoleError(plate, well, index)
                                
                            if not mit:
                                t=target[0]
                                newTraj = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label, idd)
                                newTraj.addDensity(c.center, treeC, nb_voisins_max, densities)
                                newTraj.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
                                newTraj.addDensity(t.center, nTree, nb_voisins_max, densities)
                                ensembleTraj.add(newTraj)
                                #print num, idd
                                num+=1; idd+=1
                                continue
                            elif mit:
                                newTraj = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label, idd)
                                newTraj.addDensity(c.center, treeC, nb_voisins_max, densities)
                                idSource=(idd,)
                                ensembleTraj.add(newTraj)
                                num+=1; numC=num; idd+=1; targetCourantes=[]; 
                                for t in target:
                                    newTraj = fHacktrack2.trajectoire(numC, t.center[0], t.center[1], index+1, t.label, idd)
                                    newTraj.addDensity(t.center, nTree, nb_voisins_max,densities)
                                    targetCourantes.append(idd)
                                    #print numC, idd
                                    idd+=1
                                    ensembleTraj.add(newTraj)
                                dC[idSource]=tuple(targetCourantes)

                        elif not mit:
                            first = True
                            traj=incomingTraj[0]
                            t = target.pop()
                            traj.ajoutPoint(t.center[0], t.center[1], index+1, t.label)
                            traj.addDensity(t.center, nTree, nb_voisins_max, densities)
                            #print traj.numCellule, traj.id
                        elif mit:
#                            print 'mit'
                            first = True; fusion=False; idCourantes=[]; targetCourantes=[]
                            for i, traj in enumerate(incomingTraj):
                                idCourantes.append(traj.id)
                                if first: 
                                    numC = incomingTraj[i].numCellule
                                first=False
                                fusion=fusion or traj.fusion
                                incomingTraj[i].numCellule = numC
                                
                            for t in target:
                                newTraj = fHacktrack2.trajectoire(numC, t.center[0], t.center[1], index+1, t.label, idd)
                                newTraj.addDensity(t.center, nTree, nb_voisins_max, densities)
                                newTraj.fusion = fusion
                                ensembleTraj.add(newTraj)
                                targetCourantes.append(idd)
                                #print numC, idd
                                idd+=1
                            dC[tuple(idCourantes)]=tuple(targetCourantes)
    
                    elif mer:
                        if mit:
                            raise
#                        print "FROM MORE THAN ONE"
                        t = target[0]
                        trajC = fHacktrack2.trajectoire(num, t.center[0], t.center[1], index+1, t.label, idd)
                        trajC.addDensity(t.center, nTree, nb_voisins_max, densities)
                        trajC.fusion=True
                        idTarget=(idd,)
                        #print num, idd
                        idd+=1; num+=1; ensembleTraj.add(trajC)

                        idSource=[]#; tt=[]
                        for c in source:
                            incomingTraj = ensembleTraj.findLabelFrame(c.label, index)
                            try:
                                #tt.append(incomingTraj[0])
                                idSource.append(incomingTraj[0].id)
#                                print incomingTraj[0].numCellule, incomingTraj[0].id
                            except IndexError:
                                #donc il faut creer des nouvelles trajectoires
                                if index>firstIndex and training==False: 
                                    count_HoleError+=1
                                    if count_HoleError>=10:
                                        raise HoleError(plate, well, index)
                                    
                                trajC = fHacktrack2.trajectoire(num, c.center[0], c.center[1], index, c.label, idd)
                                trajC.addDensity(c.center, treeC, nb_voisins_max, densities)
                                idSource.append(idd)
                                idd+=1; num+=1
                                ensembleTraj.add(trajC)
                                #pdb.set_trace()#dans quel cas peut-on bien tomber ici ??
                                #on tombe ici a t=0 mais est-ce qu'on y tombe apres ? oui si on calcule les trajectoires du training
                                #set pcq il y a des evenements elimines par defaut (mvt a quatre)
                        dC[tuple(idSource)]=idTarget
#                print dC

    return res, connexions, movie_length

def movementType(mom, length, verbose):
    '''
    Here we compute the movement type as indicated in this publication: Sbalzarini 2005a
    If the correlation coefficient for the coefficient diffusion regression is low, we put 0 instead.
    If one of the regressions leading to the movement type has a correlation coefficient that is under 
    0.7, then we add 1 to indicateur. It's because we are interested to know which trajectories have
    this issue, and if it's only for one regression or for all.
    
    Returns: diffusion coefficient, movement type and adequateness of diffusion
    '''
    
    indicateur = 0
    x=np.log(range(1, int(length/3)+1));y=np.log(mom)
    gamma=np.zeros(shape=(len(moments),))
    for nu in moments:
        r=linregress(x, y[nu-1])
        if verbose>5:
            print "correlation coefficient", r[2]
            print 'p-value', r[3]
        gamma[nu-1]=r[0]
        if nu==2:
            if r[2]>=0.70:
                D=np.exp(r[1])/(4)
            else:
                D=0
            corr = r[2]

        if r[2]<0.70:
            indicateur += 1
    r2=linregress(moments, gamma)

    if verbose>5:
        print "slope ", r2[0], "diffusion coefficient", D 
    return indicateur, r2[0], D, corr

def localStraightness(vecX, vecY, sh=False):
    debut = time.clock()
    initValues=linregress(vecX, vecY)
    
    def f(B, x):
        '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
        return B[0]*x + B[1]
    
    linear = odr.Model(f)
    data = odr.RealData(vecX, vecY)
    myOdr = odr.ODR(data, linear, beta0=(initValues[0], initValues[1]))
    myOutput = myOdr.run()
    temps= time.clock()-debut
    
    if sh:
        result = myOutput.beta
        import matplotlib.pyplot as p
        f=p.figure()
        ax=f.add_subplot(111)
        ax.scatter(vecX, vecY)
        ax.plot(vecX, vecX*result[0]+result[1])
        ax.set_title('Sum of squared distances to function {}'.format(myOutput.sum_square))
        p.savefig('../resultData/images/points_ODR_{}.png'.format(len(vecX)))
    return temps, myOutput.sum_square

def convexHullArea(X, Y, v):
#v for verbose mode
    def aireTriangle(a, b, c):
        return 0.5*np.absolute((b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]))
    
    result = 0
    points = np.vstack((X,Y)).transpose()
    cvxHull = ConvexHull(points)
    seenSimplices =[]
    
    a=cvxHull.simplices[0][0]
    b=cvxHull.simplices[0][1]
    seenSimplices.append(list(cvxHull.simplices[0]))
    while len(seenSimplices)<len(cvxHull.simplices):
        s =filter(lambda x: b in x and list(x) not in seenSimplices, cvxHull.simplices)[0]
        #print s, seenSimplices
        c=s[0] if np.any(s[0]!=b) else s[1] 
        if v>5:
            print a,b,c
        if a!=c: 
            r=aireTriangle(points[a], points[b], points[c])
            result+=r
        seenSimplices.append(list(s)); b=c
        
#    import matplotlib.pyplot as plt
#    for i, p in enumerate(points):  
#        plt.plot(p[0], p[1], 'o')
#        plt.text(p[0], p[1], str(i))
#    for simplex in cvxHull.simplices:
#        plt.plot(points[simplex,0], points[simplex,1], 'k-')
#    plt.show()
    return result

#def plotAllBoules(b):
#    colors = makeColorRamp(len(rayon), basic_colors)
#    f=plt.figure(figsize=(10,13))
#    for index in range(len(rayon)):
#        boules=b[index]
#        boules.plot(f, colors[index], index)
#    plt.show()
#    return

def sq_norm(X, Y):
    result=0
    for k in range(len(X)):
        try:
            dT = (X[k+1]-X[k], Y[k+1]-Y[k])
        except IndexError:
            pass
        else:
            result+=dT[0]**2+dT[1]**2
            
    return result
    

def filter_mitosis_related_moves(t,X,Y, a,d):
    if a==2:
        t=t[2:]; X=X[2:]; Y=Y[2:]
    if d in [2,3]:
        t=t[:-2]; X=X[:-2]; Y=Y[:-2]
    return t,X,Y

def simpleCorrection(listTraj, movie_start, average):
    newcoord=[]
    for traj in listTraj:
        l= sorted(traj.lstPoints.keys())
        debut = min(l)[0]; fin = max(l)[0]
        t=[k[0] for k in l]
        X=np.array([traj.lstPoints[k][0] for k in l])
        Y = np.array([traj.lstPoints[k][1] for k in l])
    
        #X=X-average[debut:fin+1, 0]; Y=Y-average[debut:fin+1,1]
        if debut ==movie_start:
            X[1:] = X[1:]-average[:fin-movie_start,0]
            Y[1:] = Y[1:]-average[:fin-movie_start,1]
        else:
            X = X-average[debut-1-movie_start:fin-movie_start,0]
            Y = Y-average[debut-1-movie_start:fin-movie_start,1]
        newcoord.append([t,X,Y])
    return newcoord

def computingHisto(traj, average, movie_start, verbose, a,d, training):
#FEATURES
    #Working with tracklets ie no split or merge in the lstPoints
    
    l= sorted(traj.lstPoints.keys())
    debut = min(l)[0]; fin = max(l)[0]
    X=np.array([traj.lstPoints[k][0] for k in l])
    Y = np.array([traj.lstPoints[k][1] for k in l])
    t=[el[0] for el in l]
    
#also keeping raw coordinates in case i want to do smt with cell visualization in cell cognition - but also deleting mitosis related moves to be able to link this info
#with features
    raw_t = np.array(l)
    raw_X = np.array(X)
    raw_Y = np.array(Y)
    raw_t, raw_X, raw_Y = filter_mitosis_related_moves(raw_t, raw_X, raw_Y,a,d)
    
#i.removing average displacement in movie to correct for plate movement
    if average is not None:
        if debut ==movie_start:
            X[1:] = X[1:]-average[:fin-movie_start,0]
            Y[1:] = Y[1:]-average[:fin-movie_start,1]
        else:
            X = X-average[debut-1-movie_start:fin-movie_start,0]
            Y = Y-average[debut-1-movie_start:fin-movie_start,1]
    else:
        if not training:
            raise AttributeError
    
#ii. keeping every coordinate above 0 so that the calculus for turning angle is ok
    X-=X.min()
    Y-=Y.min()

#iii. i remove the points that are <= 2 frames from merge or split    
    t,X,Y = filter_mitosis_related_moves(t,X,Y,a,d)
    debut=min(t); fin = max(t)
    
    compteur=np.zeros(shape=(len(moments),int(len(t)/3)), dtype=float); p=True
    count=np.zeros(shape=(TIMEDEPTH,))

    histN={nom:[] for nom in featuresHisto} #featuresHisto=['acceleration', 'displacements', 'persistence', 'straight', 'turningangle']
#Remarque pour la persistence : pour l'instant on est au niveau un #[[] for ll in range(TIMEDEPTH)]

    r={'time length':len(t),\
       'density':traj.density,\
       'total space length' : 0,\
       'moments':np.zeros(shape=(len(moments), int(len(t)/3)), dtype=float), 
       'turning angle':np.zeros(shape=(5,), dtype=float),\
       'signed turning angle':np.zeros(shape=(TIMEDEPTH,), dtype=float),\
       'effective speed':0,\
       'persistence':np.zeros(shape=(TIMEDEPTH,), dtype=float),\
       'slope logdist-logt':0,\
       'intersection logdist-logt':0,\
       'largest move':0,\
       'convex hull area': 0,\
       'entropy':np.zeros(shape=(len(rayon),), dtype=float),
       'ball number':np.zeros(shape=(len(rayon),), dtype=float)
       }
    angles=[]
#BE CAREFUL that try/except ZeroDivisionError don't work here I don't know why
    for k in range(len(X)):
        x=X[k]; y=Y[k]
#FEATURE RSS of line closest to points in time windows
        for timeWindow in windows:
            try:
                vecX = X[k:k+timeWindow]; vecY = Y[k:k+timeWindow]
            except IndexError:
                pass
            else:
                if len(vecX)==timeWindow:
                    temps, localStr= localStraightness(vecX, vecY)
                    if not np.isnan(localStr):
                        histN["straight"].append(localStr/sq_norm(vecX, vecY))
                    if verbose>1:
                        print 'TIME', temps, localStr
        try:
            dT = (X[k+1]-x, Y[k+1]-y)
        except IndexError:
            continue
        else:
            n=np.linalg.norm(dT)
#FEATURE largest move
            r['largest move']=max(n, r['largest move'])
            histN['displacements'].append(n)     
                
            try:
                histN['acceleration'].append(histN['displacements'][-1]-histN['displacements'][-2])
            except IndexError:
                pass
                
#FEATURE SIGNED TURNING ANGLE
            angles.append(np.arctan2(dT[1], dT[0]))
            
            for j in range(1,TIMEDEPTH+1):
                try:
                    dTplus = (X[k+j+1]-X[k+j],Y[k+j+1]-Y[k+j])
                except IndexError:
                    pass
                else:
                    n2=np.linalg.norm(dTplus)
                    count[j-1]+=1
                    if n!=0 and n2!=0:
                        if np.dot(dT, dTplus)/(n*n2)>1:
                            vv=0
                        elif np.dot(dT, dTplus)/(n*n2)<-1:
                            vv=np.pi
                        else:
                            vv=np.arccos(np.dot(dT, dTplus)/(n*n2))
                        persistence = np.dot(dT, dTplus)/(n*n2)
                    else:
                        vv=0
                        persistence=1
#FEATURE turning angle
                    histN['turningangle'].append(vv)
#FEATURE persistence by measuring direction correlation
                    histN['persistence'].append(persistence)
       
#FEATURE : moment of order nu at interval delta
        for nu in moments:
            for delta in range(1, int(len(t)/3)+1):
                try:
                    v=((x-X[k+delta])**2+(y-Y[k+delta])**2)**(float(nu)/2)
                except IndexError:
                    continue
                else:
                    compteur[nu-1][delta-1]+=1  
                    r['moments'][nu-1][delta-1]+=v

#FEATURE entropy in circles of radius epsilon, divided by maximum entropy (log (length traj))
    eps=4; msg=''
    for i,ll in enumerate(rayon):
        bb=Boules(ll+eps, X, Y)
        r['entropy'][i]=bb.entropy()
        r['ball number'][i]=len(bb.lstBoules)/float(sqrt(r['time length']))
        if verbose>5:
            msg+='radius {}, entropy {}'.format(ll, r['entropy'][i])
#UNNORMALIZED FEATURES
    #FEATURE : trajectory length 
    r['total space length']=r['moments'][0][0]
    #FEATURE  effective distance
    p=0; P=len(t)-1; tdT = (X[P]-X[p], Y[P]-Y[p])
    r['effective space length']=np.linalg.norm(tdT)
    #FEATURE largest move
    #FEATURE convex hull area
    
    r['convex hull area'] = convexHullArea(X, Y, verbose)
    
#NORMALIZED FEATURES
    r['moments']=r['moments']/compteur
    
    #FEATURE mean squared displacement
    r['mean squared displacement']=r['moments'][1][0]
    
    #FEATURE convex hull area divided by time length
    r['norm convex hull area'] = r['convex hull area']/float(r['time length']-1)
    
    #FEATURE signed turning angle
    r['signed turning angle']=fabs(np.arctan2(np.sum([np.sin(angles[k+1]-angles[k]) for k in range(len(angles)-1)]),
                    np.sum([np.cos(angles[k+1]-angles[k]) for k in range(len(angles)-1)])))
    
    #FEATURE : corrected straightness index
    r['corrected straightness index']=sqrt(r['time length'])*r['effective space length']/float(r['total space length'])
    
    #FEATURE sinuosity
#NON CETTE FEATURE N'EST QUE DU BRUIT

    #FEATURE movement type, diffusion coefficient
    r['mvt type adequation'], r['movement type'], r['diffusion coefficient'], r['diffusion adequation'] = movementType(r['moments'], len(t), verbose)        

    #FEATURE effective speed    
    r['effective speed']=r['effective space length']/float(sqrt(r['time length']-1))
    
#    print "finally correlation between instantaneous speed and persistence:"
#    
    if len(histN['persistence'])+1!=len(histN['displacements']):
        raise
    else:
        slope, _,corr, pValue, _ = linregress(histN['displacements'][:-1], histN['persistence'])
        r['correlation speed and persistence'] = corr
        r['slope SP'] = slope
    
    if verbose>5:
        print "tracklet length", r['time length']
        print msg
        print 'corrected straightness index', r['corrected straightness index'] 
        print 'total space length', r['total space length']
        print 'corr coeff', corr, 'p-value', pValue
 
   
    features=['entropy', 'ball number']
    for feature in features:
        for k in range(len(rayon)):
            nom = feature+str(k+1)
            r[nom]=r[feature][k]
        del r[feature]
        
    for feature in featuresHisto:
        nom = 'mean '+feature
        r[nom]=np.mean(histN[feature])
        
    r['FLATmoments']=list(r['moments'].flatten())
    del r['moments'] #ATTENTION QUE SUIVANT LA DUREE DE LA TRAJECTOIRE LES GENS N'ONT PAS LE MEME NB DE PTS DS LES MOMENTS

    return r, [t, X, Y],[raw_t, raw_X, raw_Y],histN

def driftCorrection(m_length, outputFolder, plate, well, totalTracklets, just_average_speed=False):
    #TO CORRECT for drift we take into account all displacements, even those in tracks with length smaller than 10
    filename = "avcum_P{}_{}.pkl".format(plate, well)
    movie_start = min([min(sorted(traj.lstPoints.keys())) for traj in totalTracklets])[0]
    allDisp = np.zeros(shape=(m_length, len(totalTracklets), 2), dtype=float)
    trajPerFrame = np.zeros((m_length,))
    #here we are working with tracklets ie no split or merge in the lstPoints
    i=0
    for traj in totalTracklets:
        l= sorted(traj.lstPoints.keys())
        for k in l:
            x=traj.lstPoints[k][0]; y=traj.lstPoints[k][1]
            try:
                dT = (traj.lstPoints[l[l.index(k)+1]][0]-x, traj.lstPoints[l[l.index(k)+1]][1]-y)
            except IndexError:
                continue
            else:
                allDisp[k[0]-movie_start][i]=dT
                trajPerFrame[k[0]-movie_start]+=1
        i+=1
    #average = np.zeros(shape=(m_length, 2), dtype=float)
    trajPerFrame=np.array((trajPerFrame, trajPerFrame), dtype=float).transpose()
    trajPerFrame[trajPerFrame==0]=0.00001
    
    #pdb.set_trace()
    if just_average_speed:
        return allDisp, trajPerFrame
    
    average = np.sum(allDisp, 1)/trajPerFrame
    avcum = np.cumsum(average, 0)
    
    f=open(os.path.join(outputFolder, filename), 'w')
    pickle.dump(avcum, f); f.close()
    return avcum, movie_start #taille movie length, 2

def app(track, connexions, m_start):
    #0 if the cell is there at the begining
    #2 if it comes from a division
    #1 if it just appeared
    start =min(track.lstPoints.keys())[0]
    if start == m_start:
        return 0
    elif start > m_start:
        try:
            lC = connexions[start-1]
        except KeyError:
            return 1
        else:
            if np.any([track.id in vv for vv in lC.values()]):#your code could hardly be more obscure, good job!!
                return 2
            return 1

def disp(track, connexions, m_end):
    end = max(track.lstPoints.keys())[0]
    #0 if the cell is there at the end
    #2 if it's a merge
    #3 if it's a split
    #1 if it just disappears
    if end == m_end:
        return 0
    elif end<m_end:
        try:
            lC = connexions[end]
        except KeyError:
            return 1
        else:
            try:
                who = filter(lambda x: track.id in x, lC.keys())[0]
                if len(who)==1:
                    return 2
                elif len(who)>1:
                    return 3
            except:
                return 1
            
def computingDensity(track):
    r=np.zeros(shape=(len(densities),1))
    try:
        for i in range(len(densities)):
            r[i]=[np.max(track.density[:,i])] #np.mean(track.density[:,i]),np.std(track.density[:,i]), np.min(track.density[:,i]), np.max(track.density[:,i])]
    except AttributeError:
        print "You are working with trajectories that do not have density attribute. Pgm exiting"
        sys.exit()
    
    return r

def histogramPreparationFromTracklets(dicT, connexions, outputFolder, training, verbose, movie_length, name, filtering_fusion=True,
                                      time_window=None, time_lapse=4):
    print 'histogram version'
    coord=[]; tabF={}

    if filtering_fusion:
        track_filter= (lambda x: x.fusion !=True and len(x.lstPoints)>11)
    else:
        warn("Taking into account all tracklets, even those resulting from a fusion")
        track_filter= (lambda x: len(x.lstPoints)>11)

    for plate in dicT:
        if plate not in tabF:
            #r[plate]={}
            tabF[plate]={}
        for well in dicT[plate]:
            if well not in tabF[plate]:
#                r[plate][well]=[]
                tabF[plate][well]=None
            print plate, well
            dicC = dicT[plate][well]
            movie_start = movie_length[plate][well]
            histNC={nom:[] for nom in featuresHisto}
#            histAvSC={nom:[] for nom in featuresHisto}
#            histAvMC={nom:[] for nom in featuresHisto}
#            histApSC={nom:[] for nom in featuresHisto}
            #going into all trajectories to get small tracklets
            #THEN computing all displacement vectors
            #THEN averaging that to prevent prbl from mvts of the microscope, table etc
            print 'First pass of all tracks to correct drift'

            if not training:  
                average, movie_start = driftCorrection(movie_length[plate][well]-1, outputFolder, plate, well, dicC.lstTraj)
            else:
                average=None
                movie_start=0

#ATTENTION here we throw away tracklets whose length is smaller than 10 (otherwise we don't have diffusion info because
#I use delta t in range (1, len(track)/3 +1) ie (1,2,3) if len(track)>9, and then do regression on it,
#so if we have less than three points to do linear regression it seems a bit dubious to me)

            if time_window is not None:
                print "Filtering on time window ", time_window, "for a time-lapse of ", time_lapse, "frames per hour"
                movie_end = min(time_window[1]*time_lapse, movie_length[plate][well] + movie_start-1)
                movie_start = max(time_window[0]*time_lapse, movie_start)
                print "Hence starting at ", movie_start, "and ending at ", movie_end
            else:
                movie_end = movie_length[plate][well] + movie_start-1
            print 'Second pass of all tracks to extract features of tracks longer than 11 frames'

            tabFeatures = None;
            for track in filter(track_filter,dicC.lstTraj):
                t=track.lstPoints.keys(); t.sort(); 
                if time_window is not None and (movie_start>min(t)[0] or movie_end<max(t)[0]+1):
                    track=track.copyBetween(beginning=movie_start, end=movie_end)
                    if not track_filter(track):
                        continue
                
                if verbose>5:
                    print "debut de l'extrait de trajectoire", min(t)[0], "fin ", max(t)[0]
                
        #completion info sur les trajectoires :
#            #app et disp
                if not training:
                    a, d = app(track, connexions[plate][well], movie_start), disp(track, connexions[plate][well], movie_end)
                else:
                    a=1; d=1
                if a==None or d==None:
                    raise AttributeError
                
                f, coordC,rawCoordC, histN=computingHisto(track, average, movie_start, verbose, a, d, training)
                arr=[f[k] for k in featuresSaved]
                
            #densities
                trackDensity = computingDensity(track).flatten() if not training else [4]
                arr.extend(trackDensity)
                arr=np.array(arr, dtype=float)
                coord.append(rawCoordC)
                tabFeatures = arr if tabFeatures==None else np.vstack((tabFeatures, arr))
                for nom in histNC:
                    histNC[nom].append(histN[nom])
#remarque sur la fonction histogramme de numpy : si on met le keyword normed a True cela normalise la probability density function mais pas les bins, ie on n'a pas un
#histogramme au sens de la these de Sven Siggelkow
                    #r[plate][well].append(featuredTraj(track.id, labelsSequence, min(t)[0], max(t)[0], track.numCellule, features=None))
    
            f=open(os.path.join(outputFolder, name), 'w')
            pickle.dump([tabFeatures, coord, histNC], f)
            f.close()    
            tabF[plate][well]=tabFeatures
            
    return tabF, coord
            

def WhoIsWhen(features, plate, well):
    
    global movie_length
    length = movie_length[plate][well]; fC=features[plate][well]
    r = [[] for x in range(length)]
    for track in fC:
        for l in r[track.beginning: track.end]:
            l.append(track.id)
    
    return r
    

def output(training_only=True, predict=False, new=True):
    #listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
    first = True; 
    loadingFolder = "../prediction/"

    tSol=tracking.gettingSolu(loadingFolder, training_only)
    first=False
    new_tSol = []
    for tSolu in tSol.lstSolutions: 
        new_tSol.append((tSolu, tSolu.truthVec()))
    #pdb.set_trace()
    print "Building trajectories for training set data"
    (tDicTraj, tconn, movie_length) = trackletBuilder(new_tSol, loadingFolder, training=True)
    del new_tSol;
    if predict:
        tSol=tracking.gettingSolu(loadingFolder, False, first)
        new_sol = tracking.sousProcessClassify(tSol, loadingFolder)
        print "Building trajectories for predicted data"
        (dicTraj, conn, movie_length) =trackletBuilder(new_sol, loadingFolder, training=False) 
##        del new_sol;del tSol;
    else:
        dicTraj=None; conn=None
    return dicTraj, conn, tDicTraj, tconn, movie_length


class featuredTraj():
    
    def __init__(self, ide, labelsSequence, beginning, end, num, features, to=None, fr=None):
        self.id = ide
        self.labelsSeq = labelsSequence
        self.beginning = beginning
        self.end = end
        self.numCellule = num
        self.features = features
        self.to = to
        self.fr=fr
        
    def labelsD(self):
        return dict(zip(range(self.beginning, self.end+1), self.labelsSeq))
        
    def findCorresp(self, trackletsL):
#        length=self.end-self.beginning+1
        
        r=dict(zip([track.id for track in trackletsL], [0 for x in range(len(trackletsL))] ))
        featuresR=dict(zip([track.id for track in trackletsL], [0 for x in range(len(trackletsL))] ))
        #pdb.set_trace()
        labelsD =self.labelsD() 
        for tracklet in trackletsL:
            lD = tracklet.labelsD()
            featuresR[tracklet.id]=tracklet
            for frame in labelsD:
                if frame in lD and lD[frame]== labelsD[frame]:
                    r[tracklet.id]+=1
   
        #pdb.set_trace()
        m = max(r.values())
        print "max", m
        #if m>0:pdb.set_trace()
        r=zip(r.keys(), r.values()); #featuresR=zip(featuresR.keys(), featuresR.values())
        result = np.array(filter(lambda x: x[1]==m, r))
        if np.all(result[:,1]==0):#len(correspondances)==1:
            return None, 0, 0
        elif len(result)==1:
            return result[0][0], m, featuresR[result[0][0]]
        else:
            pdb.set_trace()
            return None, 0, 0

def consecutiveList(liste):
    liste.sort()
    currR=[liste[0]]; r=[]

    for n in liste[1:]:
        if n==currR[-1]+1:
            currR.append(n)
        else:
            r=currR if len(currR)>len(r) else r
            currR=[n]
    r=currR if len(currR)>len(r) else r
    return r


        
class Boule():
    
    def __init__(self, eps, center, indexCenter=None):
        self.radius = eps
        self.center= center
        self.count = 1 if indexCenter is not None else 0
        
    def app(self, indexPoint):
        if type(indexPoint)==int:
            self.count+=1
        else:
            self.count+=len(indexPoint)   
            
class Boules():
    
    def __init__(self, eps, X, Y):
        self.radius=eps
        self.points = np.vstack((X, Y)).transpose()
        self.nbpoints=len(self.points)
        points=list(self.points)
        self.lstBoules=[]#Boule(eps, self.points[0], 0)]
        takenCareOf = [False for x in range(len(points))]
        #treeC = cKDTree([self.points[0]], leafsize=10)
#So, I decide : for every point, if it is closer than eps to the last bowl it goes in it
#if not : if it is closer than eps to any bowl it goes in it OR NOT
#en fait je pense que l'on peut calculer les deux types de descripteurs. Dans un cas on decrit si une trajectoire revient bcp sur elle-meme,
#ds l'autre juste si elle fait bcp de pauses
        j=0
        while np.bincount(takenCareOf)[0]>1:
            nb=[]
            tree = cKDTree(points, leafsize=10)
            for p in points:
                d, i = tree.query(p, k=len(points), distance_upper_bound=self.radius)
                d=list(d); i=list(i)

                if type(d)!=float:
                    l=filter(lambda x: d[i.index(x)]<100, i)
                    nb.append(consecutiveList(l))
                elif d<100:
                    nb.append([i])
                else:
                    nb.append([])
    #            who.append(filter(lambda x: d[i.index[x]]<100, i))
            nb=filter(lambda x: takenCareOf[nb.index(x)]==False,nb)
            new = [points[k] for k in range(len(points)) if takenCareOf[k]==False]
            ind = np.argmax([len(ll) for ll in nb])
            b=Boule(eps, new[ind])
            index = nb[ind]
            b.app(index) ;self.app(b)
            
            for k in index:
                points[k]=((j+1)*100000, (j+1)*100000)
                takenCareOf[k]=True
                j+=1
        try:
            last = takenCareOf.index(False)
        except ValueError:
            return
        else:
            b=Boule(eps, points[last], 0); self.app(b)
            return
                
#                
#    def plot(self, f, col, num):
#        ax=f.add_subplot(5,2, num) 
#        X = self.points[:,0]-self.points[0][0]
#        Y= self.points[:,1]-self.points[0][1]
#        ax.plot(X, Y)
#        for boule in self.lstBoules:
#            c = boule.center - self.points[0]
#            circle=plt.Circle(c,self.radius, color=col, fill=False)
#            ax.add_artist(circle)
#        ax.set_xlim((min(-40, min(X)),max(40, max(X))))
#        ax.set_ylim((min(-40, min(Y)),max(40, max(Y))))

            
        
    def app(self, boule):
        self.lstBoules.append(boule)
        
    def lstCenters(self):
        arr=self.lstBoules[0].center
        for k in range(1,len(self.lstBoules)):
            arr=np.vstack((arr, self.lstBoules[k].center))
        return arr
    
    def entropy(self):
        nbTotalPoints=float(self.nbpoints)
        e=0
        for b in self.lstBoules:
            e+=-b.count/nbTotalPoints*log(b.count/nbTotalPoints)
        e=e/log(nbTotalPoints)
        return e
        
class CorrelationError(Exception):
    def __init__(self):
        self.msg = "Not used"
        pass

class HoleError(Exception):
    def __init__(self, plate, well, frame):
        print 'Pbl of missing predictions for frame in the middle of the video : \n'
        print 'plate: {}, well: {}, frame: {}'.format(plate, well, frame)
        pass

if __name__ == '__main__':
    
    description =\
'''
%prog - Cell tracking
'''
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)

    parser.add_option("-p", "--predict", dest="predict", default = False, help="Predicting or not")

    (options, args) = parser.parse_args()
    
    training=False
    predict = int(options.predict)
    TroisD = False
    outputFolder = '../trajectories/'
    d, c, td, tc=output(True, predict); 
    verbose=4
    print movie_length
    #FOR PREDICTED DATA
    if predict:
        features, tabFeatures, coordonnees = featuresPreparationFromTracklets(d, c, outputFolder, False, verbose, tab=True)
    
    #plot(features, tabFeatures, coordonnees, TroisD, folder=None, show=True, training=False)
        #plot(features, coordonnees, TroisD, training=True, plate=None, well=None, oneByOne =None):
    #FOR TRAINING SET DATA
    tFeatures, ttabFeatures, tCoordonnees = featuresPreparationFromTracklets(td, tc, outputFolder, True,verbose, tab=True)
    #plot(tFeatures, ttabFeatures, tCoordonnees, TroisD, folder=None, show=True)
    N=5
    if not predict:
        plotFeaturesTraj(tFeatures,ttabFeatures, N, folder=None, show=True)
    else:
        plotFeaturesTraj(features, tabFeatures, N, folder=None, show=True)
    #plot(features, tabFeatures, coordonnees, TroisD, folder=None, show=True, training=False)
    #here I consider every predicted tracklet
    correspTF={}; correspF={}; N=0
    vec_persistence=[]; vec_speed=[]       
    
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
            print "nb de traj a trouver", len(tfC)
#from the point of view of predicted data
#first i see which predicted tracklet is present when
            r=WhoIsWhen(tFeatures, plate, well)
            corresp = {}
            toPlotF = []; toPlotCoord=[]
            k=0
            for tracklet in fC:
                vec_persistence.append(tracklet.features['persistence'])
                vec_speed.append(tracklet.features['mean instantaneous speed'])
                rC = r[tracklet.beginning: tracklet.end]; zz=[]
                for l in rC: zz.extend(l)
                h=np.bincount(zz)
                
                correspondance, nbCommonFrames, feat=tracklet.findCorresp(filter(lambda x: x.id in np.where(h==max(h))[0], tfC)) #filter renvoie un tuple avec une liste comme seul element
                #print correspondance, nbCommonFrames, tracklet.end-tracklet.beginning+1
                
                if correspondance is not None:
                    corresp.update({tracklet.id:correspondance})
                    correspTF[plate][well].append(feat)
                    toPlotF.append(fC[k]); toPlotCoord.append(coordonnees[plate][well][k])
                    #pdb.set_trace()
                k+=1
            correspF[plate][well]=list(toPlotF)
            N+=1
        #plotting on the same plot both training set data and predicted data
            plotBoth(tfC, toPlotF, tCoordonnees[plate][well], toPlotCoord, plate, well)
            plot(toPlotF, toPlotCoord, TroisD, training=False, plate=plate, well=well)

#from the point of view of training data
#first i see which predicted tracklet is present when
#    for plate in features:
#        for well in features[plate]:
#            fC = features[plate][well]; 
#            tfC = tFeatures[plate][well]
#            
#            r=WhoIsWhen(features, plate, well)
#            corresp = {}
#            toPlotF = []; toPlotCoord=[]
#            k=0
#            for tracklet in tfC:
#                rC = r[tracklet.beginning: tracklet.end]; zz=[]
#                for l in rC: zz.extend(l)
#                h=np.bincount(zz)
#                correspondance, nbCommonFrames, feat=tracklet.findCorresp(filter(lambda x: x.id in np.where(h==max(h))[0], fC)) #filter renvoie un tuple avec une liste comme seul element
#                #print correspondance, nbCommonFrames, tracklet.end-tracklet.beginning+1
#                
#                if correspondance is not None:
#                    corresp.update({tracklet.id:correspondance})
#                    correspTF[plate][well].append(feat)
#                    toPlotF.append(tfC[k]); toPlotCoord.append(tCoordonnees[plate][well][k])
#                else:
#                    pdb.set_trace()
#                    #pdb.set_trace()
#                k+=1
#            correspF[plate][well]=list(toPlotF)
#            N+=1
#
#            #pdb.set_trace()
#            plotBoth(fC, toPlotF, coordonnees[plate][well], toPlotCoord, plate, well)
#            plot(toPlotF, toPlotCoord, TroisD, training=False, plate=plate, well=well)
    calculus = plotComparison(correspF,correspTF, N)
    out=open("speed_persistence.pkl", "w")
    pickle.dump([vec_persistence, vec_speed], out)
    out.close()
    out=open("plotComparOutput.pkl", 'w')
    pickle.dump(calculus, out)
    out.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
