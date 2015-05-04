import os, pdb, sys, time, operator
from matplotlib import pyplot as p

from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.markers import MarkerStyle
from scipy.stats import linregress, pearsonr

import numpy as np

from trajPack import featuresSaved, featuresNumeriques
from PyPack import fHacktrack2
from dataPack.joining import FeatureException
from trajPack import featuresHisto

from util.plots import makeColorRamp, basic_colors, markers, couleurs

import brewer2mpl
from operator import itemgetter

def plotPermutationResults(res):
    f=p.figure(figsize=(24,13))
    for i,el in enumerate(res.keys()):
        ax=f.add_subplot(3,4,i)
        n,bins,patches=ax.hist(res[el], bins=30, normed=False)
        ax.set_title('{} {}'.format(el[0][:9], el[1]))
    p.title("P-values from Fisher's exact test, nb permutations = {}".format(len(res[el])))   
    p.show()

def plotCompaDiv(*args):
    for i,tab in enumerate(args[:-1]):
        t=tab[:-1]
        f=p.figure()
        ax=f.add_subplot(111)
        ax.errorbar(range(4), [np.mean(t[k]) for k in range(4)], [np.std(t[k]) for k in range(4)], fmt='o')
        ax.set_xlim(-0.5,3.5)
        p.title(tab[-1])
        f.savefig(args[-1]+'_{}.png'.format(tab[-1]))
    return

def plotHisto(histogramme, nb_bins_list):
    for feature in featuresHisto:
        lCour=[]; print feature
        for ll in histogramme[feature]:
            lCour.extend(ll)
        lCour=np.array(lCour)
        fig = p.figure()
        p.title(feature)
#        n, bins, patches = ax.hist(x, 50, normed=1, facecolor='green', alpha=0.75)
        for i,el in enumerate(nb_bins_list):
            ax = fig.add_subplot(2,3,i)
            ax.set_title(el)
            n,bins,patches = ax.hist(lCour[np.where(np.isinf(lCour)==False)], el, normed=False)
        p.show()
        
    return

def plotTrajectories(coordL, labels, result_folder):
    resultCour=[]
    how_many =np.bincount(labels)
    clusters = len(how_many)
    for coord in coordL:
        resultCour.extend(coord)

    XX=[]; YY=[]
    for k in range(len(resultCour)):
        X=np.array(resultCour[k][1]); X-=X[0]; XX.append(X)
        Y=np.array( resultCour[k][2]); Y-=Y[0]; YY.append(Y)
    minx,maxx = min([min(X) for X in XX]), max([max(X) for X in XX])
    miny,maxy = min([min(Y) for Y in YY]), max([max(Y) for Y in YY])
    
    assert(len(resultCour)==len(labels))
    
    for cluster in range(clusters):
        for i in np.where(np.array(labels)==cluster)[0][np.random.permutation(how_many[cluster])[:20]]:
            if not os.path.isdir(os.path.join(result_folder, 'images_{}'.format(clusters))): 
                os.mkdir(os.path.join(result_folder, 'images_{}'.format(clusters)))
            plotAlignedTraj(XX[i], YY[i], minx, maxx, miny, maxy, show=False, val=cluster, 
                 name=os.path.join(result_folder, 'images_{}'.format(clusters), 'cluster{}_{}.png'.format(cluster, i)))

    return 1

def plotAlignedTraj(X, Y,minx, maxx, miny, maxy,turn=False, show=False, name=None, val=None):

    X=X-X[0]; Y=Y-Y[0]
    if turn:
        r=linregress(X, Y)
        a=r[0]; b=r[1]
        a = -a
        Xb = (X-a*(Y-b))/float(sqrt(1+a**2))
        Yb=(X*a+Y-b)/float(sqrt(1+a**2))
        X=Xb-Xb[0]; Y=Yb-Yb[0]
    
        if np.mean(np.sign(X))<0:
            X=-X

    f=p.figure()
    ax=f.add_subplot(111); 
    ax.grid(color='b', linewidth=1)
    ax.plot(X, Y, color= basic_colors[1])
    ax.set_xlim((minx,maxx))
    ax.set_ylim((miny, maxy))
    ax.set_xticks(np.arange(minx, maxx+1, 10))
    ax.set_yticks(np.arange(miny, maxy+1, 10))
    if show:
        p.show()
    else: 
        fig = p.gcf()
        fig.set_size_inches(24,13)
        p.title('{} pts, val {}'.format(len(X), val))
        p.savefig(name,dpi=100)
        #f.savefig(name)
        p.close(1)
    return

def plotTraj3d(lstTraj, plate = '?', well='?', N=100, colors=None, name=None):
    #mpl.rcParams['legend.fontsize'] = 10
#    c = []
    coordonnees=[]
    for l in lstTraj:
        l=l.lstPoints
        arr = np.array(sorted(zip(np.array(l.keys())[:,0], np.array(l.values())[:,0], np.array(l.values())[:,1]), key=itemgetter(0)))
        if arr.shape[0]>40:
            coordonnees.append(arr)
    
    fig = p.figure()
    ax = fig.gca(projection='3d')
    for k in range(min(len(coordonnees),N)):
        if colors is not None:
            ax.plot(coordonnees[k][:,0], coordonnees[k][:,2], coordonnees[k][:,1], color=colors[k%len(colors)])
        else:
            ax.plot(coordonnees[k][0], coordonnees[k][1], coordonnees[k][2])
    #p.title(str(N)+' time trajectories from plate '+plate+', well '+well)
    p.ylabel('Frame number')
#pour ne pas afficher les chiffres devant les graduations des axes
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    if name is not None:
        fig.savefig(name)
    else:
        p.show()
    return

def plotRecall():    
    N = 5
    menMeans = (98.9, 88.6, 70.1,57.5 ,48.3)
    
    ind = np.arange(N)  # the x locations for the groups[0, 1, 2, 3, 4]
    width = 0.3       # the width of the bars
    
    fig = p.figure(figsize=(14,8))
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, menMeans, width, color='grey')
    
    #womenMeans = np.array([ 96.84270593,  30.52530384,  80.85932379,  83.76261076,  64.18520351]) premiere fois
    
    #womenMeans = np.array([99.31, 91.30,78.31,83.73,41.26]) #seconde fois m60_80_s100_140
    womenMeans= np.array([ 97.90629431,  61.49068323,  79.51807229,  84.93975904,  63.46456693])# m60_80_s100_120
    #ind2=[0, 1, 2, 3, 4]
    rects2 = ax.bar(ind+width, womenMeans, width, color='blue')
    
     #(99.52,88.64 , 90.54  ,  90.03 , 90.70) chiffres ms lesquels ? regarder ds le cahier probablement dist_inc2 et class_weights
    #me =(99.5911369032, 90.1515151515, 89.1891891892, 88.5630498534, 91.6279069767) #chiffres avec dist_inc3 et class_wieghts avec loss modifiee
    r=np.array([ 99.70028972,  87.42138365,  91.01796407,  86.06060606,  91.25596184])
    
    me =r
    rects3 = ax.bar(ind+2*width, me, width, color='red')
    
    # add some
    ax.set_ylabel('Recall')
    ax.set_title('Recall by event class')
    ax.set_xticks(ind+width)
    ax.set_xticklabels( ('move', 'appear', 'disappear', 'merge', 'split') )
    
    fig.legend( (rects1[0], rects2[0], rects3[0]), ('Cell Cognition', 'Cell Profiler', 'Struct. learning'))
    
    p.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    
    
    p.show()
    
def plotPrecision():
    events=('Move', 'Appear', 'Disappear', 'Merge', 'Split')
    N = 5
#PLOTTING PRECISION
    menMeans = ( 98.0884032 ,  22.7184466 ,  29.17933131,  83.78947368,80.11283498)
    #womenMeans = np.array([ 99.13367628,  99.95802186,  99.35100496,  99.38597358,  96.97432128]) #premiere fois
    #womenMeans= np.array([98.31, 27.85, 66.33, 80.55, 91.74]) #seconde fois m60_80_s100_140
    womenMeans= np.array([ 98.87124163,  55.61797753,  58.66666667,  73.74670185,  58.11088296])#m60_80_s100_120
    #98.31/27.85/66.33/80.55/91.74
    ind = np.arange(N)  # the x locations for the groups[0, 1, 2, 3, 4]
    width = 0.3       # the width of the bars
    
    fig = p.figure(figsize=(14,8))
    ax = fig.add_subplot(121)
    rects1 = ax.bar(ind, menMeans, width, color='grey', alpha=0.8)
    rects2 = ax.bar(ind+width, womenMeans, width, color='blue', alpha=0.8)
    #me = (99.50220874,  90.69767442,  82.20858896,  89.88439306, 93.40909091) chiffres ms lesquels ? regarder ds le cahier probablement dist_inc2 et class_weights
    #me =( 99.42394417,  92.96875   ,  88.,  92.31927711, 94.38116932)#chiffres avec dist_inc3 et class_wieghts avec loss modifiee
     #chiffres avec double cross-validation, class_weights ??
    me=np.array([ 99.41557363,  86.875     ,  81.72043011,  94.44444444,  95.12779553])
    rects3 = ax.bar(ind+2*width, me, width, color='red', alpha=0.8)
    #fig.legend( (rects1[0], rects2[0], rects3[0]), ('CNN', 'LAP', 'MotIW') )
    # add some
    ax.set_ylabel('Precision')
    ax.set_title('Precision by event class')
    ax.set_xticks(ind+width+0.1)
    ax.set_xticklabels(events, fontsize=15)
#PLOTTING RECALL
    menMeans = (98.9, 88.6, 70.1,57.5 ,48.3)
    womenMeans= np.array([ 97.90629431,  61.49068323,  79.51807229,  84.93975904,  63.46456693])# m60_80_s100_120
    me =np.array([ 99.70028972,  87.42138365,  91.01796407,  86.06060606,  91.25596184])
    ax = fig.add_subplot(122)
    rects1 = ax.bar(ind, menMeans, width, color='grey', alpha=0.8)
    rects2 = ax.bar(ind+width, womenMeans, width, color='blue', alpha=0.8)
    rects3 = ax.bar(ind+2*width, me, width, color='red', alpha=0.8)
    
    # add some
    ax.set_ylabel('Recall')
    ax.set_title('Recall by event class')
    ax.set_xticks(ind+width+0.1)
    ax.set_xticklabels( events,  fontsize=15)
    
    fig.legend( (rects1[0], rects2[0], rects3[0]), ('CNN', 'Jaqaman et al.', 'MotIW'),loc=1)
    #p.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    p.show()
    
    #N = 1
    #menMeans = (97.4)
    #
    #ind = np.arange(N)  # the x locations for the groups[0, 1, 2, 3, 4]
    #width = 0.35       # the width of the bars
    #
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #rects1 = ax.bar(ind, menMeans, width-0.05, color='grey')
    #
    #womenMeans = (99.2)
    #rects2 = ax.bar(ind+width, womenMeans, width-0.05, color='blue')
    #
    #me = (99.2)
    #rects3 = ax.bar(ind+2*width, me, width-0.05, color='red')
    #
    ## add some
    #ax.set_ylabel('Total recall')
    #ax.set_xticks(ind+width)
    #ax.set_xticklabels( ('move', 'appear', 'disappear', 'merge', 'split') )
    #
    #ax.legend( (rects1[0], rects2[0],rects3[0]), ('Cell Cognition', 'Lou et Hamprecht', 'us') )
    #
    #plt.show()


def plotComparison(correspF, correspTF, N):
    
    fig=p.figure()   
    num=0; plots=[]; labels=[]; first=True
    fig.suptitle('Comparison between training set and predicted data')#("Plate {0}, well {1}".format(plate, well))
    colors = makeColorRamp(N, basic_colors, hex_output=True)
    calculus={}
    for feat in features1:
        num+=1
        print '------------------', feat, '---------------------'
        calculus.update({feat:[]})
        m=100; M=-100; i=0
        ax=fig.add_subplot(4,4, num) 
        X=[]; Y=[]
        for plate in correspF:
            for well in correspF[plate]:
                print plate, well
                x=[]; y=[]
                cCorrespF = correspF[plate][well]
                cCorrespTF = correspTF[plate][well]
                if len(cCorrespF)!=len(cCorrespTF):
                    raise LengthError
                for k in range(len(cCorrespF)):
                    if cCorrespF[k].features[feat] is None or cCorrespTF[k].features[feat] is None:
                        continue
                    y.append(cCorrespF[k].features[feat])
                    x.append(cCorrespTF[k].features[feat])
                #pdb.set_trace()    
                if first:   
                    plots.append(ax.scatter(x, y, color=colors[i]))
                    labels.append('{0}_{1}'.format(plate, well))
                else:
                    ax.scatter(x, y, color=colors[i])
                #coeff de correlation pour la droite x=y
                X.extend(x); Y.extend(y)
                m=min(m, min(y)); M=max(M, max(y))
                i+=1
        X=np.array(X); Y=np.array(Y)
        coeff_corr=np.mean((X-np.mean(X))*(Y-np.mean(Y)))/float(np.sqrt(np.var(X)*np.var(Y)))
        residues=np.mean((X-Y)**2)
        print "taille de l'echantillon", len(X)
        print 'coeff corr', coeff_corr
        print 'residues', residues
        calculus[feat].extend([len(X), coeff_corr, residues])       
        m=m*1.2; M=M*1.2
        first=False
        #print m, M
        #        pdb.set_trace()
        ax.set_title(feat)
        #ax.set_xlabel('Training set data'); 
        ax.set_ylabel('Predicted data')
        ax.plot([m,M], [m, M], color='black')
        ax.set_ylim(m, M)
    fig.legend(tuple(plots), tuple(labels), 'lower right')
    p.show()
        
    return calculus    


def plotFeaturesTraj(features, tab, N, folder=None, show=False):
    
    fig=p.figure(figsize=(24,13))   
    num=0; plots=[]; labels=[]; first=True
    fig.suptitle('Evolution of features in time')#("Plate {0}, well {1}".format(plate, well))
    colors = makeColorRamp(N, basic_colors, hex_output=True)
    #calculus={}
    for feat in features1:
        num+=1
        print '------------------', feat, '---------------------'
        #calculus.update({feat:[]})
        m=100; M=-100; i=0
        ax=fig.add_subplot(4,4, num) 
        X=[]; Y=[];
        for plate in features:
            for well in features[plate]:
                print plate, well
                x=[]; y=[]
                cFeatures = features[plate][well]
                tabC = tab[plate][well]
                index=featuresListArray.index(feat)
                ff=True
                for k in range(len(cFeatures)):
                    if np.isnan(tabC[k][index]):
                        continue
                    y=tabC[k][index]
                    x=cFeatures[k].beginning
                    
                    if first and ff:   
                        plots.append(ax.scatter(x, y,s=6, color=colors[i]))
                        #ax.text(x, y, '{}_{}'.format(cFeatures[k].numCellule, cFeatures[k].id))
                        labels.append('{0}_{1}'.format(plate, well))
                    else:
                        ax.scatter(x, y,s=2, color=colors[i])
                        #ax.text(x, y, '{}'.format(cFeatures[k].numCellule))
                    ff=False
                #coeff de correlation pour la droite x=y
                    X.append(x); Y.append(y)
        
                i+=1
        m=min(Y); M=max(Y)   
        m=m*1.2; M=M*1.2
        first=False
        #print m, M
        #        pdb.set_trace()
        ax.set_title(feat)
        #ax.set_xlabel('Training set data'); 
        ax.set_ylabel('Predicted data')
        #ax.plot([m,M], [m, M], color='black')
        ax.set_ylim(m, M)
    fig.legend(tuple(plots), tuple(labels), 'lower right')
    if show:
        p.show()
    elif folder is not None:
        fig.savefig(os.path.join(folder, 'featuresTraj_{}.png'.format(well))) 
        
    return 1

def plotBoth(tfeatures, pfeatures, tcoord, pcoord, plate=None, well=None):
    N=50; TroisD=False
    r=tfeatures ; coordC = tcoord
    fig=p.figure()
    ax=fig.add_subplot(121)
    if TroisD:
        ax2=fig.add_subplot(122, projection ='3d')
    else: ax2=fig.add_subplot(122)
    maxD = 10
    maxT = 1.1#2500
    ax.set_title('Displacements'); ax2.set_title('Some features')
    fig.suptitle("Plate {0}, well {1}, red=training set, blue=predictions".format(plate, well))
    #           pdb.set_trace()
    for k in range(min(len(coordC),N)):
        if 'movement type' not in r[k].features: 
            raise FeatureException("Feature 'movement type' is not in r[{}].features, it is not normal".format(k))
        elif r[k].features['movement type'] is not None and r[k].features['diffusion coefficient'] is not None: #and r[k].features['movement type']<0:
            
            ax.plot(coordC[k][0]-coordC[k][0][0], coordC[k][1]-coordC[k][1][0], label=str(k), color = 'red')
            ax.text(coordC[k][0][-1]-coordC[k][0][0], coordC[k][1][-1]-coordC[k][1][0], '%d'%r[k].numCellule, ha='center', va='bottom')
            if TroisD :
                ax2.scatter(r[k].features['diffusion coefficient'], r[k].features['movement type'],r[k].features['direction'], color = '#'+fHacktrack2.couleurs[k])
            else: ax2.scatter(r[k].features['diffusion coefficient'], r[k].features['movement type'], color = 'red')
            maxD = max(maxD, r[k].features['diffusion coefficient']) 
            maxT = max(maxT, r[k].features['movement type']) 
            if TroisD:
                ax2.text(r[k].features['diffusion coefficient'], r[k].features['movement type'], r[k].features['direction'], '%d'%k, zdir=None)
            else:ax2.text(r[k].features['diffusion coefficient'], r[k].features['movement type'], '%d'%r[k].numCellule)
    r=pfeatures; coordC=pcoord
    for k in range(min(len(coordC),N)):
        if 'movement type' not in r[k].features: 
            raise FeatureException("Feature 'movement type' is not in r[{}].features, it is not normal".format(k))
        elif r[k].features['movement type'] is not None and r[k].features['diffusion coefficient'] is not None: #and r[k].features['movement type']<0:
            
            ax.plot(coordC[k][0]-coordC[k][0][0], coordC[k][1]-coordC[k][1][0], color = 'blue')
            ax.text(coordC[k][0][-1]-coordC[k][0][0], coordC[k][1][-1]-coordC[k][1][0], '%d'%r[k].numCellule, ha='center', va='bottom')
            if TroisD :
                ax2.scatter(r[k].features['diffusion coefficient'], r[k].features['movement type'],r[k].features['direction'], color = '#'+fHacktrack2.couleurs[k])
            else: ax2.scatter(r[k].features['diffusion coefficient'], r[k].features['movement type'], label=str(k), color = 'blue')
            maxD = max(maxD, r[k].features['diffusion coefficient']) 
            maxT = max(maxT, r[k].features['movement type']) 
            if TroisD:
                ax2.text(r[k].features['diffusion coefficient'], r[k].features['movement type'], r[k].features['direction'], '%d'%k, zdir=None)
            else:ax2.text(r[k].features['diffusion coefficient'], r[k].features['movement type'], '%d'%r[k].numCellule)
            
    #            ax.xlabel('X'); ax2.xlabel("Movement type")
    #            ax.ylabel('Y'); ax2.ylabel('Diffusion coefficient')
    #pdb.set_trace()
    ax2.set_xlabel('Diffusion coefficient'); ax2.set_ylabel('Movement type')#; ax2.set_zlabel('Mean displacement direction')
    if not TroisD: 
        ax2.set_xlim(None, maxD*1.1)
        ax2.set_ylim(None, maxT*1.1)
    ax.legend()
    ax2.legend()
    p.show()
    
    return

#TroisD variable enables to choose if we wish to get a 3D plot or not
#training says if the data comes from the training set or was predicted

def plot(features, tab, coordonnees, TroisD, folder=None, show=False, training=True, 
         plate=None, well=None, oneByOne =None):
    N=50
    indexMovementType = featuresListArray.index('movement type')
    indexDiff = featuresListArray.index('diffusion coefficient')
    if type(features)==dict:
        for plate in features:
            for well in features[plate]:
                r=features[plate][well] ; coordC = coordonnees[plate][well]; tabC=tab[plate][well]
                fig=p.figure(figsize=(24,13))
                ax=fig.add_subplot(121)
                if TroisD:
                    ax2=fig.add_subplot(122, projection ='3d')
                else: ax2=fig.add_subplot(122)
                maxD = 10
                maxT = 1.1#2500
                ax.set_title('Displacements'); ax2.set_title('Some features')
                fig.suptitle("Plate {0}, well {1}, training={2}".format(plate, well, training))
                rr=range(min(len(coordC),N, len(fHacktrack2.couleurs))) if oneByOne is None else [oneByOne]
                for k in rr:
                    if not np.isnan(tabC[k][indexMovementType]) and not np.isnan(tabC[k][indexDiff]): 
        #plotting 2D trajectories
                        ax.plot(coordC[k][0]-coordC[k][0][0], coordC[k][1]-coordC[k][1][0], label=str(k), color = '#'+fHacktrack2.couleurs[k])
                        ax.text(coordC[k][0][-1]-coordC[k][0][0], coordC[k][1][-1]-coordC[k][1][0], '%d'%r[k].numCellule,
                                ha='center', va='bottom')
        #plotting features either in 3D or 2D
                        if TroisD :
                            ax2.scatter(tabC[k][indexDiff],r[k].beginning,
                                          tabC[k][indexMovementType], color = '#'+fHacktrack2.couleurs[k])
                            
                        else: ax2.scatter(tabC[k][indexDiff], tabC[k][indexMovementType],s=5, color = '#'+fHacktrack2.couleurs[k])
                        maxD = max(maxD, tabC[k][indexDiff]) 
                        maxT = max(maxT, tabC[k][indexMovementType]) 
                        if TroisD:
                            ax2.text(tabC[k][indexDiff], r[k].beginning, tabC[k][indexMovementType], '%d'%r[k].numCellule, zdir=None)
                        else:ax2.text(tabC[k][indexDiff], tabC[k][indexMovementType], '{}_{}'.format(r[k].numCellule, r[k].id))
                        
    #            ax.xlabel('X'); ax2.xlabel("Movement type")
    #            ax.ylabel('Y'); ax2.ylabel('Diffusion coefficient')
                #pdb.set_trace()
                ax2.set_xlabel('Diffusion coefficient'); ax2.set_ylabel('Movement type')#; ax2.set_zlabel('Mean displacement direction')
                if not TroisD: 
                    ax2.set_xlim(None, maxD*1.1)
                    ax2.set_ylim(None, maxT*1.1)
                ax.legend();
                if show:
                    p.show()
                elif folder is not None:
                    fig.savefig(os.path.join(folder, 'coordFeatures_{}.png'.format(well)))
                else:
                    return 0
                return 1
    else:
        r=features ; coordC = coordonnees; tabC=tab
        fig=p.figure(figsize=(24,13))
        ax=fig.add_subplot(121)
        if TroisD:
            ax2=fig.add_subplot(122, projection ='3d')
        else: ax2=fig.add_subplot(122)
        maxD = 10
        maxT = 1.1#2500
        ax.set_title('Displacements'); ax2.set_title('Some features')
        fig.suptitle("Plate {0}, well {1}, training={2}".format(plate, well, training))
        #           pdb.set_trace()
        rr=range(min(len(coordC),N, len(fHacktrack2.couleurs))) if oneByOne is None else [oneByOne]
        for k in rr:
            if not np.isnan(tabC[k][indexMovementType]) and not np.isnan(tabC[k][indexDiff]): 
#plotting 2D trajectories
                ax.plot(coordC[k][0]-coordC[k][0][0], coordC[k][1]-coordC[k][1][0], label=str(k), color = '#'+fHacktrack2.couleurs[k])
                ax.text(coordC[k][0][-1]-coordC[k][0][0], coordC[k][1][-1]-coordC[k][1][0], '%d'%r[k].numCellule,
                        ha='center', va='bottom')
#plotting features either in 3D or 2D
                if TroisD :
                    ax2.scatter(tabC[k][indexDiff],r[k].beginning,
                                  tabC[k][indexMovementType], color = '#'+fHacktrack2.couleurs[k])
                    
                else: ax2.scatter(tabC[k][indexDiff], tabC[k][indexMovementType], color = '#'+fHacktrack2.couleurs[k])
                maxD = max(maxD, tabC[k][indexDiff]) 
                maxT = max(maxT, tabC[k][indexMovementType]) 
                if TroisD:
                    ax2.text(tabC[k][indexDiff], r[k].beginning, tabC[k][indexMovementType], '%d'%r[k].numCellule, zdir=None)
                else:ax2.text(tabC[k][indexDiff], tabC[k][indexMovementType], '{}_{}'.format(r[k].numCellule, r[k].id))
                
#            ax.xlabel('X'); ax2.xlabel("Movement type")
#            ax.ylabel('Y'); ax2.ylabel('Diffusion coefficient')
        #pdb.set_trace()
        ax2.set_xlabel('Diffusion coefficient'); ax2.set_ylabel('Movement type')#; ax2.set_zlabel('Mean displacement direction')
        if not TroisD: 
            ax2.set_xlim(None, maxD*1.1)
            ax2.set_ylim(None, maxT*1.1)
        ax.legend();
        if show:
            p.show()
            return 1
        elif folder is not None:
            fig.savefig(os.path.join(folder, 'coordFeatures_{}.png'.format(well)))
            return 1
        else:
            return 0  
        
        
def plotClusteringIndices(x, *args):
    #exemple d'appel : plots.plotClusteringIndices(range(2,16), [BIC, RBIC, 'BIC'], [PEC, RPEC, 'PEC'], [FS, RFS, 'FS'])
    n=len(args)
    f=p.figure(figsize=(24,13))
    for i, tab in enumerate(args):
        ax=f.add_subplot(1, n, i)
        if type(tab)==list:#means we have both random and real data values
            t=np.array(tab[0]); Rt=np.array(tab[1])
            if t.shape[0]==len(x):
                ax.plot(x, t,  linestyle="dashed", marker="o", color="red", label='real data')

            else:
                ax.plot(x, np.mean(t, 0),  linestyle="dashed", marker="o", color="red", label='real data - 18 features')
                #ax.plot(x, np.mean(t, 0),  linestyle="dashed", marker="o", color="red", label='10 neighbours')
                ax.errorbar(x, np.mean(t, 0), np.std(t, 0), linestyle='None', marker="None", color="red", label='Stdev on {} samples'.format(t.shape[0]))
            
            ax.plot(x, np.mean(Rt, 0),  linestyle="dashed", marker="o", color="black", label='all real data - 13 features')
            #ax.plot(x, np.mean(Rt, 0),  linestyle="dashed", marker="o", color="black", label='random data (uniform)')
            #ax.plot(x, np.mean(Rt, 0),  linestyle="dashed", marker="o", color="black", label='15 neighbours')
            ax.errorbar(x, np.mean(Rt, 0), np.std(Rt, 0), linestyle='None', marker="None", color="black", label='Stdev on {} samples'.format(Rt.shape[0]))
            
            if type(tab[2])!=str:
                Rgt=np.array(tab[2])
                ax.plot(x, np.mean(Rgt, 0),  linestyle="dashed", marker="o", color="green", label='random data (noisy Gaussian)')
                #ax.plot(x, np.mean(Rgt, 0),  linestyle="dashed", marker="o", color="green", label='20 neighbours')
                ax.errorbar(x, np.mean(Rgt, 0), np.std(Rgt, 0), linestyle='None', marker="None", color="green", label='Stdev on {} samples'.format(Rgt.shape[0]))
            
            ax.legend();ax.set_title(tab[-1])
    p.show()
    return

def plotClustInd(folder, x, labels=None,neighbours=None,show=False, *args):
    #exemple d'appel : plots.plotClusteringIndices(range(2,16), [BIC, RBIC, 'BIC'], [PEC, RPEC, 'PEC'], [FS, RFS, 'FS'])
    n=len(args)
    f=p.figure(figsize=(24,13))
    colors = [ 'red','blue', 'grey', 'green', 'yellow']#fHacktrack2.couleurs
    lignes = []
    for i, tab in enumerate(args):
        ax=f.add_subplot(1, n, i)
        ax.grid(True)
        if type(tab)==list:#means we have both random and real data values
            zz=tab
            if type(tab[-1])==str:
                zz=tab[:-1]
            for k, t in enumerate(zz):
                if t is None:
                    print 'tab {} is None !'.format(k)
                    continue
                t=np.array(t)
                lab = '{}'.format(k)
                if neighbours is not None:
                    lab ='{}, neighbours {}'.format(labels[k/len(neighbours)], neighbours[k%len(neighbours)])

                elif labels is not None: #then labels contain the labels
                    lab=labels[k]
                if t.shape[0]==len(x):
                    ax.plot(x, t,  linestyle="dashed", marker=markers[k], color=colors[k], label=lab)
    
                else:
#                    ax.plot(x, np.mean(t, 0),  linestyle="dashed", marker=markers[k], color='#'+colors[k], label=lab)
                    #ax.plot(x, np.mean(t, 0),  linestyle="dashed", marker="o", color="red", label='10 neighbours')
                    ax.errorbar(x, np.mean(t, 0), np.std(t, 0), linestyle="dashed", marker=markers[k], color=colors[k], label=lab+', {} iterations'.format(t.shape[0]))
                from matplotlib.ticker import MultipleLocator
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                ax.xaxis.grid(True,'minor')
                ax.yaxis.grid(True,'minor')
#                ax.plot(x, np.mean(Rt, 0),  linestyle="dashed", marker="o", color="black", label='all real data, only 13 features')
#                #ax.plot(x, np.mean(Rt, 0),  linestyle="dashed", marker="o", color="black", label='random data (uniform)')
#                #ax.plot(x, np.mean(Rt, 0),  linestyle="dashed", marker="o", color="black", label='15 neighbours')
#                ax.errorbar(x, np.mean(Rt, 0), np.std(Rt, 0), linestyle='None', marker="None", color="black", label='Stdev on {} samples'.format(Rt.shape[0]))
#                
#                if type(tab[2])!=str:
#                    Rgt=np.array(tab[2])
#                    ax.plot(x, np.mean(Rgt, 0),  linestyle="dashed", marker="o", color="green", label='random data (noisy Gaussian)')
#                    #ax.plot(x, np.mean(Rgt, 0),  linestyle="dashed", marker="o", color="green", label='20 neighbours')
#                    ax.errorbar(x, np.mean(Rgt, 0), np.std(Rgt, 0), linestyle='None', marker="None", color="green", label='Stdev on {} samples'.format(Rgt.shape[0]))
#                
            if type(tab[-1])==str:
                ax.set_title(tab[-1])
#                 if 'Silhouette' in tab[-1]:
#                     ax.set_ylim(0,0.4)
    legend = ax.legend();
    #ax.set_ylim(0.3,1)
    for label in legend.get_texts():
        label.set_fontsize('small')

    p.savefig(os.path.join(folder, 'outputClustInd{}.png'.format(time.time())))
    if show:
        p.show()
    return

def plotKMeansPerFilm(m, ctrlStatus,labls, colors, Td=False, ctrlcov=None, ctrlcolors=None, phenocolors=None, mp=None):
#    who=labels#donc la en cliquant sur un point on aura le nom du gene. Sinon il faut donner le nom de l'experience dans les labels
    global labels
    labels=labls
#    print labels.shape
    f=p.figure()
    ax=f.add_subplot(111, projection ='3d') if Td is True else f.add_subplot(111)
    if not Td and ctrlcov is not None:
        X, Y = np.mgrid[0:100:10000j, 0:100:10000j]
        zz = np.c_[X.ravel(), Y.ravel()]
        robdist = ctrlcov.mahalanobis(zz)
        robdist = robdist.reshape(X.shape)
        ax.contour(X,Y, np.sqrt(robdist), cmap=p.cm.PuBu_r, linestyles='dashed')
        ax.imshow(np.rot90(np.sqrt(robdist)), cmap=p.cm.gist_earth_r,extent=[0, 100, 0, 100])
#dealing with colors
    cols = []; xbL=[]
    if colors=='gene':
        for i,el in enumerate(labls):
            if el=='no target':cols.append("black")
            elif el=='ctrl': cols.append("red")
            else: 
                if 'ctrl' in labls:
                    cols.append("green")
                else:
                    if el.split(' ')[0] not in xbL:
                        xbL.append(el.split(' ')[0])
                    cols.append(couleurs[xbL.index(el.split(' ')[0])])
    elif colors=='ctrlStatus':
        for i, el in enumerate(ctrlStatus):
            if el==0:
                cols.append('red')
            elif el==1:
                cols.append("green")
            else: cols.append("orange")
            
    elif type(colors)==dict or type(colors)==list:
        cols=colors
    else:
        raise
#pour que le picker marche en fait il faut faire ax.scatter(x,y, picker=True)
    if ctrlcolors !=None:
        mc=m[np.where(np.array(ctrlStatus)==0)]
        ax.scatter(mc[:,0], mc[:,1], mc[:,2], color=ctrlcolors, picker=True)
        
        ax.scatter(mp[:,0], mp[:,1], mp[:,2], color=phenocolors, picker=True)
    else:
        if Td is False:
            ax.scatter(m[:,0], m[:,1], color=cols, picker=True)
    #            if labels[k]!='ctrl':
    #                ax.text(m[k,0], m[k,1], labels[k], ha='center', va='bottom')
        else:
            ax.scatter(m[:,0], m[:,1], m[:,2], color=cols, picker=True)
            for i, label in enumerate(labls):
                ax.text(m[i,0], m[i,1], m[i,2], label)
            
    f.canvas.mpl_connect('pick_event', onpick)
    p.show()

def onpick(event):
    global labels
    ind = event.ind
    for k in range(len(ind)):
        print ind[k], labels[ind[k]]

def plotMovies(folder, arr, filename = 'patternMovies'):
    nbPheno = arr.shape[0]
    X,Y=np.meshgrid(range(arr.shape[1]+1),range(arr.shape[2]+1))
    print 'plotting movies'
    for k in range(nbPheno):
        p.figure(figsize=(12,12))
        cmap = brewer2mpl.get_map('RdBu', 'diverging', 3).mpl_colormap
        im = p.pcolormesh(X,Y,arr[k].transpose(), cmap=cmap)
        p.colorbar(im, orientation='horizontal')
        p.show()
        #p.savefig(os.path.join(folder, '{}_{}.png'.format(filename, k)))
    return 1

def featureFeatureCorrelation(data, feature_names = featuresNumeriques, sh=True, method = 'ward', metric='euclidean'):
    '''
    This function plots the correlation between features
    '''
    correlation_mat = np.zeros(shape=(len(feature_names), len(feature_names)))
    for i in range(len(feature_names)):
        correlation_mat[i,i]=1
        for j in range(i+1, len(feature_names)):
            correlation_mat[i,j]=pearsonr(data[:,i], data[:,j])[0]
            correlation_mat[j,i]=correlation_mat[i,j]
    if sh:
        heatmap(np.absolute(correlation_mat), feature_names, feature_names, normalization =False, row_method=method,
            column_method=method, row_metric=metric, column_metric=metric, color_gradient = 'red_black_sky', filename = 'correlation_features')
    else:
        return correlation_mat
    

def featureTimeCorrelation(time_length, data, feature_names, sh=True):
    ''' 
    This function plots the correlation of given features in list feature_names, with temporal trajectory length
    '''
    f=p.figure(figsize=(24,13))
    for i,feature in enumerate(feature_names):
        ax=f.add_subplot(4,4,i)
        zou=zip(time_length, data[:,featuresSaved.index(feature)])
        zou=np.array(sorted(zou, key=operator.itemgetter(0)))
        ax.scatter(zou[:,0], zou[:,1])
        ax.set_title('{} \n corr {}'.format(feature, pearsonr(zou[:,0], zou[:,1])[0]))
        
    if sh: p.show()
    else:
        p.savefig('../TIMEFEATURES.png')
    return

def featureLogTransformation(vector, vectorC, vectorP, name, nbins=500):
    '''
    This function compares a feature distribution with and without logarithm transformation by plotting them
    '''
    f=p.figure(figsize=(24,13))
    ax=f.add_subplot(221); histogramme(ax, vector, name, nbins)
    ax=f.add_subplot(222); histogramme(ax, vectorP, name, nbins, None, vectorC)
    
    ax=f.add_subplot(223); histogramme(ax, vector, name, nbins, True)
    ax=f.add_subplot(224); histogramme(ax, vectorP, name, nbins, True, vectorC)
    
    #p.title(name+', green = pheno, red = ctrl')
    p.savefig('../resultData/dist_feature_logOuPas_{}.png'.format(name))

def histogramme(ax, vector1, name, nbins =500, transfo=None, vector2=None):
    if transfo is not None:
        print "transforming by putting distribution min to 1"
        vector1=vector1-np.min(vector1)+1
        vector1=np.log(vector1)
        name+=', log transformed'
    name+=', nbins = {}'.format(nbins)
    if vector2==None:
        print 'all'
        n, bins, patches = ax.hist(vector1,nbins, normed=1, histtype='bar', cumulative=False)
        p.setp(patches, 'facecolor', 'b', 'alpha', 0.5)
        ax.set_title(name)
#        p.show()
    else:
        print "pheno ctrl"
        if transfo is not None:
            vector2=vector2-np.min(vector2)+1
            vector2=np.log(vector2)
        n, bins, patches = ax.hist(vector1,bins=np.linspace(0, max(np.max(vector1),np.max(vector2)),nbins), normed=1, histtype='bar', cumulative=False,label="Ctrl", color='green', alpha=0.5)
        #p.setp(patches, 'facecolor', 'g', 'alpha', 0.5)
        n2, bins2, patches2 = ax.hist(vector2,bins=np.linspace(0, max(np.max(vector1),np.max(vector2)),nbins), normed=1, histtype='bar', cumulative=False,label='Exp', color='red', alpha=0.5)
        #p.setp(patches2, 'facecolor', 'r', 'alpha', 0.5)
        ax.set_title(name); ax.legend()
#        p.show()


class LengthError(Exception):
    def __init__(self):
        print "Length error between objects", sys.exc_info()
        pass
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
