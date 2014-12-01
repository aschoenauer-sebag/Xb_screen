import os, sys, pdb, string
import cPickle as pickle
#import prettyplotlib as ppl
import brewer2mpl
import matplotlib.pyplot as p
import numpy as np
import vigra.impex as vi
from warnings import warn
from itertools import repeat
from matplotlib.backends.backend_pdf import PdfPages

from analyzer import cecog_features, CONTROLS
#from interface import fromPlateToPython
from tracking import test
from tracking.trajPack import featuresSaved
from tracking.importPack import FEATURE_NUMBER
from tracking.plots import markers
from tracking.trajPack.trajFeatures import driftCorrection

from plateSetting import readPlateSetting
from util.plots import plotBarPlot, couleurs, linestyles
from util.listFileManagement import expSi, siEntrez, correct_from_Nan
from util.sandbox import histLogTrsforming
from tracking.PyPack.fHacktrack2 import sortiesAll, initXml, ecrireXml, finirXml

def plot_special_case(ll):
    timep=(0,40,96,144,184)
    timep1=(0,40,96,144)

    if ll[2][0].shape[1]<140:
        print "Modifying shapes"
        ll[0][0]=np.hstack((ll[0][0][:,:77], np.zeros(shape=(ll[0][0].shape[0], 9)), ll[0][0][:,77:]))
        ll[0][1]=np.hstack((ll[0][1][:,:77], np.zeros(shape=(ll[0][1].shape[0], 9)), ll[0][1][:,77:]))
        
        ll[2][0]=np.hstack((ll[2][0][:,:77], np.zeros(shape=(ll[2][0].shape[0], 9)), ll[2][0][:,77:]))
        ll[2][1]=np.hstack((ll[2][1][:,:77], np.zeros(shape=(ll[2][1].shape[0], 9)), ll[2][1][:,77:]))

    f=p.figure()
    ax=f.add_subplot(111)
    for k in range(5):
        ax.scatter(range(4), ll[0][0][k,(0,40,96,144)], color = couleurs[0], label='130814')#label='050914__1')
    for k in range(9):
        ax.scatter(np.array(range(4))+(k+1)*0.05, ll[0][1][k,timep1], color = couleurs[k+1])
#    for k in range(3):
#        ax.scatter(range(5), ll[1][0][k,timep], color = couleurs[0], label='071114__1')
#    for i in range(7):
#        for k in range(3):
#            ax.scatter(np.array(range(5))+(i+1)*0.05, ll[1][i+1][k,timep], color = couleurs[i+1])
#    for i in range(2):
#        for k in range(2):
#            ax.scatter(np.array(range(5))+(8+i)*0.05, ll[1][8+i][k,timep], color = couleurs[8+i])
    for k in range(5):
        ax.scatter(range(5), ll[1][0][k,timep], color = couleurs[0],  label='150814')
    for k in range(9):    
        ax.scatter(np.array(range(5))+(k+1)*0.05, ll[1][1][k,timep], color = couleurs[k+1])
    for k in range(6):    
        ax.scatter(range(4), ll[2][0][k,timep1], color = couleurs[0], marker='+', label='050914__1')
    for k in range(9):    
        ax.scatter(np.array(range(4))+(k+1)*0.05, ll[2][1][k,timep1], color = couleurs[k+1], marker='+')
    for k in range(3):
        ax.scatter(range(5), ll[3][0][k,timep], color = couleurs[0], marker='+', label='071114__1')
    for i in range(7):
        for k in range(3):
            ax.scatter(np.array(range(5))+(i+1)*0.05, ll[3][i+1][k,timep], color = couleurs[i+1], marker='+')
    for i in range(2):
        for k in range(2):
            ax.scatter(np.array(range(5))+(8+i)*0.05, ll[3][8+i][k,timep], color = couleurs[8+i], marker='+')
            
    ax.legend(); ax.grid(True)
    p.show()


def comparDistributions(filename, labels,values, features=None, 
                        outputFolder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions',
                        mode='distribution', normed=True):
    if features==None:
        features=featuresSaved
        
    if mode=='distribution':
        for k in range(min(len(features),19)):
            f,axes=p.subplots(1,len(values), figsize=(26,8), sharex=True, sharey=True)
            n,b,pat=axes[0].hist(values[0][0][:,featuresSaved.index(features[k])], bins=50, color=couleurs[0], normed=normed, alpha=0.5,label=labels[0])
            for i, tab in enumerate(values):
                if i>0:
                    n,b,pat=axes[i].hist(tab[0][:,featuresSaved.index(features[k])], bins=b, color=couleurs[i], alpha=1,normed=normed, label=labels[i])
                    n,b,pat=axes[i].hist(values[0][0][:,featuresSaved.index(features[k])], bins=50, color=couleurs[0], alpha=0.5,normed=normed,label=labels[0])
                    axes[i].grid(True)
                    axes[i].legend()
                    axes[i].set_title('{}+-{}'.format(round(np.mean(tab[0][:,featuresSaved.index(features[k])]), 2), round(np.std(tab[0][:,featuresSaved.index(features[k])]),2)))
            axes[0].set_title('{}+-{}'.format(round(np.mean(values[0][0][:,featuresSaved.index(features[k])]), 2), round(np.std(values[0][0][:,featuresSaved.index(features[k])]),2)))
            p.savefig(os.path.join(outputFolder, '{}_{}.png'.format(filename, features[k])))
    elif mode=='scatter':
        for k in range(min(len(features),19)):
            f,axes=p.subplots(1,len(values), figsize=(26,8), sharex=True, sharey=True)
            tableC=np.vstack((values[0][-1], values[0][0][:,featuresSaved.index(features[k])])).T
            tableC=tableC[np.argsort(tableC[:,0])]
            axes[0].scatter(tableC[:,0], tableC[:,1], color=couleurs[0],label=labels[0])
            for i, tab in enumerate(values):
                if i>0:
                    axes[i].scatter(tableC[:,0], tableC[:,1], color=couleurs[0],label=labels[0])
                    table=np.vstack((tab[-1], tab[0][:,featuresSaved.index(features[k])])).T
                    table=table[np.argsort(table[:,0])]
                    axes[i].scatter(table[:,0], table[:,1], color=couleurs[i], alpha=1, label=labels[i])
                    axes[i].grid(True)
                    axes[i].legend()
                    axes[i].set_title('{}+-{}'.format(round(np.mean(tab[0][:,featuresSaved.index(features[k])]), 2), round(np.std(tab[0][:,featuresSaved.index(features[k])]),2)))
            axes[0].set_title('{}+-{}'.format(round(np.mean(values[0][0][:,featuresSaved.index(features[k])]), 2), round(np.std(values[0][0][:,featuresSaved.index(features[k])]),2)))
            p.savefig(os.path.join(outputFolder, 'Time_{}_{}.png'.format(filename, features[k])))
    else:
        range_={'diffusion coefficient':[0,4], 'mean squared displacement':[0,4],
                'effective speed':[0,10], 'mean displacements':[0,6], 'norm convex hull area':[0,50]}
        features=range_.keys()
        for k in range(min(len(features),19)):
            f,axes=p.subplots(1,len(values), figsize=(26,8), sharex=True, sharey=True)
            tableC=np.vstack((values[0][-1], values[0][0][:,featuresSaved.index(features[k])])).T
            #tableC=tableC[np.argsort(tableC[:,0])]
            c,x,y,im=axes[0].hist2d(tableC[:,0], tableC[:,1], range=[[0,190],range_[features[k]]],bins=50, normed=normed)#, color=couleurs[0],label=labels[0])
            for i, tab in enumerate(values):
                if i>0:
                    #axes[i].scatter(tableC[:,0], tableC[:,1], color=couleurs[0],label=labels[0])
                    table=np.vstack((tab[-1], tab[0][:,featuresSaved.index(features[k])])).T
                    #table=table[np.argsort(table[:,0])]
                    axes[i].hist2d(table[:,0], table[:,1],bins=[x,y], normed=normed) #color=couleurs[i], alpha=1, label=labels[i])
                    axes[i].grid(True)
                    axes[i].legend()
                    axes[i].set_title('{}+-{}'.format(round(np.mean(tab[0][:,featuresSaved.index(features[k])]), 2), round(np.std(tab[0][:,featuresSaved.index(features[k])]),2)))
            axes[0].set_title('{}+-{}'.format(round(np.mean(values[0][0][:,featuresSaved.index(features[k])]), 2), round(np.std(values[0][0][:,featuresSaved.index(features[k])]),2)))
            p.savefig(os.path.join(outputFolder, 'Time_{}_{}.png'.format(filename, features[k])))
    

def xbConcatenation(folder, exp_list, xb_list='processedDictResult_P{}.pkl', qc=None, filename = 'hist_tabFeatures_{}.pkl', verbose=0, perMovie = False):
    who=[]; length=[]; r=[]; X=[]; Y=[]; ctrlStatus = []; xb=[]
    time_length=[]
    
    processed={}

    if qc is not None:
        yqualDict=expSi(qc)

    for i, exp in enumerate(exp_list):
        print i,
        pl,w=exp
        if pl not in processed:
            f=open(os.path.join(folder, xb_list.format(pl)))
            processed[pl]=pickle.load(f); f.close()
    #TO ADD WHEN THERE IS A QC
#        if pl[:9]+'--'+w[2:5] not in yqualDict:
##i. checking if quality control passed
#            sys.stderr.write("Quality control not passed {} {} \n".format(pl[:9], w[2:5]))
#            continue   
        try:
#iii.loading data
            f=open(os.path.join(folder, 'track_predictions', pl, filename.format(w)), 'r')
            arr, coord, _= pickle.load(f)
            f.close()
        except IOError:
            sys.stderr.write("No file {}\n".format(os.path.join(pl, filename.format(w))))
        except EOFError:
            sys.stderr.write("EOFError with file {}\n".format(os.path.join(pl, filename.format(w))))
        else:
            if arr==None:
                sys.stderr.write( "Array {} is None\n".format(os.path.join(pl, filename.format(w))))
                continue
            elif len(arr.shape)==1 or arr.shape[0]<20:
                sys.stderr.write("Array {} has less than 20 trajectories. One needs to investigate why. \n".format(os.path.join(pl, filename.format(w))))
                continue
            else:
                try:
                    arr, toDel = correct_from_Nan(arr, perMovie)

                    r= arr if r==[] else np.vstack((r, arr))
#                    if hist:
#                        for nom in featuresHisto:
#                            histN[nom]=np.delete(np.array(histN[nom]), toDel)
#                            histNtot[nom].extend(list(histN[nom]))
                    ll=arr.shape[0]
    
                except (TypeError, ValueError, AttributeError):
                    sys.stderr.write("Probleme avec le fichier {}".format(os.path.join(pl, filename.format(w))))
                else:   
                    time_length.extend([coord[k][0][len(coord[k][0])/2] for k in filter(lambda x: x not in toDel, range(len(coord)))])
                    xbCourant = processed[pl][int(w.split('_')[0])]['Xenobiotic']
                    xb.append(xbCourant)
                    who.append((pl, w))
                    length.append(ll)
                    
                    if xbCourant in CONTROLS.values():
                        ctrlStatus.append(0)
                    else:
                        ctrlStatus.append(1)

    if r ==[]:
        return None

#log trsforming data
    r2 = histLogTrsforming(r, verbose=verbose)        

    warn('The data was not normalized. Please check that it will be done before applying any algorithm.')
    
    return r2[:,:-1], who,ctrlStatus, length, xb, time_length

def plottingBar(folder, plate, xenobioticL, solvent,sh=True, target='speed', average_over_xb=False):
    timepoints={}
    timepoints['050914']=(0,40,78,136)
    timepoints['071114']=(0,48,96,144,184)
    timepoints['150814']=(0,40,96,144,184)
    timepoints['130814']=(0,40,96,144)
    timecorres=[0,'10h', '24h', '36h','48h']
    
    saved=[]
    
    timepoints=timepoints[plate]; timecorres=timecorres[:len(timepoints)]
    
    if type(xenobioticL)==str:
        OLD_plottingBar(folder, plate, [xenobioticL], solvent,sh, target)
    else:
        average_solv, proliferation_solv = featureEvolOverTime(None, None, folder, True, 192, solvent,True, True, False, [plate])
        
        if target=='proliferation':
            average_solv=proliferation_solv
        elif target !='speed':
            raise AttributeError
        
        if plate !='130814':
            av = np.mean(average_solv, 0)[np.newaxis, :]
            std=np.std(average_solv, 0)[:,np.newaxis][timepoints,:]
            saved.append(average_solv)
            
        for xenobiotic in xenobioticL:
            average_xb, proliferation_xb = featureEvolOverTime(None, None, folder, True, 192, xenobiotic,True, True, False, [plate])
            saved.append(average_xb)
            if target=='proliferation':
                average_xb=proliferation_xb
            elif target !='speed':
                raise AttributeError
        
            if plate == '130814':
                arr=np.vstack((a[:149] for a in average_solv))
                print arr.shape
                saved=[arr]
                av = np.mean(arr, 0)[np.newaxis, :]
                
                average_xb=np.vstack((a[:149] for a in average_xb))
                saved.append(average_xb)
                std=np.std(arr, 0)[:,np.newaxis][(0,40,96,144),:]
                std=[std[:,0].T]
                
                
                if xenobiotic == 'TGF beta 1':
                    average_xb = np.mean(average_xb,0)[np.newaxis, :]
                    std.append(np.vstack((np.std(average_xb, 0)[:,np.newaxis][(0,40,96,144),:], np.zeros(shape=(1,1)))))
                else:
                    std.extend([None for el in range(10)])
                av = np.vstack((av, average_xb))
                av = np.hstack((av, np.zeros(shape=(av.shape[0],100))))
                
            elif plate=='150814':
                if xenobiotic == 'TGF beta 1':
                    average_xb = np.mean(average_xb,0)[np.newaxis, :]
                    std.append(np.std(average_xb, 0)[:,np.newaxis][(0,40,96,144, 184),:])
                else:
                    std=[std[:,0].T]
                    std.extend([None for el in range(10)])
        
                av=np.vstack((av, average_xb))
            else:
                if average_over_xb:
                    av=np.vstack((av, np.mean(average_xb,0)[np.newaxis, :]))
                    std=np.hstack((std, np.std(average_xb, 0)[:,np.newaxis][timepoints,:]))
                else:
                    av=np.vstack((av, average_xb))
                    zz=np.empty(shape=(len(timepoints),average_xb.shape[0])); zz.fill(None)
                    std=np.hstack((std, zz))
    if plate in ['071114', '050914']:
        std=std.T
    av=av[:,timepoints]
                
    f=open(os.path.join(folder, '{}_data{}_{}.pkl'.format(target, plate, xenobiotic)), 'w')
    pickle.dump(saved,f); f.close()
    plotBarPlot(av, '{} average per well'.format(target), ['D {}'.format(k) for k in range(av.shape[0])], timecorres, 
                '{} average over time, {}, {}'.format(target, xenobiotic, plate), target=target,
                stds = std, name='{}_{}'.format(plate,xenobiotic), folder=folder, sh=sh)


def OLD_plottingBar(folder, plate, xenobioticL, solvent,sh=True, target='speed', average_over_xb=False):
    timepoints={}
    timepoints['050914']=(0,40,78,136)
    timepoints['071114']=(0,40,96,144,184)
    timepoints['150814']=(0,40,96,144,184)
    timepoints['130814']=(0,40,96,144)
    timecorres=[0,'10h', '24h', '36h','48h']
    
    saved=[]
    
    timepoints=timepoints[plate]; timecorres=timecorres[:len(timepoints)]
    
    if type(xenobioticL)==str:
        OLD_plottingBar(folder, plate, [xenobioticL], solvent,sh, target)
    else:
        average_solv, proliferation_solv = featureEvolOverTime(None, None, folder, True, 192, solvent,True, True, False, [plate])
        
        if target=='proliferation':
            average_solv=proliferation_solv
        elif target !='speed':
            raise AttributeError
        
        if plate !='130814':
            av = np.mean(average_solv, 0)[np.newaxis, :]
            std=np.std(average_solv, 0)[:,np.newaxis][timepoints,:]
            saved.append(average_solv)
            
        for xenobiotic in xenobioticL:
            average_xb, proliferation_xb = featureEvolOverTime(None, None, folder, True, 192, xenobiotic,True, True, False, [plate])
            saved.append(average_xb)
            if target=='proliferation':
                average_xb=proliferation_xb
            elif target !='speed':
                raise AttributeError
        
            if plate == '130814':
                arr=np.vstack((a[:149] for a in average_solv))
                print arr.shape
                saved=[arr]
                av = np.mean(arr, 0)[np.newaxis, :]
                
                average_xb=np.vstack((a[:149] for a in average_xb))
                saved.append(average_xb)
                std=np.std(arr, 0)[:,np.newaxis][(0,40,96,144),:]
                std=[std[:,0].T]
                
                
                if xenobiotic == 'TGF beta 1':
                    average_xb = np.mean(average_xb,0)[np.newaxis, :]
                    std.append(np.vstack((np.std(average_xb, 0)[:,np.newaxis][(0,40,96,144),:], np.zeros(shape=(1,1)))))
                else:
                    std.extend([None for el in range(10)])
                av = np.vstack((av, average_xb))
                av = np.hstack((av, np.zeros(shape=(av.shape[0],100))))
                
            elif plate=='150814':
                if xenobiotic == 'TGF beta 1':
                    average_xb = np.mean(average_xb,0)[np.newaxis, :]
                    std.append(np.std(average_xb, 0)[:,np.newaxis][(0,40,96,144, 184),:])
                else:
                    std=[std[:,0].T]
                    std.extend([None for el in range(10)])
        
                av=np.vstack((av, average_xb))
            else:
                if average_over_xb:
                    av=np.vstack((av, np.mean(average_xb,0)[np.newaxis, :]))
                    std=np.hstack((std, np.std(average_xb, 0)[:,np.newaxis][timepoints,:]))
                else:
                    av=np.vstack((av, average_xb))
                    zz=np.empty(shape=(len(timepoints),average_xb.shape[0])); zz.fill(None)
                    std=np.hstack((std, zz))
    if plate in ['071114', '050914']:
        std=std.T
    av=av[:,timepoints]
                
    f=open(os.path.join(folder, '{}_data{}_{}.pkl'.format(target, plate, xenobiotic)), 'w')
    pickle.dump(saved,f); f.close()
    plotBarPlot(av, '{} average per well'.format(target), ['D {}'.format(k) for k in range(av.shape[0])], timecorres, 
                '{} average over time, {}, {}'.format(target, xenobiotic, plate), target=target,
                stds = std, name='{}_{}'.format(plate,xenobiotic), folder=folder, sh=sh)
        
def plotSpeedHistogramEvolution(plate, well, speed, name, outputFolder):
    if not os.path.isdir(os.path.join(outputFolder, name)):
        os.mkdir(os.path.join(outputFolder,name))
    for frame in range(speed.shape[0]-6):
        f=p.figure()
        ax=f.add_subplot(111)
        arr = speed[frame:frame+6].flatten()
        n,bins,patches = ax.hist(arr, bins=50,normed=False, range=(0.5,30))
        ax.set_xlim(0,30)
        ax.set_ylim(0,75)
        f.savefig(os.path.join(outputFolder, name, 'hist_speed_{}_{:>05}.png'.format(well,frame)))

    return

def writeTracks(plate, inputFolder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions__settings2', outputFolder='tracking_annotations/annotations', all_=False):
    '''
    To produce annotated tracklets, either all tracklets (option all_=True) or only those used for feature extraction (option all=False)
    '''
    
    if 'tracking_annotations' not in os.listdir(inputFolder):
        os.mkdir(os.path.join(inputFolder, 'tracking_annotations'))
    if 'annotations' not in os.listdir(os.path.join(inputFolder, 'tracking_annotations')):
        os.mkdir(os.path.join(os.path.join(inputFolder, 'tracking_annotations', 'annotations')))
    
    if all_:
        well_files = filter(lambda x: 'traj_noF_densities' in x, os.listdir(os.path.join(inputFolder, plate)))
        rang=3
    else:
        well_files = filter(lambda x: 'hist_tabFeatures' in x, os.listdir(os.path.join(inputFolder, plate)))
        rang=2
    for file_ in well_files:
        well = '{:>05}_01'.format(int(file_.split('_')[rang][1:]))
        nomFichier = "PL"+plate+"___P"+well+"___T00001.xml"
        nomFichier = os.path.join(inputFolder, outputFolder, nomFichier)

        if all_:
            f=open(os.path.join(inputFolder, plate, file_), 'r')
            d=pickle.load(f); f.close()
            ensTraj = d[plate][well]
            coord = sortiesAll(nomFichier, ensTraj)
        else:
            f=open(os.path.join(inputFolder, plate, file_), 'r')
            _,coord,_=pickle.load(f); f.close()
            
            compteur =1
            fichierX = open(nomFichier, "w")
            fichierX.write(initXml())
            
            for lstPoints in coord:
                frameL=lstPoints[0]; X=lstPoints[1]; Y=lstPoints[2]
                d={(frameL[k],0):(X[k], Y[k]) for k in range(len(frameL))}
                txt, coord = ecrireXml(compteur,d, True)
                fichierX.write(txt)
                compteur+=1
                if compteur>1600:
                    break
        
            fichierX.write(finirXml())
            fichierX.close()
            print 'finally ', compteur, 'trajectories'
        
    return

def featureEvolOverTime(dicT, outputFolder, verbose, movie_length,xenobiotic, plot=False,
                        load=False, plates=None):
    tabF={}
    filename='speedOverTime{}.pkl'.format(xenobiotic)
    seen_plates=[]
    
    MIN_LENGTH_MOVIES=np.min([np.min(el.values()) for el in movie_length.values()])-1
    if load:
        f=open(os.path.join(outputFolder, filename), 'r')
        tabF = pickle.load(f)
        f.close()  

    else:
        for xb in dicT:
            for dose in dicT[xb]:
                if dose not in tabF:
                    tabF[dose]={'raw':[], 'speedavg':None, 'proliferation':None, 'speedarr':None, 'trajPerFrame':None}
                print dicT[xb][dose].keys()
                for plate in dicT[xb][dose]:        
                    for well in dicT[xb][dose][plate]:
    
                        dicC = dicT[xb][dose][plate][well]
                        allDisp, trajPerFrame = driftCorrection(movie_length[plate][well]-1, outputFolder, plate, well, dicC.lstTraj, just_average_speed=True)
                        X = allDisp[:,:,0]; Y=allDisp[:,:,1]
                        speed = np.sqrt(X**2+Y**2)
                                
                        tabF[dose]['raw'].append({(plate,well):(speed, trajPerFrame)})
                                            
                        averageCour = np.sum(speed,1)/trajPerFrame[:,0]
                        averageCour = np.array([np.mean(averageCour[frame:frame+6]) for frame in range(MIN_LENGTH_MOVIES-6)], dtype=float)
                        
                        proliferationCour = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(MIN_LENGTH_MOVIES-6)])
                        proliferationCour/=float(np.mean(trajPerFrame[:6,0]))

                        tabF[dose]['speedarr']=speed[:MIN_LENGTH_MOVIES] if tabF[dose]['speedarr'] is None else np.hstack((tabF[dose]['speedarr'], speed[:MIN_LENGTH_MOVIES]))    
                        tabF[dose]['trajPerFrame']=trajPerFrame[:MIN_LENGTH_MOVIES] if tabF[dose]['trajPerFrame'] is None else tabF[dose]['trajPerFrame']+ trajPerFrame[:MIN_LENGTH_MOVIES]           
                        
                        tabF[dose]['speedavg'] = averageCour[:MIN_LENGTH_MOVIES] if tabF[dose]['speedavg'] is None else np.vstack((tabF[dose]['speedavg'], averageCour[:MIN_LENGTH_MOVIES]))
        
    #                    proliferationCour = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(trajPerFrame.shape[0]-6)])                
    #                    proliferationCour/=float(np.mean(trajPerFrame[:6,0]))
                        tabF[dose]['proliferation'] = proliferationCour[:MIN_LENGTH_MOVIES] if tabF[dose]['proliferation'] is None else np.vstack((tabF[dose]['proliferation'], proliferationCour[:MIN_LENGTH_MOVIES]))
                        #plotSpeedHistogramEvolution(plate, well, speed,name,outputFolder)
    if (1 in tabF and  len(tabF[1]['raw'])==0) or (15 in tabF and  len(tabF[15]['raw'])==0):
    #meaning we are dealing with a control
        func=(lambda x:type(x)==int)
    else:
        func=(lambda x:x>0 and type(x)==int)
    for dose in filter(func, tabF):
        try:
            if plot:
                fig=p.figure(figsize=(24,13))
                ax1=fig.add_subplot(221)
                ax2 = fig.add_subplot(222)
                for el in tabF[0]['raw']:
                    pl,w=el.keys()[0]
                    if pl not in seen_plates:
                        seen_plates.append(pl)
                    speed, trajPerFrame=el[(pl,w)]
                    avCour = np.sum(speed,1)/trajPerFrame[:,0]
                    avCour = np.array([np.mean(avCour[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
                    
                    ax1.plot(range(avCour.shape[0]), avCour, color=couleurs[0], label='D {}_{} {}'.format(0, pl, w), ls=linestyles[seen_plates.index(pl)])
                    m=np.argmax(avCour)
                    ax1.text(m, avCour[m], '{} {} {}'.format(m, pl,w))
                    proliferationCour = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(trajPerFrame.shape[0]-6)])
                    proliferationCour/=float(np.mean(trajPerFrame[:6,0]))
                    
                    ax2.plot(range(proliferationCour.shape[0]), proliferationCour,color=couleurs[0], label='D {}_{} {}'.format(0, pl, w),ls=linestyles[seen_plates.index(pl)])
                    
                for el in tabF[dose]['raw']:
                    pl,w=el.keys()[0]
                    speed, trajPerFrame=el[(pl,w)]
                    avCour = np.sum(speed,1)/trajPerFrame[:,0]
                    avCour = np.array([np.mean(avCour[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
                    
                    ax1.plot(range(avCour.shape[0]), avCour, color=couleurs[dose], label='D {}_{} {}'.format(dose, pl, w),ls=linestyles[seen_plates.index(pl)])
                    
                    proliferationCour = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(trajPerFrame.shape[0]-6)])
                    proliferationCour/=float(np.mean(trajPerFrame[:6,0]))
                    
                    ax2.plot(range(proliferationCour.shape[0]), proliferationCour,color=couleurs[dose], label='D {}_{} {}'.format(dose, pl, w),ls=linestyles[seen_plates.index(pl)])
    
                
                ax21=fig.add_subplot(223, sharex=ax1, sharey=ax1)
                ax22 = fig.add_subplot(224, sharex=ax2, sharey=ax2)
                for d in [0, dose]:
                    average=tabF[d]['speedavg'].T; 
                    proliferation=tabF[d]['proliferation'].T
                    average_pooling_cells = np.sum(tabF[d]['speedarr'], 1)/tabF[d]['trajPerFrame'][:MIN_LENGTH_MOVIES,0]
                    #stds=np.sqrt((np.sum(speedT**2,1)/trajPerFrameT[:,0] - average**2))
                    
                    average_pooling_cells = np.array([np.mean(average_pooling_cells[frame:frame+6]) for frame in range(MIN_LENGTH_MOVIES-6)], dtype=float)
                    #stds = np.array([np.mean(stds[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
                    
                    ax21.plot(range(average.shape[0]), average_pooling_cells, label=xenobiotic, color=couleurs[d])
                    
                    try:
                        ax22.plot(range(proliferation.shape[0]), np.mean(proliferation,1), label=xenobiotic, color=couleurs[d])
                        ax21.errorbar(range(average.shape[0]), average_pooling_cells, np.std(average,1))
                        ax22.errorbar(range(proliferation.shape[0]), np.mean(proliferation,1), np.std(proliferation,1))
                    except ValueError:
                        ax22.plot(range(proliferation.shape[0]), proliferation, label=xenobiotic, color=couleurs[d])
                        print "One well hence no error bars"
                    
                ax1.grid(True);ax21.grid(True)
                ax2.grid(True);ax22.grid(True)
                ax1.set_ylim(1.5,6.5); 
                ax2.set_ylim(0.7, 1.7)
                ax1.set_xlim(-5, 220)
                ax2.set_xlim(-5, 220)
                ax2.legend(loc=2); ax1.set_title('Speed average over time {} for dose {}'.format(xenobiotic, dose))
                ax2.set_title("Proliferation over time")    
                p.savefig(os.path.join(outputFolder, '{}_speed_D{}.png'.format(xenobiotic, dose)))
        except AttributeError:
            print dose
            pass        
    #            ax1.plot(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(average[frames[1]:],1))
#            ax1.errorbar(range(frames[1], frames[1]+average[frames[1]:].shape[0]),np.mean(average[frames[1]:],1), np.std(average[frames[1]:],1))
#            
#            ax2.plot(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(proliferation[frames[1]:],1))
#            ax2.errorbar(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(proliferation[frames[1]:],1), np.std(proliferation[frames[1]:],1))
            
    if not load:
        filename='speedOverTime{}.pkl'.format(xenobiotic)
        if filename in os.listdir(outputFolder):
            f=open(os.path.join(outputFolder, filename), 'r')
            tab=pickle.load(f)
            tab.update(tabF)
        else:
            tab=tabF
            
        f=open(os.path.join(outputFolder, filename), 'w')
        pickle.dump(tab, f)
        f.close()  
    
    return tabF

def OLD_featureEvolOverTime(dicT, outputFolder, verbose, movie_length, name, plot=False,
                        load=False, averagingOverWells=False, plates=None):
    tabF={}
    filename='speedOverTime{}.pkl'.format(name)

    if load:
        f=open(os.path.join(outputFolder, filename), 'r')
        tabF = pickle.load(f)
        f.close()  
        
        platesL = tabF.keys()
    else:
        platesL = dicT.keys()
        
    if plates is not None:
        platesL=plates
        
    if plot:
        fig=p.figure(figsize=(24,13))
        ax1=fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        
    for plate in platesL:
        averages=[]; proliferations=[]
        average = None
        proliferation = None
        speedT=None; trajPerFrameT=None
            
        if plate not in tabF:
            tabF[plate]={}
            wells = sorted(dicT[plate].keys())
        else:
            wells = sorted(tabF[plate].keys())    
            
        if plate=='050914':
            w=wells[0]; wells=wells[1:]; wells.append(w)        
            
        for k,well in enumerate(wells):
            k=k%9
            if well not in tabF[plate]:
                tabF[plate][well]=None
            print plate, well
            
            if not load:
                dicC = dicT[plate][well]
                allDisp, trajPerFrame = driftCorrection(movie_length[plate][well]-1, outputFolder, plate, well, dicC.lstTraj, just_average_speed=True)
                X = allDisp[:,:,0]; Y=allDisp[:,:,1]
                speed = np.sqrt(X**2+Y**2)
                
                tabF[plate][well]=(speed, trajPerFrame)
            else:
                speed, trajPerFrame = tabF[plate][well]
            if plot and not averagingOverWells:
                average = np.sum(speed,1)/trajPerFrame[:,0]
                stds=np.sqrt((np.sum(speed**2,1)/trajPerFrame[:,0] - average**2))
                
                average = np.array([np.mean(average[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
                stds = np.array([np.mean(stds[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
                
                ax1.plot(range(average.shape[0]), average, color=couleurs[k], label='D {}'.format(wells.index( well)%9))
                #ax1.errorbar(range(average.shape[0]), average, stds,color=couleurs[k])
                
                proliferation = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(trajPerFrame.shape[0]-6)])
                proliferation/=float(np.mean(trajPerFrame[:6,0]))
                
                ax2.plot(range(proliferation.shape[0]), proliferation,color=couleurs[k], label='D {}'.format(wells.index( well)%9))
                
                ax1.text(average.shape[0]+1, average[-1], 'D {}'.format(wells.index( well)%9))
                ax2.text(average.shape[0]+1, proliferation[-1], 'D {}'.format(wells.index( well)%9))
                
                averages.append(average); proliferations.append(proliferation)
                
            if averagingOverWells:
                averageCour = np.sum(speed,1)/trajPerFrame[:,0]
                averageCour = np.array([np.mean(averageCour[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
                speedT=speed if speedT is None else np.hstack((speedT, speed))    
                trajPerFrameT=trajPerFrame if trajPerFrameT is None else trajPerFrameT+ trajPerFrame           
                
                average = averageCour if average is None else np.vstack((average, averageCour))

                proliferationCour = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(trajPerFrame.shape[0]-6)])                

                proliferationCour/=float(np.mean(trajPerFrame[:6,0]))
                proliferation = proliferationCour if proliferation is None else np.vstack((proliferation, proliferationCour))
            #plotSpeedHistogramEvolution(plate, well, speed,name,outputFolder)
        
        if plot and averagingOverWells:
            average=average.T; 
            proliferation=proliferation.T
            average_pooling_cells = np.sum(speedT, 1)/trajPerFrameT[:,0]
            #stds=np.sqrt((np.sum(speedT**2,1)/trajPerFrameT[:,0] - average**2))
            
            average_pooling_cells = np.array([np.mean(average_pooling_cells[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
            #stds = np.array([np.mean(stds[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
            
            ax1.plot(range(average.shape[0]), average_pooling_cells, label=name)
            ax1.errorbar(range(average.shape[0]), average_pooling_cells, np.std(average,1))
            
            ax2.plot(range(proliferation.shape[0]), np.mean(proliferation,1), label=name)
            ax2.errorbar(range(proliferation.shape[0]), np.mean(proliferation,1), np.std(proliferation,1))
            
#            ax1.plot(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(average[frames[1]:],1))
#            ax1.errorbar(range(frames[1], frames[1]+average[frames[1]:].shape[0]),np.mean(average[frames[1]:],1), np.std(average[frames[1]:],1))
#            
#            ax2.plot(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(proliferation[frames[1]:],1))
#            ax2.errorbar(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(proliferation[frames[1]:],1), np.std(proliferation[frames[1]:],1))
            
        if not load:
            filename='speedOverTime{}.pkl'.format(name)
            if filename in os.listdir(outputFolder):
                f=open(os.path.join(outputFolder, filename), 'r')
                tab=pickle.load(f)
                tab.update(tabF)
            else:
                tab=tabF
                
            f=open(os.path.join(outputFolder, filename), 'w')
            pickle.dump(tab, f)
            f.close()  
    if plot:
        ax1.grid(True)
        ax2.grid(True)
        ax1.set_ylim(1.5,6.5)
        ax2.set_ylim(0.7, 1.7)
        ax1.set_xlim(-5, 200)
        ax2.set_xlim(-5, 200)
        ax1.legend(loc=2); ax1.set_title('Speed average over time {}'.format(name))
        #ax1.set_xticklabels(['24', '36.5', '49', '61.5', '74'])
        #ax2.set_xticklabels(['24', '36.5', '49', '61.5', '74'])
        ax2.set_title("Proliferation over time")
    if averagingOverWells:
        p.savefig(os.path.join(outputFolder, '{}_speed_{}.png'.format(plate,name+'_average')))
    else:
        p.savefig(os.path.join(outputFolder, '{}_speed_{}.png'.format(plate, name)))
    return np.array(averages), np.array(proliferations)


def featureEvolOverTime_august(dicT, connexions, outputFolder, verbose, movie_length, name, plot=False, load=False, averagingOverWells=False):
    tabF={}
    filename='speedOverTime{}.pkl'.format(name)

    if load:
        f=open(os.path.join(outputFolder, filename), 'r')
        tabF = pickle.load(f)
        f.close()  
        
        platesL = tabF.keys()
    else:
        platesL = dicT.keys()
        
    for plate in platesL:
        averages=[]; proliferations=[]
        average = None
        proliferation = None
        
        if plate =='150814':
            frames=[41,129,185]
            
        elif plate=='130814':
            frames=[2,75,149]
            
        if plate not in tabF:
            tabF[plate]={}
            wells = sorted(dicT[plate].keys())
        else:
            wells = sorted(tabF[plate].keys())

        if plot:
            import matplotlib.pyplot as p
            fig=p.figure(figsize=(24,13))
            ax1=fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            
            
        for k,well in enumerate(wells):
            if well not in tabF[plate]:
                tabF[plate][well]=None
            print plate, well
            
            if not load:
                dicC = dicT[plate][well]
                allDisp, trajPerFrame = driftCorrection(movie_length[plate][well]-1, outputFolder, plate, well, dicC.lstTraj, just_average_speed=True)
                X = allDisp[:,:,0]; Y=allDisp[:,:,1]
                speed = np.sqrt(X**2+Y**2)
                
                tabF[plate][well]=(speed, trajPerFrame)
            else:
                speed, trajPerFrame = tabF[plate][well]
                
            if plot and not averagingOverWells:
                average = np.sum(speed,1)/trajPerFrame[:,0]
                average = np.array([np.mean(average[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
                
                ax1.plot(range(frames[0]), average[:frames[0]], color=couleurs[k], label='D {}'.format(wells.index( well)))
                ax1.plot(range(frames[1],frames[1]+average[frames[1]:].shape[0]), average[frames[1]:], color=couleurs[k])
                
                proliferation = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(trajPerFrame.shape[0]-6)])
                proliferation/=float(np.mean(trajPerFrame[:6,0]))
                
                ax2.plot(range(frames[0]), proliferation[:frames[0]],color=couleurs[k], label='D {}'.format(wells.index( well)))
                ax2.plot(range(frames[1], frames[1]+average[frames[1]:].shape[0]), proliferation[frames[1]:],color=couleurs[k])
                
                ax1.text(frames[2]+2, average[-1], 'D {}'.format(wells.index( well)))
                ax2.text(frames[2]+2, proliferation[-1], 'D {}'.format(wells.index( well)))
                
                averages.append(average); proliferations.append(proliferation)
                
            if averagingOverWells:
                averageCour = np.sum(speed,1)/trajPerFrame[:,0]
                if plate=='130814':
                    averageCour = np.array([np.mean(averageCour[frame:frame+6]) for frame in range(149)], dtype=float)
                else:
                    averageCour = np.array([np.mean(averageCour[frame:frame+6]) for frame in range(trajPerFrame.shape[0]-6)], dtype=float)
                average = averageCour if average is None else np.vstack((average, averageCour))

                if plate=='130814':
                    proliferationCour = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(149)])
                else:
                    proliferationCour = np.array([np.mean(trajPerFrame[frame:frame+6,0]) for frame in range(trajPerFrame.shape[0]-6)])                

                proliferationCour/=float(np.mean(trajPerFrame[:6,0]))
                proliferation = proliferationCour if proliferation is None else np.vstack((proliferation, proliferationCour))
            #plotSpeedHistogramEvolution(plate, well, speed,name,outputFolder)
        
        if plot and averagingOverWells:
            average=average.T; proliferation=proliferation.T
            ax1.plot(range(frames[0]), np.mean(average[:frames[0]],1), label=name)
            ax1.errorbar(range(frames[0]), np.mean(average[:frames[0]],1), np.std(average[:frames[0]],1))
            
            ax2.plot(range(frames[0]), np.mean(proliferation[:frames[0]],1), label=name)
            ax2.errorbar(range(frames[0]), np.mean(proliferation[:frames[0]],1), np.std(proliferation[:frames[0]],1))
            
            ax1.plot(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(average[frames[1]:],1))
            ax1.errorbar(range(frames[1], frames[1]+average[frames[1]:].shape[0]),np.mean(average[frames[1]:],1), np.std(average[frames[1]:],1))
            
            ax2.plot(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(proliferation[frames[1]:],1))
            ax2.errorbar(range(frames[1], frames[1]+average[frames[1]:].shape[0]), np.mean(proliferation[frames[1]:],1), np.std(proliferation[frames[1]:],1))
            
            
        if not load:
            filename='speedOverTime{}.pkl'.format(name)
            f=open(os.path.join(outputFolder, filename), 'w')
            pickle.dump(tabF, f)
            f.close()  
        if plot:
            ax1.grid(True)
            ax2.grid(True)
            ax1.set_ylim(1.5,6.5)
            ax2.set_ylim(0.7, 1.7)
            ax1.set_xlim(-5, 200)
            ax2.set_xlim(-5, 200)
            ax1.legend(loc=2); ax1.set_title('Speed average over time {}'.format(name))
            #ax1.set_xticklabels(['24', '36.5', '49', '61.5', '74'])
            #ax2.set_xticklabels(['24', '36.5', '49', '61.5', '74'])
            ax2.set_title("Proliferation over time")
            if averagingOverWells:
                p.savefig(os.path.join(outputFolder, 'speed_{}.png'.format(name+'_average')))
            else:
                p.savefig(os.path.join(outputFolder, 'speed_{}.png'.format(name)))
        return np.array(averages), np.array(proliferations)

def compaSegmentation(hdf5Folder, plate, sh=False, number=False): 
    '''
    Outputs the average number of well-segmented objects per frame and the percentage of well-segmented objects per frame
    
    '''
    result={}
    
    for well in filter(lambda x: 'hdf5' in x, os.listdir(hdf5Folder)):
        print well
        path = "/sample/0/plate/"+plate+"/experiment/"+well[:5]+"/position/1/feature/primary__primary/object_classification/prediction"
        classification = np.array(vi.readHDF5(os.path.join(hdf5Folder, well), path), dtype=int)
        #num_objects = len(classification)
        
        pathObjects= "/sample/0/plate/"+plate+"/experiment/"+well[:5]+"/position/1/object/primary__primary"
        objects = vi.readHDF5(os.path.join(hdf5Folder, well), pathObjects)
        frameList = np.array(objects, dtype=int)
        num_frames = np.max(frameList)    
        classif_globale = np.sum(classification)/float(num_frames)    
#        for k in range(1, num_frames+1):
#            if len(np.where(frameList==k)[0])==0:
#                print k
#                pdb.set_trace()
        classif_perFrame = [np.sum(classification[np.where(frameList==k)])/float(len(np.where(frameList==k)[0])) for k in range(1, num_frames+1)] if not number else [np.sum(classification[np.where(frameList==k)]) for k in range(1, num_frames+1)]

        result[well]=(classif_globale, classif_perFrame)
        
    if sh:
        from util.plots import couleurs
        import matplotlib.pyplot as p
        f=p.figure(); ax=f.add_subplot(111)
        for k,w in enumerate(sorted(result.keys())):
            ax.scatter(range(174), result[w][1], color=couleurs[k], label=int(w[:5]))
            ax.text(180, result[w][1][-1], w[:5])
            ax.legend(); ax.grid(True); ax.set_xlim(0,200)
            if not number:
                ax.set_ylim(0.3, 1)
            else:
                ax.set_ylim(80, 350)
        p.show()
    return result

def allDataExtraction(dataFolder):
    newFrameLot = None
    
    listD = filter(lambda x: os.path.isdir(os.path.join(dataFolder, x)), os.listdir(dataFolder))

    for plate in listD:
        listW = filter(lambda x: '.hdf5' in x, os.listdir(os.path.join(dataFolder, plate, 'hdf5')))
        for filename in listW:
            well=filename[:-5]
                
            filename = os.path.join(dataFolder, plate,"hdf5", filename)
            
#ajout du frameLot et du tabF
            frameLotC, _ = test.gettingRaw(filename, None, plate, well, secondary=True)
            if frameLotC==None:
                sys.stdout.write("File {} containing data for plate {}, well {} does not contain all necessary data".format(filename, plate, well))
                continue
            if newFrameLot == None:
                newFrameLot = frameLotC 
            else: newFrameLot.addFrameLot(frameLotC)
    return newFrameLot 

def featureComparison(table, index):
    primary = table[:,index]; secondary =table[:,FEATURE_NUMBER+index]+0.00001
    return np.divide(primary,secondary)

def targetedDataExtraction(frameLot, features):
    R={}
    for plate in frameLot.lstFrames:
        R[plate]={}
        for well in frameLot.lstFrames[plate]:
            well2=well[:-3]
            result={'cell_count':[], 'red_only_dist':[]}
            result.update({'{}_ch{}'.format(arg[0], arg[1]+1):[] for arg in features})
            
            for frame_nb in frameLot.lstFrames[plate][well]:
                frame = frameLot.lstFrames[plate][well][frame_nb]
                result["cell_count"].append(frame.centers.shape[0])
                result["red_only_dist"].append(featureComparison(frame.features, cecog_features.index("roisize")))
                for arg in features:
                    if arg[0] not in cecog_features:
                        print "Argument non understood"
                        continue
                    if arg[0]=='irregularity':
                        result['{}_ch{}'.format(arg[0], arg[1]+1)].append(frame.features[:,arg[1]*FEATURE_NUMBER+cecog_features.index(arg[0])]+1)
                        continue
                    result['{}_ch{}'.format(arg[0], arg[1]+1)].append(frame.features[:,arg[1]*FEATURE_NUMBER+cecog_features.index(arg[0])])
            R[plate][well2]=result
                    
    return R

def quickPlot(dict_, bins=20):
    for plate in dict_:
        for well in sorted(dict_[plate]):
            result = dict_[plate][well]
        
            f=p.figure()
            movie_length =len(result['cell_count'])
            ax=f.add_subplot(2,2,0)
            ax.scatter(range(movie_length), result['cell_count'])
            ax.set_title('cell_count')
            
            i=1
            for el in result:
                if el is not 'cell_count':
                    ax2=f.add_subplot(2,2,i)
                    ax2.set_title('{}'.format(el))
                    n,bins, patches=ax2.hist([result[el][int(movie_length*point)] for point in [0,0.25, 0.50, 0.75]], range=(0,2),
                                             bins=bins, normed=False, 
                                             label=['t={}'.format(point) for point in [0,0.25, 0.50, 0.75] ])
                    i+=1
            p.legend()
            p.suptitle('Plate {}, well {}'.format(plate, well))
            p.show()    
                
def comparPlot(dict_, feature, bins=20, folderPl="/media/lalil0u/New/Documents/Toxicology/Plans_plaques"):
    timepoints = [0,0.25, 0.50, 0.75]
    
    for plate in dict_:
        l=sorted(dict_[plate])
        wellLabels = readPlateSetting([plate], folderPl)
        counts = []
        feat_list=[]
        labels=[]
        movie_length=0
        for well in l:
            result = dict_[plate][well]
            if feature not in result:
                continue
            if folderPl is not None:
                labels.append(wellLabels[int(well)])
            else:
                labels.append( well)
            counts.append(result['cell_count'])
            movie_length = max(movie_length, len(result['cell_count']))
            feat_list.append(result[feature])
            
    f=p.figure()
    ax=f.add_subplot(3,2,0)
    ax.set_title('cell_count')
    axes=[]; i=1 
    for point in timepoints:
        zz=f.add_subplot(3,2,i)
        zz.set_title('t={}'.format(point))
        axes.append(zz)
        i+=1
        
    for i, el in enumerate(counts):
        if len(el)<movie_length:
            el=el+list(repeat(0, movie_length-len(el)))
        
        ax.scatter(range(movie_length), el, color ='#'+ couleurs[i], label =labels[i])
    
    for k,point in enumerate(timepoints):
        axes[k].hist([feat_list[i][int(movie_length*point)] for i in range(len(counts))], color=['#'+couleurs[i] for i in range(len(counts))], 
                     label=labels, bins=bins, normed=True)
                
    p.legend()
    p.suptitle('Comparison between wells of {}'.format(feature))
    p.show()
    
def plotPlateAsHeatmap(data):
    #data.shape = (nb_row, nb_col)
    red_purple = brewer2mpl.get_map('RdPu', 'Sequential', 9).mpl_colormap
    
    fig, ax = p.subplots(1)
    
    ax.pcolormesh(data.T,
                  cmap=red_purple)
    xticks = np.array(range(data.shape[0]))
    yticks = np.array(range(data.shape[1]))
    
    yticklabels=list(string.uppercase[:data.shape[1]]); yticklabels.reverse()
    xticklabels=range(1,data.shape[0]+1)
    
    ax.set_xticks(xticks+0.5)
    ax.set_xticklabels(xticklabels)
    
    ax.set_yticks(yticks+0.5)
    ax.set_yticklabels(yticklabels)
    ax.set_title("zou")
    ax.legend()
#    fig.savefig('pcolormesh_prettyplotlib_labels_other_cmap_sequential.png')
