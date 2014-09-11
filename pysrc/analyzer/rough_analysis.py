import os, sys, pdb, string
import cPickle as pickle
#import prettyplotlib as ppl
import brewer2mpl
import matplotlib.pyplot as p
import numpy as np
import vigra.impex as vi

from itertools import repeat

from analyzer import cecog_features
#from interface import fromPlateToPython
from tracking import test
from tracking.importPack import FEATURE_NUMBER
from tracking.plots import markers
from tracking.trajPack.clustering import couleurs
from tracking.trajPack.trajFeatures import driftCorrection

from plateSetting import readPlateSetting


def plotSpeedHistogramEvolution(plate, well, speed, name, outputFolder):
    if not os.path.isdir(os.path.join(outputFolder, name)):
        os.mkdir(os.path.join(outputFolder,name))
    import matplotlib.pyplot as p
    for frame in range(speed.shape[0]-6):
        f=p.figure()
        ax=f.add_subplot(111)
        arr = speed[frame:frame+6].flatten()
        n,bins,patches = ax.hist(arr, bins=50,normed=False, range=(0.5,30))
        ax.set_xlim(0,30)
        ax.set_ylim(0,75)
        f.savefig(os.path.join(outputFolder, name, 'hist_speed_{}_{:>05}.png'.format(well,frame)))

    return

def featureEvolOverTime(dicT, connexions, outputFolder, verbose, movie_length, name, plot=False, load=False, averagingOverWells=False):
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
        return average, proliferation

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
