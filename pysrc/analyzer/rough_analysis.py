import os, sys, pdb, string

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

from plateSetting import readPlateSetting

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
