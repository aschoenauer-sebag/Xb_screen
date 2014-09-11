import os, pdb
import vigra
import vigra.impex as vi
import numpy as np
import cPickle as pickle
from collections import Counter
import matplotlib.pyplot as p
from util.plots import couleurs
from tracking.importPack import imp 

def checkNormalization(folder, plate, lengthMovie = 192):
    '''
    This function checks the max intensity on both channels on all images for all wells of a given plate.
    Maxima should be around 1,000 or 3,000 depending on the channel.
    
    If images have a maximum under 200 it is black (from experience).
    '''
    arr=np.zeros(shape=(len(os.listdir(folder)), lengthMovie, 2))
    
    if 'CR_normalization.pkl' not in os.listdir(folder):
        for w, well in enumerate(os.listdir(folder)):
            imageName = '%s--%s--P0001_t00{:>03}_c0000{}.tif'%(plate, well)
            for timepoint in range(1,lengthMovie+1):
                for ch in [1,2]:
                    im = vi.readImage(os.path.join(well, imageName.format(timepoint, ch)))
                    arr[w, timepoint-1, ch-1]=np.max(im)
        f=open(os.path.join(folder, 'CR_normalization.pkl'), 'w')
        pickle.dump(arr, f); f.close()
    else:
        f=open(os.path.join(folder, 'CR_normalization.pkl'))
        arr=pickle.load(f); f.close()
        
    print 'Compte rendu, plate ', plate
    total_missing = []
    for k in range(arr.shape[0]):
        print '------Well ', k
        l=sorted(Counter(np.where(arr[k]<200)[0]).keys())
        total_missing.extend(l)
        print 'Frames where no images? \n', l
        ll=filter(lambda x: x not in l, range(lengthMovie))
        print 'Else maximun average', np.mean(arr[k, ll,:],0), 'std ', np.std(arr[k, ll,:],0)
    
    summary = np.bincount(total_missing, minlength=lengthMovie)
    p.scatter(range(lengthMovie), summary )
    p.grid(True)
    p.show()
    print summary
    return summary

def imageNormalization(wellL, rawFolder):
    '''
    This function makes float images integer images
    '''
    for puits in wellL:
        for image in sorted(os.listdir(os.path.join(rawFolder, puits))):
            im = vi.readImage(os.path.join(rawFolder, puits, image))
            normImage = vigra.VigraArray(im.shape, dtype=np.dtype('uint16'))
            normImage[:,:,0]=np.array(im[:,:,0], dtype=int)
            vi.writeImage(normImage,os.path.join(rawFolder, puits, image))
            
def quantileObjects(plate, hdfFolder="/media/lalil0u/New/projects/Xb_screen/plates/", wellL=None, sh=False):
    '''
    Here the goal is to get the intensity distributions in objects to better normalize the images
    
    '''
    if wellL is None:
        wellL = sorted(os.listdir(os.path.join(hdfFolder, plate, 'hdf5')))
        wellL=[el[:-8] for el in wellL]
    print wellL
        
    r={}
    for well in wellL:
        print well
        feature = 'n2_avg'
        filename = os.path.join(hdfFolder, plate, 'hdf5', '%s_01.hdf5'%well)
        result = imp.importTargetedFromHDF5(filename, plate, '%s_01'%well,[feature], secondary=True, secondary_outside=False)
        result= result[1].lstFrames[plate]['%s_01'%well]#une liste de frames, les informations sont dans result[frame].features
        rCy3 = []; rGFP=[]
        for frame in result:
            rCy3.extend(result[frame].features[:,0])
            try:
                rGFP.extend(result[frame].features[:,1])
            except IndexError:
                pass
            
        r[well]=[rCy3, rGFP]
    import matplotlib.pyplot as p

    f=p.figure()
    ax=f.add_subplot(211)
    for k,w in enumerate(wellL):
        n,bins,patches = ax.hist(r[w][0], bins=100, range=(0,255), color=couleurs[k%len(couleurs)], alpha=0.5, label = w)
    ax.set_xlim((0,255))
    ax.set_title('Cy3 pixel intensity distributions')
    ax.legend()
    ax=f.add_subplot(212)
    for k,w in enumerate(wellL):
        try:
            n,bins,patches = ax.hist(r[w][1], bins=100, range=(0,255), color=couleurs[k%len(couleurs)], alpha=0.5, label = w)
        except IndexError:
            pass
    ax.set_xlim((0,255))
    ax.set_title('GFP pixel intensity distributions')
    ax.legend()
    if sh:
        p.show()
    else:
        p.savefig(os.path.join(hdfFolder, 'intensity_dist_{}_{}.png'.format(plate, well)))
        
    return r
    
