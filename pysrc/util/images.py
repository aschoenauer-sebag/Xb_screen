import os, pdb
import vigra
import vigra.impex as vi
import numpy as np

from util.plots import couleurs
from tracking.importPack import imp 

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
    
