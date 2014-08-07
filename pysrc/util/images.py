import os
import vigra
import vigra.impex as vi
import numpy as np

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
            
def quantileObjects(plate, rawFolder="/media/lalil0u/New/data/Xb_screen/Screen_trials/", wellL=None):
    
    if wellL == None:
        wellL = sorted(os.listdir(os.path.join(rawFolder, plate)))
    
