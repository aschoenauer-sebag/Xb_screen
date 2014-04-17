import vigra.impex as vi
import vigra
import os, pdb
import numpy as np
from PIL import Image

def makeMovieFromExpDict(idDict, tempDir=None, inDir='/media/lalil0u/New/workspace2/Tracking/cluster', \
                         outDir='/media/lalil0u/New/workspace2/Tracking/cluster'):
    if tempDir is None:
        tempDir = os.path.join(outDir, 'temp')
    
    for gene in idDict:
        for exp in idDict[gene]:
            imgInDir = os.path.join(inDir, gene, exp[0]+'_'+exp[1][2:5])
            makeMovie(imgInDir, outDir, gene, exp[0], exp[1], tempDir)
    return

def makeMovie(imgDir, outDir,gene, plate, well, tempDir=None):
    def normWrite(img, filename):
        img=(img-2**15)*(2**8-1)/(2**12-1)
        vi.writeImage(img, filename)
        return 1
    # temp directory
    if tempDir is None:
        tempDir = os.path.join(outDir, 'temp')
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
        
    lstImageNames=os.listdir(imgDir)
    
    for imageName in lstImageNames:
        img = Image.open(os.path.join(imgDir, imageName))#vi.readImage(os.path.join(imgDir, imageName), dtype='DOUBLE')
        suffix = imageName.split('.')[-1]
#pour une image en deux couleurs il faut enregistrer un canal sur une couleur, un autre sur une autre et avoir une image RGB que l'on sauve
        if np.max(img)<256:
            img.save(os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg')))#vi.writeImage(img, os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg')))
        else:
            print "normalization"
            normWrite(img, os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg')))
    
    # movie filename
    movieName = '{}_P{}_W{}.avi'.format(gene, plate[:9], well)
    
    # make output directory
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    
    # encode command
    encode_command = 'mencoder "mf://%s/*.jpg" -mf fps=3 -o %s -ovc xvid -oac copy -xvidencopts fixed_quant=2.5'
    encode_command %= (tempDir, os.path.join(outDir, movieName))
    print encode_command
    os.system(encode_command)
    
    # cleaning up temporary directory
    shell_command = 'rm %s/*.jpg' % tempDir
    print shell_command
    os.system(shell_command)

    return

def makeMovieMultiChannels(imgDir, outDir,plate, well, channels=[1,2], tempDir=None):
    '''
    From a list of images containing information for two channels in different images, make 
    a movie. It is assumed that the name of a file is something like
     '[conditions]_--W[well]--P0001_t[frame number]_c[channel number].tif'
    
    It is not possible I think that there are more than three channels ? However
    it is possible to have 2 or three channels
    '''
    def normWrite(img, filename):
        img=(img-2**15)*(2**8-1)/(2**12-1)
        vi.writeImage(img, filename)
        return 1

    # temp directory
    if tempDir is None:
        tempDir = os.path.join(outDir, 'temp')
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
        
    lstImage={ch:filter(lambda x: int(x.split('_')[-1][1:6])==ch, os.listdir(imgDir)) for ch in channels}#os.listdir(imgDir)
    
#first I look at the min pixel values on all channels for subsequent background substraction
#    min_={ch : 0 for ch in channels}
#    for imageName in lstImage[channels[0]]:
#        img = vi.readImage(os.path.join(imgDir, imageName))
#        
#        if np.min(img)<min_[channels[0]]:
#            min_[channels[0]]=np.min(img)
#        
#        suffix = imageName.split('_')[-1]
#        for ch in channels[1:]:
#            imageName2 = os.path.basename(imageName).replace(suffix, 'c{:>05}.tif'.format(ch))
#            im = vi.readImage(os.path.join(imgDir, imageName2))
#            if np.min(img)<min_[channels[0]]:
#                min_[ch]=np.min(im)
                
    
    for imageName in lstImage[channels[0]]:        
        img = vi.readImage(os.path.join(imgDir, imageName))
        shape = img.shape
        
 #this is the image that will be saved : 3 for RGB
        colorImage = vigra.VigraArray((shape[0],shape[1], 3))#, dtype=np.int8)
        colorImage[:,:,0] = img[:,:,0]*(2**8-1)/(2**12-1)
        
        suffix = imageName.split('_')[-1]
        for i,ch in enumerate(channels[1:]):
            imageName2 = os.path.basename(imageName).replace(suffix, 'c{:>05}.tif'.format(ch))
            im = vi.readImage(os.path.join(imgDir, imageName2))
            colorImage[:,:,i+1] = im[:,:,0]*(2**8-1)/(2**12-1)
            
        suffix = imageName.split('.')[-1]
        vi.writeImage(colorImage, os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg')))#, dtype = 'UINT8')

    
    # movie filename
    movieName = 'P{}_W{}.avi'.format(plate, well)
    
    # make output directory
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    
    # encode command
    encode_command = 'mencoder "mf://%s/*.jpg" -mf fps=3 -o %s -ovc xvid -oac copy -xvidencopts fixed_quant=2.5'
    encode_command %= (tempDir, os.path.join(outDir, movieName))
    print encode_command
    os.system(encode_command)
    
    # cleaning up temporary directory
    shell_command = 'rm %s/*.jpg' % tempDir
    print shell_command
    os.system(shell_command)

    return