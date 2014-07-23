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
            imgInDir = filter(lambda x: x[:3] == exp[1][2:5], os.listdir(os.path.join(inDir, exp[0])))[0]
            imgInDir = os.path.join(inDir, exp[0], imgInDir)
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
    
#    for imageName in lstImageNames:
#        print imageName
#        img = vi.readImage(os.path.join(imgDir, imageName))
#        pdb.set_trace()
#        normImage = vigra.VigraArray(img.shape, dtype=np.dtype('uint8'))
##WARNING if you only do normImage = (img - etc then we have a flickering effect. Apparently vigra decides to do its normalization on every image as it pleases
#        suffix = imageName.split('.')[-1]
#        vi.writeImage(normImage, os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg')), dtype = np.dtype('uint8'))
    
    # movie filename
    movieName = '{}_P{}_W{}.avi'.format(gene, plate[:9], well)
    
    # make output directory
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    
    # encode command
    #encode_command = 'mencoder "mf://%s/*.jpg" -mf fps=3 -o %s -ovc xvid -oac copy -xvidencopts fixed_quant=2.5'
    encode_command = "mencoder mf://%s/*.png -mf w=800:h=600:fps=25:type=png -ovc copy -oac copy -o %s"
    encode_command %= (imgDir, os.path.join(outDir, movieName))
    print encode_command
    os.system(encode_command)
    
    # cleaning up temporary directory
    shell_command = 'rm %s/*.jpg' % tempDir
    print shell_command
    os.system(shell_command)

    return

def makeMovieMultiChannels(imgDir, outDir,plate, well, channels=[2,1], tempDir=None, redo=True, doChannel1Only=True):
    '''
    From a list of images containing information for two channels in different images, make 
    a movie. It is assumed that the name of a file is something like
     '[conditions]_--W[well]--P0001_t[frame number]_c[channel number].tif'
    
    It is not possible I think that there are more than three channels ? However
    it is possible to have 2 or three channels
    '''

    # movie filename
    movieName = 'P{}_W{}'.format(plate, well)
    
    if not redo:
        if os.path.isdir(outDir) and movieName+'_1.avi' in os.listdir(outDir):
            return
    
    # make output directory
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    # temp directory
    if tempDir is None:
        tempDir = os.path.join(outDir, 'temp')
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
    if len(channels)>1:
        lstImage={ch:filter(lambda x: 'tif' in x and int(x.split('_')[-1][1:6])==ch, os.listdir(imgDir)) for ch in channels}#os.listdir(imgDir)
    else:
        lstImage={channels[0]:os.listdir(imgDir)}
    
#first I look at the min pixel values on all channels for subsequent background substraction
    min_={ch : 500 for ch in channels}
    max_ = {ch : 0 for ch in channels}
    for imageName in lstImage[channels[0]]:
        img = vi.readImage(os.path.join(imgDir, imageName))
        min_[channels[0]]=min(np.min(img), min_[channels[0]])
        max_[channels[0]]=max(np.max(img), max_[channels[0]])
        
        suffix = imageName.split('_')[-1]
        for ch in channels[1:]:
            imageName2 = os.path.basename(imageName).replace(suffix, 'c{:>05}.tif'.format(ch))
            im = vi.readImage(os.path.join(imgDir, imageName2))
            min_[ch]=min(np.min(im), min_[ch])
            max_[ch]=max(np.max(im), max_[ch])
            
    if doChannel1Only:
        #making movie for channel 1 only (nucleus)
        for imageName in lstImage[channels[0]]:
            img = vi.readImage(os.path.join(imgDir, imageName))
            normImage = vigra.VigraArray(img.shape, dtype=np.dtype('uint8'))
#WARNING if you only do normImage = (img - etc then we have a flickering effect. Apparently vigra decides to do its normalization on every image as it pleases
            normImage[:,:,0] = (img[:,:,0]-min_[channels[0]])*(2**8-1)/(max_[channels[0]]-min_[channels[0]])
            suffix = imageName.split('.')[-1]
            vi.writeImage(normImage, os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg')), dtype = np.dtype('uint8'))
    
        # encode command
        encode_command = 'mencoder "mf://%s/*.jpg" -mf fps=3 -o %s -ovc xvid -oac copy -xvidencopts fixed_quant=2.5'
        encode_command %= (tempDir, os.path.join(outDir, '{}_1.avi'.format(movieName)))
        print encode_command
        os.system(encode_command)
        
        # cleaning up temporary directory
        shell_command = 'rm %s/*.jpg' % tempDir
        print shell_command
        os.system(shell_command)
        
    #making movie for both channels    
    for imageName in lstImage[channels[0]]:        
        img = vi.readImage(os.path.join(imgDir, imageName))
        shape = img.shape
        colorImage = vigra.VigraArray((shape[0],shape[1], 3), dtype=np.dtype('uint8'))

        colorImage[:,:,0] = (img[:,:,0]-min_[channels[0]])*(2**8-1)/(max_[channels[0]]-min_[channels[0]])
        
        suffix = imageName.split('_')[-1]
        for i,ch in enumerate(channels[1:]):
            imageName2 = os.path.basename(imageName).replace(suffix, 'c{:>05}.tif'.format(ch))
            im = vi.readImage(os.path.join(imgDir, imageName2))

            colorImage[:,:,i+1] = (im[:,:,0]-min_[ch])*(2**8-1)/(max_[ch]-min_[ch])
            
        suffix = imageName.split('.')[-1]
        vi.writeImage(colorImage, os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg')), dtype = np.dtype('uint8'))

    
    # encode command
    encode_command = 'mencoder "mf://%s/*.jpg" -mf fps=3 -o %s -ovc xvid -oac copy -xvidencopts fixed_quant=2.5'
    encode_command %= (tempDir, os.path.join(outDir, '{}_2.avi'.format(movieName)))
    print encode_command
    os.system(encode_command)
    
    # cleaning up temporary directory
    shell_command = 'rm %s/*.jpg' % tempDir
    print shell_command
    os.system(shell_command)

    return