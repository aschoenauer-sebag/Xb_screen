import vigra.impex as vi
import vigra
import os, pdb, shutil
import numpy as np

def makeClassifMovieFromExpDict(idDict, tempDir = None, inDir = '/share/data40T/aschoenauer/drug_screen',\
                                outDir = "/cbio/donnees/aschoenauer/public_html/interface_screen/plates/static/movies",\
                                clef=lambda x:int(x.split('_')[2][1:-4]),folderName='analyzed/{:>05}_01/images/primary_classification_primary',\
                                extension="jpg", redo=True, rename_pl=lambda x:x):
    
    if tempDir is None:
        tempDir = os.path.join(outDir, 'temp')
        
    if idDict==None:
        idDict={3:[]}
        for pl in os.listdir(inDir):
            idDict[3].extend([(pl, el) for el in range(1,384)])
    
    for gene in idDict:
        for pl,w in idDict[gene]:
            print pl,w
#             if len(w.split('_'))==1 or w.split('_')[1]!='01':
#                 w=w+'_01'
#            imgInDir = os.path.join(inDir, pl, filter(lambda x: w[2:5] ==x[:3], os.listdir(os.path.join(inDir, pl)))[0])
            imgInDir=os.path.join(inDir, pl, folderName.format(w))
            try:
                makeMovieWithoutRenorm(imgInDir, outDir, gene, rename_pl(pl),w, clef, tempDir,extension=extension, redo=redo)
            except OSError:
                print "No folder ", imgInDir
                continue
    return


def makeRawMovieFromExpDict(idDict, tempDir=None, inDir='/share/data20T/mitocheck/compressed_data', \
                         outDir='/cbio/donnees/aschoenauer/cluster', clef=lambda x:int(x.split('--')[3][1:])):
    if tempDir is None:
        tempDir = os.path.join(outDir, 'temp')
    
    for gene in idDict:
        for exp in idDict[gene]:
            print exp[:9], exp[11:]
            folder = filter(lambda x: x[:9] == exp[:9], os.listdir(inDir))[0]
            imgInDir = filter(lambda x: x[:3] == exp[13:], os.listdir(os.path.join(inDir, folder)))[0]
            imgInDir = os.path.join(inDir, folder, imgInDir)
            makeMovieWithoutRenorm(imgInDir, outDir, gene, exp[:9], exp[13:], clef, tempDir)
    return

def makeMovieWithoutRenorm(imgDir, outDir,gene, plate, well, clef, tempDir=None, extension="png", redo=True):

    # movie filename
    movieName = 'P{}_W{}_{}.avi'.format(plate, well, gene)
    # make output directory
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    elif not redo and movieName in os.listdir(outDir):
        print "Movie {} already done".format(movieName)
        return
    
    # temp directory
    if tempDir is None:
        tempDir = os.path.join(outDir, 'temp')
    
    lstImageNames=os.listdir(imgDir)
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
    lstImageNames.sort(key=clef)
    for el in lstImageNames:
        shutil.copy(os.path.join(imgDir, el), tempDir)   
        
    # encode command
    encode_command = "mencoder mf://%s/*.%s -mf w=800:h=600:fps=3:type=%s -ovc copy -oac copy -o %s"
    encode_command %= (tempDir,extension, extension, os.path.join(outDir, movieName))
    print encode_command
    os.system(encode_command)
    
    # cleaning up temporary directory
    shell_command = 'rm %s/*.%s' %(tempDir, extension)
    print shell_command
    os.system(shell_command)

    return

def makeMovie(imgDir, outDir, plate, well,clef = lambda x:int(x.split('_')[2]), filtre=None, tempDir=None, offset=0):
#    def normWrite(img, filename):
#        img=(img-2**15)*(2**8-1)/(2**12-1)
#        vi.writeImage(img, filename)
#        return 1
    # temp directory
    if tempDir is None:
        tempDir = os.path.join(outDir, 'temp')
    # make output directory
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    # movie filename
    movieName = 'P{}_W{}_1.avi'.format(plate, well)
    
    #checking for existence
    if movieName in os.listdir(outDir):
        print 'Done already'
        return
    if filtre == None:
        lstImageNames=os.listdir(imgDir)
    else:
        lstImageNames = filter(filtre, os.listdir(imgDir))
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
    f=open(os.path.join(tempDir, 'list.txt'), 'w')
    
    lstImageNames.sort(key=clef)
    
    for imageName in lstImageNames:
        print imageName
        img = vi.readImage(os.path.join(imgDir, imageName))
        normImage = vigra.VigraArray(img.shape, dtype=np.dtype('uint8'))
        normImage[:,:,0] = (img[:,:,0]-offset)*(2**8-1)/(2**12-1)
###WARNING if you only do normImage = (img - etc then we have a flickering effect. Apparently vigra decides to do its normalization on every image as it pleases
        suffix = imageName.split('.')[-1]
        vi.writeImage(normImage, os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg')), dtype = np.dtype('uint8'))
        f.write(os.path.join(tempDir, os.path.basename(imageName).replace(suffix, 'jpg'))+'\n')    

    f.close()
    # encode command
    encode_command = "mencoder mf://@%s/list.txt -mf fps=3:type=jpg -o %s -ovc xvid -oac copy -xvidencopts fixed_quant=2.5"
    #encode_command = "mencoder mf://@list.txt -mf w=800:h=600:fps=3:type=png -ovc copy -oac copy -o %s"
    encode_command %= (tempDir, os.path.join(outDir, movieName))
    print encode_command
    os.system(encode_command)
    # cleaning up temporary directory
    shell_command = 'rm %s/*' % tempDir
    print shell_command
    os.system(shell_command)
    return

def makeMovieMultiChannels(imgDir, outDir,plate, well, channels=[2,1], tempDir=None, doChannel1Only=True, 
                           ranges=None,secondary = True):
    '''
    From a list of images containing information for two channels in different images, make 
    a movie. It is assumed that the name of a file is something like
     '[conditions]_--W[well]--P0001_t[frame number]_c[channel number].tif'
    
    It is not possible I think that there are more than three channels ? However
    it is possible to have 2 or three channels
    
    Ranges indicate how to do the normalization for the images. It should be [(min_CH1, max_CH1), (min_CH2, max_CH2)]
    300714 90% [(50,750),(300,1000)]
    230714 90% [(150,400),(300,1100)]
    '''

    # movie filename
    movieName = 'P{}_W{}'.format(plate, well)
    
    
    if os.path.isdir(outDir) and movieName+'_1.avi' in os.listdir(outDir):
        print "Done already"
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
    
#If ranges is not provided I look at the min pixel values on all channels for subsequent background substraction
    
    
    if ranges == None:
        min_={ch : 500 for ch in channels}
        max_ = {ch : 0 for ch in channels}
        for imageName in lstImage[channels[0]]:
            img = vi.readImage(os.path.join(imgDir, imageName))
            min_[channels[0]]=min(np.min(img), min_[channels[0]])
            max_[channels[0]]=max(np.max(img), max_[channels[0]])
            
            suffix = imageName.split('_')[-1]
            for ch in channels[1:]:
                imageName2 = os.path.basename(imageName).replace(suffix, 'c{:>05}.tif'.format(ch))
                if secondary:
                    im = vi.readImage(os.path.join(imgDir, imageName2))
                    min_[ch]=min(np.min(im), min_[ch])
                    max_[ch]=max(np.max(im), max_[ch])
    else:
        min_={ch:ranges[channels.index(ch)][0] for ch in channels}
        max_={ch:ranges[channels.index(ch)][1] for ch in channels}
    print min_, max_
    if doChannel1Only:
        #making movie for channel 1 only (nucleus)
        for imageName in lstImage[channels[0]]:
            try:
                img = vi.readImage(os.path.join(imgDir, imageName))
            except OSError:
                normImage = vigra.VigraArray(img.shape, dtype=np.dtype('uint8'))
            else:
                normImage = vigra.VigraArray(img.shape, dtype=np.dtype('uint8'))
    #WARNING if you only do normImage = (img - etc then we have a flickering effect. Apparently vigra decides to do its normalization on every image as it pleases
                im = np.array(img[:,:,0])
                im[np.where(im<min_[channels[0]])]=min_[channels[0]]
                im[np.where(im>max_[channels[0]])]=max_[channels[0]]
                normImage[:,:,0] = (im-min_[channels[0]])*(2**8-1)/(max_[channels[0]]-min_[channels[0]])
                #normImage[np.where(normImage>1)]=1
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
        
    #making movie for both channels if the images exist
    if secondary:    
        for imageName in lstImage[channels[0]]:        
            try:
                img = vi.readImage(os.path.join(imgDir, imageName))
            except OSError:
                colorImage = vigra.VigraArray((shape[0],shape[1], 3), dtype=np.dtype('uint8'))
            else:
                shape = img.shape
                colorImage = vigra.VigraArray((shape[0],shape[1], 3), dtype=np.dtype('uint8'))
                im = np.array(img[:,:,0])
                im[np.where(im<min_[channels[0]])]=min_[channels[0]]
                im[np.where(im>max_[channels[0]])]=max_[channels[0]]
                colorImage[:,:,0] = (im-min_[channels[0]])*(2**8-1)/(max_[channels[0]]-min_[channels[0]])

                suffix = imageName.split('_')[-1]
                for i,ch in enumerate(channels[1:]):
                    imageName2 = os.path.basename(imageName).replace(suffix, 'c{:>05}.tif'.format(ch))
                    im = vi.readImage(os.path.join(imgDir, imageName2))
                    im = np.array(im[:,:,0])
                    im[np.where(im<min_[channels[ch]])]=min_[channels[ch]]
                    im[np.where(im>max_[channels[ch]])]=max_[channels[ch]]
                    colorImage[:,:,i+1] = (im-min_[ch])*(2**8-1)/(max_[ch]-min_[ch])
                    
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