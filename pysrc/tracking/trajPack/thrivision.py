import os, vigra, pdb, sys

import numpy as np
import cPickle as pickle
import vigra.impex as vi

from collections import defaultdict
from util import settings
from optparse import OptionParser
from util.listFileManagement import usable_MITO
from util import jobSize, progFolder, scriptFolder, path_command, pbsArrayEnvVar, pbsErrDir, pbsOutDir

def scriptThrivision(exp_list, baseName='thrivision'):
    fileNumber = int(len(exp_list)/float(jobSize))+1
    
    head = """#!/bin/sh
cd %s""" %progFolder
    
    for k in range(fileNumber):
        cmd = ''
        for pl, w in exp_list[jobSize*k:jobSize*(k+1)]:
            temp_cmd = """
python tracking/trajPack/division_in_3.py -p %s -w %s"""
            temp_cmd %= (
                        pl,
                        w
                        )
            cmd += temp_cmd
    
        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, '%s%i.sh' % (baseName, k+1))
        script_file = file(script_name, "w")
        script_file.write(head + cmd)
        script_file.close()

        # make the script executable (without this, the cluster node cannot call it)
        os.system('chmod a+x %s' % script_name)
        
        # write the main script
    array_script_name = '%s.sh' % os.path.join(scriptFolder, baseName)
    main_script_file = file(array_script_name, 'w')
    main_content = """#!/bin/sh
%s
#$ -o %s
#$ -e %s
%s$%s.sh
""" % (path_command,
       pbsOutDir,  
       pbsErrDir, 
       os.path.join(scriptFolder, baseName),
       pbsArrayEnvVar)

    main_script_file.write(main_content)
    os.system('chmod a+x %s' % array_script_name)

    # the submission commando is:
    #sub_cmd = 'qsub -o %s -e %s -t 1-%i %s' % (self.oBatchSettings.pbsOutDir,  
    #                                           self.oBatchSettings.pbsErrDir, 
    #                                           jobCount, array_script_name)
    sub_cmd = 'qsub -t 1-%i %s' % (fileNumber, array_script_name)

    print 'array containing %i jobs' % fileNumber
    print sub_cmd
    return 1



class featureExtraction(object):
    def __init__(self, settings_file):
        self.settings=settings.Settings(settings_file, globals())
        
    def _getElements(self, loadingFolders):
        result={}
        for folder in loadingFolders:
            element_list = filter(lambda x: 'crop' in x and int(x.split('_')[-1][2:-4])!=0, os.listdir(folder))
            for el in element_list:
                print el
                decomp=el.split('_')
                cell_id = int(decomp[-1][2:-4])
                frame = int(decomp[-2][1:])
                well = "{}_{}".format(decomp[-4][1:], decomp[-3])
                pl=""
                for x in decomp[1:-4]:
                    pl+="{}_".format(x)
                pl=pl[1:-1]
                
                print pl, well, frame, cell_id
                if pl not in result:
                    result[pl]={well:defaultdict(list)}
                elif well not in result[pl]:
                    result[pl][well]=defaultdict(list)
                result[pl][well][frame].append(cell_id)
            
        return  result
    
    def _getFeatures(self, elements):
        if self.settings.new_h5:
            file_=os.path.join(self.settings.hdf5Folder, "{}", 'hdf5', "{}.ch5")
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary3"
            path_features="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary3/object_features"
        else:
            file_=os.path.join(self.settings.hdf5Folder, "{}", 'hdf5', "{}.hdf5")
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary"
            path_features="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary3/object_features"
        result=None
        for plate in elements:
            for well in elements[plate]:
                objects = vi.readHDF5(file_.format(plate, well), path_objects.format(plate, well))
                features=vi.readHDF5(file_.format(plate, well), path_features.format(plate, well))
                
                for frame in elements[plate][well]:
                    for cell_id in elements[plate][well][frame]:
                        try:
                            line = np.where((objects['time_idx']==frame)&(objects['obj_label_id']))[0]
                        except:
                            pdb.set_trace()
                        else:
                            result=features[line] if result is None else np.vstack((result, features[line]))
                            
        return result
                    
    def _saveResults(self, matrix, filename):
        f=open(os.path.join(self.settings.outputFolder, filename))
        pickle.dump(matrix, f); f.close()
        
        return 1
    
    def __call__(self, loadingFolders=None, filename=None):
        if loadingFolders ==None:
            loadingFolders=[os.path.join(self.settings.outputFolder, plate) for plate in filter(lambda x: os.path.isdir(x) and 'LT' in x, os.listdir(self.settings.outputFolder))]
            
        elements = self._getElements(loadingFolders)
        
        feature_matrix = self._getFeatures(elements)
        if filename is not None:
            self._saveResults(feature_matrix, filename)
        
        return feature_matrix
    
        
class trainingFeatureExtraction(featureExtraction):
    def __init__(self, settings_file):
        
        super(trainingFeatureExtraction, self).__init__(settings_file)
        self.settings_file=settings_file
        
    def __call__(self):
        #i. loading data for positive examples
        extractor = featureExtraction(self.settings_file)
        matrix_TRUE=extractor(loadingFolders = [os.path.join(self.settings.outputFolder, 'True')], filename =None)
        matrix_FALSE=extractor(loadingFolders = [os.path.join(self.settings.outputFolder, 'False')], filename =None)

        matrix_FALSE=np.hstack((matrix_FALSE, np.ones(shape=(matrix_FALSE.shape[0],1))))
        matrix_TRUE=np.hstack((matrix_TRUE, np.ones(shape=(matrix_TRUE.shape[0],1))))
        matrix=np.vstack((matrix_FALSE, matrix_TRUE))
        self._saveResults(matrix, filename=self.settings.outputTrainingFilename)
        
        return 1

class thrivisionExtraction(object):
    def __init__(self, settings_file, plate, well):
        self.settings = settings.Settings(settings_file, globals())
        self.plate = plate
    #NB here the wells are expected in format 00***_01
        self.well = well
    
    def _usable(self):
        if not self.settings.new_h5:
            return usable_MITO(self.settings.trackingFolder, [(self.plate, self.well)], self.settings.qc_file, self.settings.mitocheck_file, self.settings.trackingFilename)[0]
        else:
            f=open(self.settings.qc_file, 'r')
            visual_d=pickle.load(f); f.close()
            
            f=open(self.settings.qc_file2, 'r')
            flou_d=pickle.load(f); f.close()
            
            if self.plate in visual_d and int(self.well.split('_')[0]) in visual_d[self.plate]:
                sys.stderr.write("Visual quality control not passed {} {} \n".format(self.plate, self.well))
                return False   
            if self.plate in flou_d and int(self.well.split('_')[0]) in flou_d[self.plate]:
                sys.stderr.write("Flou quality control not passed {} {} \n".format(self.plate, self.well))
                return False
            return True
            
    def load(self):
        '''
        Here it is important that we also record the tracklet dictionary because that's where the cell ids are, to find back the bounding box in h5 files
        '''
        try:
            f=open(os.path.join(self.settings.trackingFolder, self.plate, self.settings.trackingFilename.format(self.well)))
            d=pickle.load(f); f.close()
        except IOError:
            try:
                f=open(os.path.join(self.settings.trackingFolder, self.plate, self.settings.trackingFilename2.format(self.well)))
                d=pickle.load(f); f.close()
            except IOError:
                print "Non existing trajectory file"
                sys.exit()
        t,c, _ = d['tracklets dictionary'], d['connexions between tracklets'], d['movie_length']
    #NB this dictionary starts with image number 0
        return t[self.plate][self.well], c[self.plate][self.well]
    
    def findConnexions(self, tracklets, connexions):
        '''
        We need to find the last id of the cell that is tracked for linking with h5 file.
        After a quick look at the crops, I find that we also need a crop of the cells in the next image
        '''
        c= {im:{el: connexions[im][el] for el in connexions[im] if len(connexions[im][el])==3} for im in connexions}
        c= {im:c[im] for im in filter(lambda x: c[x]!={},c)}
        
        new_c=defaultdict(list)
        i=0#this is the id for each thrivision in the movie
        for im in sorted(c):
            for el in c[im]:
                try:
                    traj=filter(lambda x:x.id==el[0], tracklets.lstTraj)[0]
                except IndexError:
                    pdb.set_trace()
                else:
                    last_fr, last_id = sorted(traj.lstPoints.keys(), key=lambda tup:tup[0])[-1]
                    assert (last_fr==im)
                    new_c[im].append([i,last_id,-1,-1])
                r=[i]
                for outcomingEl in c[im][el]:
                    try:
                        traj=filter(lambda x:x.id==outcomingEl, tracklets.lstTraj)[0]
                    except IndexError:
                        pdb.set_trace()
                    else:
                        first_fr, first_id = sorted(traj.lstPoints.keys(), key=lambda tup:tup[0])[0]
                        assert (first_fr==im+1)
                        r.append(first_id)
                new_c[im+1].append(r)
                i+=1
        return new_c
    
    def findObjects(self, thrivisions):
        if self.settings.new_h5:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.ch5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary3".format(self.plate, self.well.split('_')[0])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary3/bounding_box".format(self.plate, self.well.split('_')[0])
        else:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.hdf5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary".format(self.plate, self.well.split('_')[0])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/bounding_box".format(self.plate, self.well.split('_')[0])
            
        objects = vi.readHDF5(file_, path_objects)
        bounding_boxes = vi.readHDF5(file_, path_boundingBox)
        boxes=defaultdict(list)
        for im in thrivisions:
            where_=np.where(np.array(objects, dtype=int)==im)[0]
            for thrivision in thrivisions[im]:
                if thrivision[-1]==-1:#we're looking at the mother cell
                    for w in where_:
                        if objects[w][1]==thrivision[1]:
                            boxes[im].append((thrivision[0], thrivision[1], np.array(list(bounding_boxes[w]))))
                            
                else:
                    local_box=np.zeros(shape=(3,4), dtype=int); k=0
                    for w in where_:
                        if np.any(thrivision[1:]==objects[w][1]):
                            local_box[k]=np.array(list(bounding_boxes[w]))
                            k+=1
                    boxes[im].append((thrivision[0],0, np.array([min(local_box[:,0]), max(local_box[:,1]), min(local_box[:,2]), max(local_box[:,3])]) ))
        return boxes
        
    def _findFolder(self):
        if self.settings.new_h5:
            folderName="W{:>05}".format(self.well.split('_')[0])
        else:
            folderName = filter(lambda x: self.well.split('_')[0][2:]==x[:3], os.listdir(os.path.join(self.settings.rawDataFolder, self.plate)))[0]
        return folderName
            
    def _newImageSize(self, crop_):
        '''
        Because if they are on the border otherwise it's not the right sizes for the crop
        '''
        X= min(self.settings.XMAX, crop_[1]+self.settings.margin)
        x = max(0, crop_[0]-self.settings.margin)
        x__=int(X-x)
        
        Y = min(self.settings.YMAX, crop_[3]+self.settings.margin)
        y = max(0, crop_[2]-self.settings.margin)
        y__=int(Y-y)
        
        return X,x,x__, Y,y,y__
        
    def crop(self, boxes):
        '''
        In the crop filename, I add the id of the cell on the image. Otherwise it won't be possible to find it again afterwards,
        for classification purposes
        '''
        folderName=self._findFolder()
        for im in boxes:
            if not self.settings.new_h5:
                #renumbering according to mitocheck image numbering
                local_im=30*im
            else:
                #renumbering according to xb screen image numbering
                local_im=im+1
            image_name=filter(lambda x: self.settings.imageFilename.format(self.well.split('_')[0], local_im) in x, \
                              os.listdir(os.path.join(self.settings.rawDataFolder, self.plate, folderName)))[0]
            image=vi.readImage(os.path.join(self.settings.rawDataFolder, self.plate, folderName, image_name))
            
            for crop_ in boxes[im]:
                id_, cell_id, crop_coordinates = crop_
                X,x,x__, Y,y,y__=self._newImageSize(crop_coordinates)
                
                croppedImage = vigra.VigraArray((x__, y__, 1), dtype=np.dtype('uint8'))
                croppedImage[:,:,0]=(image[x:X, y:Y,0]-self.settings.min_)*(2**8-1)/(self.settings.max_-self.settings.min_)  
                vi.writeImage(croppedImage, \
                              os.path.join(self.settings.outputFolder, self.plate, self.settings.outputImage.format(self.plate, self.well.split('_')[0],id_, im,  cell_id)),\
                              dtype=np.dtype('uint8'))
                
        return    
    
    def save(self, boxes):
        f=open(os.path.join(self.settings.outputFolder,self.plate, self.settings.outputFile.format(self.plate, self.well)), 'w')
        pickle.dump(boxes, f); f.close()
        return
    
    def __call__(self):
    #before anything checking that it passed the qc
        try:
            assert self._usable()
        except AssertionError:
            print "QC failed"
            return

        
        #i.load tracking info
        tracklets, connexions = self.load()
        
        #ii. find thrivisions
        thrivisions = self.findConnexions(tracklets, connexions)
        #iii. find their bounding boxes
        boxes=self.findObjects(thrivisions)

        if not os.path.isdir(self.settings.outputFolder):
            os.mkdir(self.settings.outputFolder)
        if not os.path.isdir(os.path.join(self.settings.outputFolder, self.plate)):
            os.mkdir(os.path.join(self.settings.outputFolder, self.plate))
        #iv. crop the boxes in the images
        self.crop(boxes)
        self.save(boxes)
        
        return
    
    
    
if __name__ == '__main__':
    verbose=0
    description =\
'''
%prog - Finding division in three in an experiment
Input:
- plate, well: experiment of interest
- settings file

'''
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    
    parser.add_option("-f", "--settings_file", dest="settings_file", default='tracking/settings/settings_thrivision.py',
                      help="Settings_file")

    parser.add_option("-p", "--plate", dest="plate",
                      help="The plate which you are interested in")
    
    parser.add_option("-w", "--well", dest="well",
                      help="The well which you are interested in")


    
    (options, args) = parser.parse_args()
    
    thr=thrivisionExtraction(options.settings_file, options.plate, options.well)
    thr()
    print "Done"
    