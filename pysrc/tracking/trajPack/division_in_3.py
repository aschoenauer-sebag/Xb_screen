import os, vigra

import numpy as np
import cPickle as pickle
import vigra.impex as vi

from collections import defaultdict
from util import settings


class thrivisionExtraction(object):
    def __init__(self, settings_file, plate, well):
        self.settings = settings.Settings(settings_file, globals())
        self.plate = plate
    #NB here the wells are expected in format 00***_01
        self.well = well
        
    def load(self):
        f=open(os.path.join(self.settings.trackingFolder, self.plate, self.settings.trackingFilename.format(self.well)))
        d=pickle.load(f); f.close()
        _,c, _ = d['tracklets dictionary'], d['connexions between tracklets'], d['movie_length']
    #NB this dictionary starts with image number 0
        return c[self.plate][self.well]
    
    def findConnexions(self, connexions):
        c= {im:{el: connexions[im][el] for el in connexions[im] if len(connexions[im][el])==3} for im in connexions}
        return {im:c[im] for im in filter(lambda x: c[x]!={},c)}
    
    def findObjects(self, thrivisions):
        if self.settings.new_h5:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.ch5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary".format(self.plate, self.well.split('_')[0])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/bounding_box".format(self.plate, self.well.split('_')[0])
        else:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.hdf5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary3".format(self.plate, self.well.split('_')[0])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary3/bounding_box".format(self.plate, self.well.split('_')[0])
            
        objects = vi.readHDF5(file_, path_objects)
        bounding_boxes = vi.readHDF5(file_, path_boundingBox)
        boxes=defaultdict(list)
        for im in thrivisions:
            for thrivision in thrivisions[im]:
                where_=np.where((objects[:,0]==im)&(objects[:,1]==thrivision))
                boxes[im].append(bounding_boxes[where_]+self.settings.margin)
                
        return boxes
        
    def crop(self, boxes):
        for im in boxes:
            if self.settings.renumber:
                im=30*im
            image_name=filter(lambda x: self.settings.imageFilename.format(self.well, im) in x, os.listdir(os.path.join(self.settings.rawDataFolder, self.plate)))[0]
            im=vi.readImage(os.path.join(self.settings.rawDataFolder, self.plate, image_name))
            for i,crop_ in enumerate(boxes[im]):
                x__=crop_[1]-crop_[0]; y__=crop_[3]-crop_[2]
                croppedImage = vigra.VigraArray((x__, y__, 1), dtype=np.dtype('uint8'))
                croppedImage[:,:,0]=(im[crop_[0]:crop_[1], crop_[2]:crop_[3],0]-self.settings.min_)*(2**8-1)/(self.settings.max_-self.settings.min_)  
                vi.writeImage(os.path.join(self.settings.outputFolder, self.settings.outputImage.format(self.plate, self.well, im, i)))
                
        return    
    
    def save(self, boxes):
        f=open(os.path.join(self.settings.outputFolder, self.settings.outputFile.format(self.plate, self.well)))
        pickle.dump(boxes, f); f.close()
        return
    
    def __call__(self):
        
        #i.load tracking info
        connexions = self.load()
        
        #ii. find thrivisions
        thrivisions = self.findConnexions(connexions)
        
        #iii. find their bounding boxes
        boxes=self.findObjects(thrivisions)

        if not os.path.isdir(self.settings.outputFolder):
            os.mkdir(self.settings.outputFolder)
        #iv. crop the boxes in the images
        self.crop(boxes)
        self.save(boxes)
        
        return