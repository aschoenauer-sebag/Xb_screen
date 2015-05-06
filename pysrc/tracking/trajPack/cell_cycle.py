import sys, os, pdb, vigra
import cPickle as pickle
import numpy as np

import vigra.impex as vi
from collections import defaultdict
from util import settings
from util.listFileManagement import usable_MITO
from optparse import OptionParser


class completeTrackExtraction(object):
    '''
    So this class permits to extract the coordinates of complete tracks (ie where we both have beginning after a split
    and end before a split) in an experiment. Furthermore, it gives two scores to all such complete tracks, based on the classification
    of mother and daughters : the first one for the first split, and the second to the second split. With this come crops,
    it enables to see which scores are good enough for us to take the track into consideration.
    '''
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
        We need to find the cells for which we have the first appearance after a split, and the disappearance
        After a quick look at the crops, I find that we also need a crop of the cells in the next image
        '''
        children=[]; mothers=[]
        for im in connexions:
            for el in filter(lambda x: len(connexions[im][x]) in [2,3] and x!=(-1,), connexions[im]):
                mothers.extend(el)
                children.extend(connexions[im][el])
    #ids of complete tracks
        complete_tracks = filter(lambda x: x in mothers, children)
        result = {el :[] for el in complete_tracks}; siblings={}
        
        c= {im:{el: connexions[im][el] for el in connexions[im] if el[0] in complete_tracks or 
                np.any([ww in complete_tracks for ww in connexions[im][el]])} for im in connexions}
        c= {im:c[im] for im in filter(lambda x: c[x]!={},c)}
        
        for im in sorted(c):
            for el in c[im]:
                try:
                    traj=filter(lambda x:x.id==el[0], tracklets.lstTraj)[0]
                except IndexError:
                    pdb.set_trace()
                else:
                    sorted_mom_points = sorted(traj.lstPoints.keys(), key=lambda tup:tup[0])
#                     if el[0] in complete_tracks:
#                         assert(len(result[el[0]])==1)
#                         result[el[0]].extend((sorted_mom_points[0], sorted_mom_points[-1]))
                
                sorted_child_points= {outEl: sorted(filter(lambda x:x.id==outEl, tracklets.lstTraj)[0].lstPoints.keys(), 
                                                    key=lambda tup:tup[0]) for outEl in c[im][el]}
                for outEl in c[im][el]:
                    if el[0] in complete_tracks:
                        result[el[0]].append(sorted_child_points[outEl][0])
                    if outEl in complete_tracks:
                        result[outEl].append(sorted_mom_points[-1])
                        result[outEl].extend((sorted_child_points[outEl][0], sorted_child_points[outEl][-1]))
                    #ici je mets l'id des siblings de outEl pour pouvoir faire le bon crop
                        siblings[outEl]=[sorted_child_points[sibling][0][1] for sibling in c[im][el] if sibling !=outEl]
                            
        return result, siblings
    
    def findObjects(self, splits, siblings):
        if self.settings.new_h5:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.ch5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary3".format(self.plate, self.well.split('_')[0])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary3/bounding_box".format(self.plate, self.well.split('_')[0])
        else:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.hdf5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary".format(self.plate, self.well.split('_')[0])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/bounding_box".format(self.plate, self.well.split('_')[0])
            path_classif="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/object_classification/prediction".format(self.plate, self.well.split('_')[0])
            
        objects = vi.readHDF5(file_, path_objects)
        bounding_boxes = vi.readHDF5(file_, path_boundingBox)
        classification = vi.readHDF5(file_, path_classif)
        
        boxes=defaultdict(list)
        score={t:[0,0] for t in splits}
        
        for track_id in splits:
            
            local_box=None
            
            for k in range(len(splits[track_id])):
                im, cell_id= splits[track_id][k]
                where_=np.where((objects['time_idx']==im)&(objects['obj_label_id']==cell_id))
                classif = classification[where_]['label_idx'][0]
                #looking at mom or me
                if k==0 or k==2:
                    boxes[im].append((track_id, k, bounding_boxes[where_]))
                    if k==0:
                        #looking at mom
                        score[track_id][0]+= int(classif in [6,8,9])
                    elif k==2:
                        #looking at me being a mom
                        score[track_id][1]+= int(classif in [6,8,9])
                elif k==1:
                    #looking at me, being born
                    score[track_id][0]+= int(classif ==7)
                    currSibBox = bounding_boxes[where_]
                    for sibling in siblings[track_id]:
                        currSibBox=np.vstack((currSibBox, bounding_boxes[np.where((objects['time_idx']==im)&(objects['obj_label_id']==sibling))]))
                    
                    boxes[im].append((track_id, k, np.array([min(currSibBox['left']), max(currSibBox['right']), min(currSibBox['top']), max(currSibBox['bottom'])])))
                    
                #looking at daughters
                else:
                    local_box=bounding_boxes[where_] if local_box is None else np.vstack((local_box, bounding_boxes[where_]))
                    score[track_id][1]+= int(classif ==7)
                print score[track_id]
            boxes[im].append((track_id,3, np.array([min(local_box['left']), max(local_box['right']), min(local_box['top']), max(local_box['bottom'])]) ))
            
        return boxes, score
        
    def _findFolder(self):
        if self.settings.new_h5:
            folderName="W{:>05}".format(self.well.split('_')[0])
        else:
            folderName = filter(lambda x: self.well.split('_')[0][2:]==x[:3], os.listdir(os.path.join(self.settings.rawDataFolder, self.plate)))[0]
        return folderName
            
    def _newImageSize(self, crop_):
        '''
        Because if they are on the border otherwise it's not the right sizes for the crop
        Crop is left, right, top, bottom points
        '''
        try:
            X= min(self.settings.XMAX, crop_['right']+self.settings.margin)
            x = max(0, crop_['left']-self.settings.margin)
            x__=int(X-x)
            
            Y = min(self.settings.YMAX, crop_['bottom']+self.settings.margin)
            y = max(0, crop_['top']-self.settings.margin)
            y__=int(Y-y)
            
        except ValueError:
            X= min(self.settings.XMAX, crop_[1]+self.settings.margin)
            x = max(0, crop_[0]-self.settings.margin)
            x__=int(X-x)
            
            Y = min(self.settings.YMAX, crop_[3]+self.settings.margin)
            y = max(0, crop_[2]-self.settings.margin)
            y__=int(Y-y)
            
        return X,x,x__, Y,y,y__
        
    def crop(self, boxes, scores):
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
                              os.path.join(self.settings.outputFolder, self.plate, 
                                           self.settings.outputImage.format(scores[id_], self.plate, self.well.split('_')[0],id_, im,  cell_id)),\
                              dtype=np.dtype('uint8'))
                
        return    
    
    def save(self, boxes):
        f=open(os.path.join(self.settings.outputFolder,self.plate, self.settings.outputFile.format(self.plate[:10], self.well)), 'w')
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
        splits, siblings = self.findConnexions(tracklets, connexions)
        #iii. find their bounding boxes
        boxes, scores=self.findObjects(splits, siblings)

        if not os.path.isdir(self.settings.outputFolder):
            os.mkdir(self.settings.outputFolder)
        if not os.path.isdir(os.path.join(self.settings.outputFolder, self.plate)):
            os.mkdir(os.path.join(self.settings.outputFolder, self.plate))
        #iv. crop the boxes in the images
        self.crop(boxes, scores)
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
    #Surreal: apparently I need to do this because I changed the architecture of my code between the computation of the trajectories and now that I need to reopen the files
    from tracking import PyPack
    sys.modules['PyPack.fHacktrack2']=PyPack.fHacktrack2
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    
    parser.add_option("-f", "--settings_file", dest="settings_file", default='tracking/settings/settings_cell_cycle.py',
                      help="Settings_file")

    parser.add_option("-p", "--plate", dest="plate",
                      help="The plate which you are interested in")
    
    parser.add_option("-w", "--well", dest="well",
                      help="The well which you are interested in")

    (options, args) = parser.parse_args()
    
    hh=completeTrackExtraction(options.settings_file, options.plate, options.well)
    hh()
    print "Done"
    
    
    
    