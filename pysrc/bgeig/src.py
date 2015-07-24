import os, sys
import numpy as np
import cPickle as pickle
import vigra.impex as vi
from optparse import OptionParser


from tracking.PyPack.fHacktrack2 import initXml, ecrireXml, finirXml
from util.settings import Settings
from collections import defaultdict, Counter

def imageRenaming(folder):
    currFolder=os.getcwd()
    os.chdir(folder)
    plates=os.listdir(os.getcwd())
    for plate in plates:                               
        wells=os.listdir(plate)    
        for well in wells:                                                      
            for image in filter(lambda x: 'TIF' in x, os.listdir(os.path.join(plate, well))):
                frame = int(image.split('.')[1])
                if 'phase' in image:
                    os.rename(os.path.join(plate, well, image), os.path.join(plate, well, '--W{:>04}--P00001_t{:>05}_c1.tif'.format(int(well[1:]), frame)))
                elif 'calibrate' in image:
                    os.rename(os.path.join(plate, well, image), os.path.join(plate, well, '--W{:>04}--P00001_t{:>05}_c0.tif'.format(int(well[1:]), frame)))
                    
    for plate in plates:                               
        wells=os.listdir(plate)    
        for well in wells:                                                               
            os.rename(os.path.join(plate, well), os.path.join(plate, 'W{:>04}'.format(int(well[1:]))))
    os.chdir(currFolder)
    return

#then for predicting trajectories, writing xml files and computing features:
# for w in ['{:>05}_01.ch5'.format(k) for k in range(21, 25)]:
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic3 -w $w -c 0 --cecog_file 1   
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic2 -w $w -c 0 --cecog_file 1
# for w in ['{:>05}_01.ch5'.format(k) for k in range(21, 25)]:
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic3 -w $w -c 1 --cecog_file 1   
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic2 -w $w -c 1 --cecog_file 1
# for w in ['{:>05}_01.ch5'.format(k) for k in range(21, 25)]:
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic3 -w $w -c 2 --cecog_file 1   
#     %run tracking/trajPack/efficient_traj_pred -f tracking/settings/settings_trajFeatures_Geiger.py -p 14_05_15_Mic2 -w $w -c 2 --cecog_file 1

def findingSharpMovementDistribution(setting_file='bgeig/settings/settings_bgeig.py', measure='turningangle'):
    settings=Settings(setting_file, globals())
    
    plates=os.listdir(settings.dataFolder)
    result=[]
    for plate in plates:
        wells=['{:>05}_01'.format(k) for k in range(21, 25)]
        
        for well in wells:
            f=open(os.path.join(settings.outputFolder, plate, settings.feature_filename.format(well)))
            arr, coord, hist=pickle.load(f); f.close()
        
            for el in hist[measure]:
                result.extend(el)
                
    return result

def findingSharpMovements(setting_file='bgeig/settings/settings_bgeig.py', measure='turningangle', speed_filter=10):
    settings=Settings(setting_file, globals())
    
    f=open(settings.dict_corresp_nuclei_cyto)
    dict_corresp=pickle.load(f); f.close()
    
    plates=os.listdir(settings.dataFolder)
    for plate in plates:
        print '-----', plate, '<br>'
        wells=['{:>05}_01'.format(k) for k in range(21, 25)]
        for well in wells:
            print well, '<br>'
            contact=0
            not_contact=0
            f=open(os.path.join(settings.outputFolder, plate, settings.feature_filename.format(well)))
            arr, coord, hist=pickle.load(f); f.close()
            
            nomFichier = "PL"+plate+"___P"+well+"___T00001.xml"
            nomFichier = os.path.join(settings.outputFolder,'sharp_mov', 'annotations', nomFichier)
            compteur =1
            fichierX = open(nomFichier, "w")
            fichierX.write(initXml())
            
            for i,el in enumerate(hist[measure]):
                if np.any(np.array(el)>0.70):
                    wh_=np.where((np.array(el)>0.70)&(np.array(hist['displacements'][i][:-1])>speed_filter))[0]
                    if wh_.shape[0]==0:
                        continue
                    d={}
                    
                    currcoord=np.array(zip(*coord[i]))

                    for index in wh_:
                        el=currcoord[index]
                        im, cell_id=el[0]
                        try:
                            previous_im, previous_id = currcoord[index-1][0]
                        except IndexError:
                            previous_im=im
                            previous_id=cell_id
                        try:
                            next_im, next_id = currcoord[index+1][0]
                        except IndexError:
                            next_im=im
                            next_id=cell_id
                        
                        corresp_cytoplasm = dict_corresp[plate][well][im][cell_id]
                        previous_cytoplasm = dict_corresp[plate][well][previous_im][previous_id]
                        next_cytoplasm=dict_corresp[plate][well][next_im][next_id]
                        
                        if Counter(dict_corresp[plate][well][im].values())[corresp_cytoplasm]==1\
                            and Counter(dict_corresp[plate][well][previous_im].values())[previous_cytoplasm]==1\
                            and Counter(dict_corresp[plate][well][next_im].values())[next_cytoplasm]==1:
                            d[tuple(el[0])]=(el[1], el[2])
                            not_contact+=1
                        else:
                            contact+=1
                            
                    txt, _ = ecrireXml(compteur,d, True)
                    fichierX.write(txt)
                    compteur+=1
            
            fichierX.write(finirXml())
            fichierX.close()
            print 'Cells with sharp movements + contact', contact, '<br>'
            print 'Cells with sharp movements + alone', not_contact, '<br>'
                
    return

def findingCellsInContact(setting_file='bgeig/settings/settings_bgeig.py'):
    settings=Settings(setting_file, globals())
    
    plates=os.listdir(settings.allDataFolder)
    result={}
    for plate in plates:
        wells=os.listdir(os.path.join(settings.allDataFolder, plate, 'hdf5'))
        result[plate]={}
        for well in wells:
            file_=os.path.join(settings.allDataFolder, plate, 'hdf5', well)
            well=well.split('.')[0]
            result[plate][well]=defaultdict(dict)
            
            path_nuclei_object="/sample/0/plate/{}/experiment/{}/position/{}/object/secondary__primary3".format(plate,well.split('_')[0],
                                                                                            well[-1])
            path_nuclei_centers="/sample/0/plate/{}/experiment/{}/position/{}/feature/secondary__primary3/center".format(plate, well.split('_')[0],
                                                                                            well[-1])
            path_cytoplasmic_object="/sample/0/plate/{}/experiment/{}/position/{}/object/primary__primary3".format(plate,well.split('_')[0],
                                                                                            well[-1])
            path_cytoplasmic_bb="/sample/0/plate/{}/experiment/{}/position/{}/feature/primary__primary3/bounding_box".format(plate, well.split('_')[0],
                                                                                            well[-1])
            
            nuclei_object=vi.readHDF5(file_, path_nuclei_object)
            nuclei_centers = vi.readHDF5(file_, path_nuclei_centers)
            
            cytoplasmic_object=vi.readHDF5(file_, path_cytoplasmic_object)
            cytoplasmic_bb=vi.readHDF5(file_, path_cytoplasmic_bb)
            
            frames=sorted(list(set(nuclei_object['time_idx'])))
            
            for frame in frames:
                cytoplasm_curr_bb = cytoplasmic_bb[np.where(cytoplasmic_object['time_idx']==frame)]
                cytoplasm_curr_id = cytoplasmic_object['obj_label_id'][np.where(cytoplasmic_object['time_idx']==frame)]
                
                for i,bb in enumerate(cytoplasm_curr_bb):
                    wh_=np.where((nuclei_object['time_idx']==frame) & (nuclei_centers['x']>bb['left']) & (nuclei_centers['x']<bb['right'])
                                    & (nuclei_centers['y']>bb['top']) & (nuclei_centers['y']<bb['bottom']))[0]
                    for nuclei_id in nuclei_object['obj_label_id'][wh_]:
                        result[plate][well][frame].update({nuclei_id:cytoplasm_curr_id[i]})
    return result

class geigTrackExtraction(object):
    '''
    So this class permits to extract BGeiger's features with respect to cells, features that we'd like to link to
    sharp changes in cell movement direction.
    '''
    def __init__(self, settings_file, plate, well):
        self.settings = Settings(settings_file, globals())
        self.plate = plate
    #NB here the wells are expected in format 00***_01
        self.well = well
        self.outputFolder=self.settings.outputFolder
        
        f=open(self.settings.dict_corresp_nuclei_cyto)
        self.dict_corresp_nuclei_cyto=pickle.load(f)[plate][well.split('.')[0]]
        f.close()
            
    def load(self):
        '''
        Here it is important that we also record the tracklet dictionary because that's where the cell ids are, to find back the bounding box in h5 files
        '''
        try:
            f=open(os.path.join(self.settings.outputFolder, self.plate, self.settings.traj_filename.format(self.well)))
            d=pickle.load(f); f.close()
        except IOError:
            try:
                f=open(os.path.join(self.settings.outputFolder, self.plate, self.settings.traj_filename2.format(self.well)))
                d=pickle.load(f); f.close()
            except IOError:
                print "Non existing trajectory file"
                sys.exit()
        t,c, m = d['tracklets dictionary'], d['connexions between tracklets'], d['movie_length']
    #NB this dictionary starts with image number 0
    
        #self.movie_length = m[self.plate][self.well]
        return t[self.plate][self.well.split('.')[0]], c[self.plate][self.well.split('.')[0]]
    
#     def _completeConnexions(self, complete_tracks, tracklets, connexions):
#         '''
#         Here we're interested in getting objects ids when they're complete tracks.
#         
#         In result we are going to have {object_id : [mom_id, me_id at first frame, me_id at last frame, children id (1,2 or 3) ]
#         In siblings we are going to have {object_id : siblings id (1 or 2), appearing also at my first frame}.
#         '''
#         result = {el :[] for el in complete_tracks}; siblings={}
#         
#         c= {im:{el: connexions[im][el] for el in connexions[im] if el[0] in complete_tracks or 
#                 np.any([ww in complete_tracks for ww in connexions[im][el]])} for im in connexions}
#         c= {im:c[im] for im in filter(lambda x: c[x]!={},c)}
#         
#         for im in sorted(c):
#             for el in c[im]:
#                 try:
#                     traj=filter(lambda x:x.id==el[0], tracklets.lstTraj)[0]
#                 except IndexError:
#                     pdb.set_trace()
#                 else:
#                     sorted_mom_points = sorted(traj.lstPoints.keys(), key=lambda tup:tup[0])
# #                     if el[0] in complete_tracks:
# #                         assert(len(result[el[0]])==1)
# #                         result[el[0]].extend((sorted_mom_points[0], sorted_mom_points[-1]))
#                 
#                 sorted_child_points= {outEl: sorted(filter(lambda x:x.id==outEl, tracklets.lstTraj)[0].lstPoints.keys(), 
#                                                     key=lambda tup:tup[0]) for outEl in c[im][el]}
#                 for outEl in c[im][el]:
#                     if el[0] in complete_tracks:
#                         result[el[0]].append(sorted_child_points[outEl][0])
#                     if outEl in complete_tracks:
#                         result[outEl].append(sorted_mom_points[-1])
#                         result[outEl].extend((sorted_child_points[outEl][0], sorted_child_points[outEl][-1]))
#                     #ici je mets l'id des siblings de outEl pour pouvoir faire le bon crop
#                         siblings[outEl]=[sorted_child_points[sibling][0][1] for sibling in c[im][el] if sibling !=outEl]
#                             
#         return result, siblings
    
    def findConnexions(self, tracklets, connexions):
        '''
        Maybe we'll use this info, dunno yet
        
        '''
        children=[]; mothers=[]
        for im in connexions:
            for el in filter(lambda x: len(connexions[im][x]) in [2,3] and x!=(-1,), connexions[im]):
                mothers.extend(el)
                children.extend(connexions[im][el])
    #ids of complete tracks
        complete_tracks = filter(lambda x: x in mothers, children)
        c_result, c_siblings = self._completeConnexions(complete_tracks, tracklets, connexions)
        
    #looking at other tracks if we want to, according to settings file
        if self.settings.not_ending_track:
            incomplete_tracks = filter(lambda y: sorted(filter(lambda x:x.id==y, tracklets.lstTraj)[0].lstPoints.keys(), 
                                                        key=lambda tup:tup[0])[-1][0]==self.movie_length-1, children)
            i_result, i_siblings = self._incompleteConnexions(incomplete_tracks, tracklets, connexions)
            
        else:
            i_result, i_siblings=None, None
            
        return c_result, c_siblings, i_result, i_siblings
        
    def getObjective(self,objectives, tracklets):
        '''
        This function extracts information dealing with tracklets.
        
        One chooses which information by setting objective in the setting file. This should be under the form :
        objective = [
            {'name':[name of the info to extract], 
            'channel':[name of the channel on which to extract the info as contained in the hdf5 file],
            'function':[function taking as argument objective in case one wants to divide or do smt],
            'features':[list of feature names that are going to be used in objective['function']}
            ]
        
        '''
        file_=os.path.join(self.settings.allDataFolder, self.plate, 'hdf5', "{}".format(self.well))
        path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/{{}}".format(self.plate, self.well.split('_')[0])
        path_features="/sample/0/plate/{}/experiment/{}/position/1/feature/{{}}/object_features".format(self.plate, self.well.split('_')[0])
        path_feature_names = "definition/feature/primary__primary3/object_features"
        path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/{{}}/bounding_box".format(self.plate, self.well.split('_')[0])
        path_centers="/sample/0/plate/{}/experiment/{}/position/1/feature/{{}}/center".format(self.plate, self.well.split('_')[0])

        feature_names = vi.readHDF5(file_, path_feature_names)
        result = defaultdict(dict)
        
        for objective in filter(lambda x:x['name']!='nuclear_pos_ratio', objectives):
            if 'function' not in objective:
#simply a question of extracting the right table
                objects = vi.readHDF5(file_, path_objects.format(objective['channel']))
                features =vi.readHDF5(file_, path_features.format(objective['channel']))
                tab=features[:,np.where(feature_names['name']==objective['feature'])[0]]
                for track in tracklets:  
                    result[objective['name']][track.id]=self._getObjective(objective, objects, tab, track, 
                                channel=objective['channel'])
            else:
                tabs={}
#meaning we're going to actually do calculations on extracted tables
                for el in objective['feature']:
#                     if el in feature_names['name']:
#                         tabs[el]=features[:,np.where(feature_names['name']==el)[0]]
                    if el=='center':
                        tabs[el]=vi.readHDF5(file_, path_centers.format('secondary__primary3'))
                    elif el=='bounding_box':
                        tabs[el]=vi.readHDF5(file_, path_boundingBox.format('primary__primary3'))
                    else:
                        raise ValueError
                tab = objective['function'](tabs)
        
        for track_id in result['nuclear_pos_min']:
            result['nuclear_pos_ratio'][track_id]=result['nuclear_pos_min'][track_id]/float(result['nuclear_pos_max'][track_id])

        return result
    
    def _getObjective(self, objective, objects, tab, traj, channel):
        result=[]
        
        if channel=='primary__primary3' or channel=='tertiary__expanded':
            #meaning we're looking at the cytoplasm, so we need to look for the right cytoplasm
            for im, cell_id in sorted(traj.lstPoints.keys(), key=lambda tup:tup[0]):
                result.append(tab[np.where((objects['time_idx']==im)&(objects['obj_label_id']==self.dict_corresp_nuclei_cyto[im][cell_id]))])
        elif channel=='secondary__primary3':             
            for im, cell_id in sorted(traj.lstPoints.keys(), key=lambda tup:tup[0]):
                result.append(tab[np.where((objects['time_idx']==im)&(objects['obj_label_id']==cell_id))])
        else:
            raise ValueError
            
        return np.array(result)
        
    def _findFolder(self):
        if self.settings.new_h5:
            folderName="W{}".format(self.well)
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
        
    def crop_to_sequence(self, boxes, scores):
        '''
        In the crop filename, I add the number of the image in the galerie so that they're by default ordered in time
        '''
        if scores==None:
            scores=defaultdict(int)
        
        folderName=self._findFolder()
        for im in boxes:
            if not self.settings.new_h5:
                #renumbering according to mitocheck image numbering
                local_im=30*im
                splitting_index=0
            else:
                #renumbering according to xb screen/PCNA image numbering
                local_im=im+1
                splitting_index=1
                
            image_name=filter(lambda x: self.settings.imageFilename.format(self.well.split('_')[splitting_index], local_im) in x, \
                              os.listdir(os.path.join(self.settings.rawDataFolder, self.plate, folderName)))[0]
            image=vi.readImage(os.path.join(self.settings.rawDataFolder, self.plate, folderName, image_name))
            
            for crop_ in boxes[im]:
                id_, num, crop_coordinates = crop_
                X,x,x__, Y,y,y__=self._newImageSize(crop_coordinates)
                
                croppedImage = vigra.VigraArray((x__, y__, 1), dtype=np.dtype('uint8'))
                croppedImage[:,:,0]=(image[x:X, y:Y,0]-self.settings.min_)*(2**8-1)/(self.settings.max_-self.settings.min_)  
                vi.writeImage(croppedImage, \
                              os.path.join(self.outputFolder, self.plate, 
                                           self.settings.outputImage.format(scores[id_], self.plate, self.well.split('_')[0],id_, im,  num)),\
                              dtype=np.dtype('uint8'))
                
        return
    
    def crop_to_single_image(self, boxes, scores=None):
        '''
        In the crop filename, I add the number of the image in the galerie so that they're by default ordered in time
        '''
        if scores==None:
            scores=defaultdict(int)
            
        new_boxes=defaultdict(list)
        
        folderName=self._findFolder()
        for im in boxes:
            if not self.settings.new_h5:
                #renumbering according to mitocheck image numbering
                local_im=30*im
                splitting_index=0
            else:
                #renumbering according to xb screen/PCNA image numbering
                local_im=im+1
                splitting_index=1
                
            image_name=filter(lambda x: self.settings.imageFilename.format(self.well.split('_')[splitting_index], local_im) in x, \
                              os.listdir(os.path.join(self.settings.rawDataFolder, self.plate, folderName)))[0]
            image=vi.readImage(os.path.join(self.settings.rawDataFolder, self.plate, folderName, image_name))
            
            for crop_ in boxes[im]:
                id_, num, crop_coordinates = crop_
                X,x,x__, Y,y,y__=self._newImageSize(crop_coordinates)
                
                croppedImage = vigra.VigraArray((x__, y__, 1), dtype=np.dtype('uint8'))
                croppedImage[:,:,0]=(image[x:X, y:Y,0]-self.settings.min_)*(2**8-1)/(self.settings.max_-self.settings.min_)
                
                new_boxes[id_].append(croppedImage)
                
        for id_ in new_boxes:
            y_max = int(np.max([el.shape[1] for el in new_boxes[id_]]))
            x_size = int(np.sum([el.shape[0] for el in new_boxes[id_]]))
            newImage = vigra.VigraArray((x_size, y_max, 1), dtype=np.dtype('uint8'))
            currX=0
            for croppedImage in new_boxes[id_]:
                newImage[currX:currX+croppedImage.shape[0], :croppedImage.shape[1], 0]=croppedImage[:,:,0]
                currX+=croppedImage.shape[0]
        
            vi.writeImage(newImage, \
                              os.path.join(self.outputFolder, self.plate, 
                                           self.settings.outputImage.format(scores[id_], self.plate, self.well.split('_')[0],id_, 0,  0)),\
                              dtype=np.dtype('uint8'))
                
        return    
    
    
    def saveBoxes(self, boxes):
        f=open(os.path.join(self.outputFolder,self.plate, self.settings.outputFile.format(self.plate[:10], self.well)), 'w')
        pickle.dump(boxes, f); f.close()
        return
    
    def saveResults(self, result):
        if self.settings.outputFile.format(self.well) in os.listdir(os.path.join(self.settings.outputFolder, self.plate)):
            print "Re-opening existing file"
            f=open(os.path.join(self.settings.outputFolder, self.plate, self.settings.outputFile.format(self.well)), 'r')
            d=pickle.load(f); f.close()
        else:
            d={}
        d.update(result)
        f=open(os.path.join(self.settings.outputFolder, self.plate, self.settings.outputFile.format(self.well)), 'w')
        pickle.dump(d, f); f.close()
        
        return
    
    def loadTracks(self):
        try:
            f=open(os.path.join(self.settings.outputFolder, self.plate, self.settings.outputFile.format(self.well)), 'r')
            d=pickle.load(f); f.close()
        except:
            print "Not able to load pre-existing tracks"
            return -1
        else:
            return d['length'].keys()
    
    def __call__(self):
#before anything there's no qc so we can't check that it was passed
        
        #i.load tracking info
        tracklets, connexions = self.load()
        
        noting_objective = self.getObjective(self.settings.objective, tracklets)
        
        self.saveResults(noting_objective)
        return 
    
    def exportGalerieImages(self, outputFolder):
        '''
        Exporting galerie images.
        '''
        track_ids=-1
        if not self.settings.redo:
            #If fail to open existing track file, will return -1
            track_ids = self.loadTracks()
        if track_ids==-1:
            track_ids = self.findObjective()
            
        tracklets, _ = self.load()
        
        tracklets=filter(lambda x: x.id in track_ids and len(x.lstPoints)>5, tracklets.lstTraj)
        
        if outputFolder!=None:
            self.outputFolder= outputFolder

        boxes=self.findGaleries(tracklets)
        if not os.path.isdir(self.outputFolder):
            os.mkdir(self.outputFolder)
        if not os.path.isdir(os.path.join(self.outputFolder, self.plate)):
            os.mkdir(os.path.join(self.outputFolder, self.plate))
        self.crop_to_single_image(boxes)
        
        return
    
if __name__ == '__main__':
    verbose=0
    description =\
'''
%prog - Finding division in three in an experiment
Input:
- plate, well: experiment of interest. Well name without the suffix hdf5/ch5
- settings file

'''
    #Surreal: apparently I need to do this because I changed the architecture of my code between the computation of the trajectories and now that I need to reopen the files
    from tracking import PyPack
    sys.modules['PyPack.fHacktrack2']=PyPack.fHacktrack2
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    
    parser.add_option("-f", "--settings_file", dest="settings_file", default='bgeig/settings/settings_bgeig.py',
                      help="Settings_file")

    parser.add_option("-p", "--plate", dest="plate",
                      help="The plate which you are interested in")
    
    parser.add_option("-w", "--well", dest="well",
                      help="The well which you are interested in")
    
#     parser.add_option("-o", dest="outputFolder",default=None,
#                       help="Output folder")
#     
#     parser.add_option('--action', type=str, default=None,
#                       help='Leave None for finding settings.objective, galerie for extracting images of complete tracks')

    (options, args) = parser.parse_args()
    
    hh=geigTrackExtraction(options.settings_file, options.plate, options.well)
    hh()
#     
#     if options.action=='galerie':
#         #This will compute galerie images for complete tracks of the plate,well given
#         hh.exportGalerieImages(options.outputFolder)
#         
#     elif options.action=='objective':
#         #This will compute length + intensity or roisize of the plate,well given
#         hh.findObjective()
#     else:
#         raise AttributeError('This is not a valid action yet')
    print "Done"




