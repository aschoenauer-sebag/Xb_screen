import os, sys
import numpy as np
import cPickle as pickle
import vigra.impex as vi
from optparse import OptionParser
import matplotlib.pyplot as p
from operator import itemgetter
from vigra import VigraArray
from skimage import draw

from tracking.PyPack.fHacktrack2 import initXml, ecrireXml, finirXml
from util.settings import Settings
from collections import defaultdict, Counter
from _collections import deque

couleurs=['black', 'orange', 'blue', 'red', 'green', 'purple', 'pink']
y_lim={'nucleus_area':(0,500),
       'cell_area':(0,5000),
       'nuclear_pos_ratio':(0,1),
       'nuclear_axial_ratio':(0,1),
       'cell_axial_ratio':(0,1),
       'nuclear_pos_max':(0,200),
       'nuclear_pos_min':(0,100),
       }

def visualization(setting_file='bgeig/settings/settings_bgeig.py'):
    settings=Settings(setting_file, globals())
    f=open(os.path.join(settings.outputFolder, 'dict_nuclei_sharpmov_10pixelmin.pkl'))
    sharp=pickle.load(f)
    f.close()
    
    f=open(os.path.join(settings.outputFolder, 'dict_track_contact.pkl'))
    contact=pickle.load(f); f.close()
    
    plates=os.listdir(settings.dataFolder)
    wells=['{:>05}_01'.format(k) for k in range(21, 25)]
    
    for plate in plates:
        for well in wells:
            currSharp=sharp[plate][well]
            currContact=contact[plate][well]
            
            f=open(os.path.join(settings.outputFolder, plate, settings.outputFile.format(well+'.ch5')))
            d=pickle.load(f)
            f.close()
            
            for track_id in currSharp:
                print track_id,
                debut = currSharp[track_id][0]
                
                f,axes=p.subplots(2,4,figsize=(24,12))
                
                for i,el in enumerate(d.keys()):
                    axes.flatten()[i].plot(range(debut, debut+len(d[el][track_id])),d[el][track_id], color=couleurs[i], label=el)
                    axes.flatten()[i].set_ylim(y_lim[el])
                    
                    for frame in currContact[track_id]:
                        axes.flatten()[i].axvline(x=frame, color='red', ls='-.', linewidth=3)
                    for s in currSharp[track_id][1:]:
                        axes.flatten()[i].axvline(x=s, color='green', ls='--')
                        
                    axes.flatten()[i].legend()
                axes.flatten()[i].set_title("{} {} {}".format(plate, well, track_id))
                p.savefig(os.path.join(settings.outputFolder, plate, settings.figname.format(well, track_id)))
                p.close('all')
                
    return
                
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
            _, _, hist, _=pickle.load(f); f.close()
        
            for el in hist[measure]:
                result.extend(el)
                
    return result

def findingSharpMovements(setting_file='bgeig/settings/settings_bgeig.py', measure='turningangle', 
                          speed_filter=10, nocontact=True):
    settings=Settings(setting_file, globals())
    
    f=open(settings.dict_corresp_nuclei_cyto)
    dict_corresp=pickle.load(f); f.close()
    
    plates=os.listdir(settings.dataFolder)
    sharp={}; contact={}
    
    for plate in plates:
        print '-----', plate, '<br>'
        wells=['{:>05}_01'.format(k) for k in range(21, 25)]
        sharp[plate]={}; contact[plate]={}
        
        for well in wells:
            print well, '<br>'
            sharp[plate][well]=defaultdict(list)
            contact[plate][well]=defaultdict(list)
            
            contact_count=0; not_contact=0
            
            f=open(os.path.join(settings.outputFolder, plate, settings.feature_filename.format(well)))
            _, coord, hist, id_list=pickle.load(f); f.close()
            
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
                    id_=id_list[i]
                    currcoord=np.array(zip(*coord[i]))
                    sharp[plate][well][id_].append(currcoord[0][0][0])
                    
                    assert(len(el)==len(currcoord)-2)
                    
                    for index in range(len(el)):
                        el=currcoord[index]
                        im, cell_id=el[0]
                        corresp_cytoplasm = dict_corresp[plate][well][im][cell_id]
                        if Counter(dict_corresp[plate][well][im].values())[corresp_cytoplasm]>1:
                            contact[plate][well][id_].append(im)
                            
                        if index in wh_:                        
                            sharp[plate][well][id_].append(im)
                            
                            if nocontact:
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
                                
                                previous_cytoplasm = dict_corresp[plate][well][previous_im][previous_id]
                                next_cytoplasm=dict_corresp[plate][well][next_im][next_id]
                                
                                test = Counter(dict_corresp[plate][well][im].values())[corresp_cytoplasm]==1\
                                and Counter(dict_corresp[plate][well][previous_im].values())[previous_cytoplasm]==1\
                                and Counter(dict_corresp[plate][well][next_im].values())[next_cytoplasm]==1
                            else:
                                test=True
                            
                            if test:
                                d[tuple(el[0])]=(el[1], el[2])
                                not_contact+=1
                            else:
                                contact_count+=1
                            
                    txt, _ = ecrireXml(compteur,d, True)
                    fichierX.write(txt)
                    compteur+=1
            
            fichierX.write(finirXml())
            fichierX.close()
            print 'Cells with sharp movements + contact', contact_count#, '<br>'
            print 'Cells with sharp movements + alone', not_contact#, '<br>'
#So if I'm looking for the track ids that contain at least one sharp movement 
#I should be able to find it in the dictionary sharp which is saved in projects/Geiger/results/dict_nuclei_sharpmov_10pixelmin.pkl
    return sharp, contact

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
                    result[objective['name']][track.id]=self._getObjective(objects, tab, track, 
                                channel=objective['channel'])[:,0,0]
            elif objective['name'] in ['nuclear_pos_min', 'nuclear_pos_max']:
                tabs=defaultdict(dict)
#meaning we're going to actually do calculations on extracted tables
                for el in objective['feature']:
                    if el=='center':
                        
                        objects = vi.readHDF5(file_, path_objects.format('secondary__primary3'))
                        center_=vi.readHDF5(file_, path_centers.format('secondary__primary3'))
                        
                        for track in tracklets:
                            tabs[track.id][el]=self._getObjective(objects, center_, track, channel='secondary__primary3')
                            
                    elif el=='bounding_box':
                        objects = vi.readHDF5(file_, path_objects.format('primary__primary3'))
                        bb_=vi.readHDF5(file_, path_boundingBox.format('primary__primary3'))
                        
                        for track in tracklets:
                            tabs[track.id][el]=self._getObjective(objects, bb_, track, channel='primary__primary3')
                        
                    else:
                        raise ValueError
                
                result[objective['name']] = {el :objective['function'](tabs[el]) for el in sorted(tabs.keys())} 
        
        for track_id in result['nuclear_pos_min']:
            result['nuclear_pos_ratio'][track_id]=result['nuclear_pos_min'][track_id]/result['nuclear_pos_max'][track_id]

        return result
    
    def _getObjective(self, objects, tab, traj, channel):
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
        
#     def _findFolder(self):
#         if self.settings.new_h5:
#             folderName="W{}".format(self.well)
#         else:
#             folderName = filter(lambda x: self.well.split('_')[0][2:]==x[:3], os.listdir(os.path.join(self.settings.rawDataFolder, self.plate)))[0]
#         return folderName
            
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
        
    def crop_to_sequence(self, track_list, currSharp):
        '''
        In the crop filename, I add the number of the image in the galerie so that they're by default ordered in time
        
        I need the bounding box of the membrane and also to use whitefield images.
        '''
        
        for track in track_list:
            id_=track.id
            lstFrames=sorted(track.lstPoints.keys(), key=itemgetter(0))
            rr=[]; cc=[]; val=[]; nextCoord=None
            for i, el in enumerate(lstFrames):
                im, cell_id=el
                coordonnees=track.lstPoints[(im, cell_id)] if nextCoord==None else nextCoord
                try:
                    nextCoord=track.lstPoints[lstFrames[i+1]]
                except IndexError:
                    continue
                else:
                    r,c,v=draw.line_aa(coordonnees[0], coordonnees[1],nextCoord[0], nextCoord[1])
                    rr.extend(r); cc.extend(c); val.extend(v)
                    
            for im, cell_id in lstFrames:
                #renumbering according to xb screen/PCNA image numbering
                local_im=im+1
                
                #draw a dot on the cell which is followed
                cell_x, cell_y=track.lstPoints[(im, cell_id)]
                dot_rr, dot_cc=draw.circle(cell_x, cell_y, radius=3)
                    
                image_name= self.settings.imageFilename.format(self.well, local_im)
                image=vi.readImage(os.path.join(self.settings.allDataFolder, self.plate, 'analyzed', self.well, 'images/tertiary_contours_expanded', image_name))
                
                #X,x,x__, Y,y,y__=self._newImageSize(crop_coordinates)
                x__=self.settings.XMAX; y__=self.settings.YMAX; x=0; y=0; X=self.settings.XMAX; Y=self.settings.YMAX
                
                croppedImage = VigraArray((x__, y__, 3), dtype=np.dtype('float32'))
                croppedImage=image[x:X, y:Y]  
                croppedImage[rr,cc,0]=np.array(val)*255
                
        #If there is a sharp movement, the cell center is pinky red
                if im in currSharp[id_]:
                    croppedImage[dot_rr, dot_cc, 0]=242
                    croppedImage[dot_rr, dot_cc, 1]=21
                    croppedImage[dot_rr, dot_cc, 2]=58
                else:
        #If not, it is green
                    croppedImage[dot_rr, dot_cc, 1]=255
                
                vi.writeImage(croppedImage, \
                              os.path.join(self.outputFolder, self.plate, 'galerie',
                                           self.settings.outputImage.format(self.plate, self.well.split('_')[0],id_, im)),\
                              dtype=np.dtype('uint8'))
                
        return
#     
#     def crop_to_single_image(self, boxes, scores=None):
#         '''
#         In the crop filename, I add the number of the image in the galerie so that they're by default ordered in time
#         '''
#         if scores==None:
#             scores=defaultdict(int)
#             
#         new_boxes=defaultdict(list)
#         
#         folderName=self._findFolder()
#         for im in boxes:
#             if not self.settings.new_h5:
#                 #renumbering according to mitocheck image numbering
#                 local_im=30*im
#                 splitting_index=0
#             else:
#                 #renumbering according to xb screen/PCNA image numbering
#                 local_im=im+1
#                 splitting_index=1
#                 
#             image_name=filter(lambda x: self.settings.imageFilename.format(self.well.split('_')[splitting_index], local_im) in x, \
#                               os.listdir(os.path.join(self.settings.rawDataFolder, self.plate, folderName)))[0]
#             image=vi.readImage(os.path.join(self.settings.rawDataFolder, self.plate, folderName, image_name))
#             
#             for crop_ in boxes[im]:
#                 id_, num, crop_coordinates = crop_
#                 X,x,x__, Y,y,y__=self._newImageSize(crop_coordinates)
#                 
#                 croppedImage = vigra.VigraArray((x__, y__, 1), dtype=np.dtype('uint8'))
#                 croppedImage[:,:,0]=(image[x:X, y:Y,0]-self.settings.min_)*(2**8-1)/(self.settings.max_-self.settings.min_)
#                 
#                 new_boxes[id_].append(croppedImage)
#                 
#         for id_ in new_boxes:
#             y_max = int(np.max([el.shape[1] for el in new_boxes[id_]]))
#             x_size = int(np.sum([el.shape[0] for el in new_boxes[id_]]))
#             newImage = vigra.VigraArray((x_size, y_max, 1), dtype=np.dtype('uint8'))
#             currX=0
#             for croppedImage in new_boxes[id_]:
#                 newImage[currX:currX+croppedImage.shape[0], :croppedImage.shape[1], 0]=croppedImage[:,:,0]
#                 currX+=croppedImage.shape[0]
#         
#             vi.writeImage(newImage, \
#                               os.path.join(self.outputFolder, self.plate, 
#                                            self.settings.outputImage.format(scores[id_], self.plate, self.well.split('_')[0],id_, 0,  0)),\
#                               dtype=np.dtype('uint8'))
#                 
#         return    
    
    
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
        f=open(os.path.join(self.settings.outputFolder, self.settings.sharpMovementFile), 'r')
        sharp=pickle.load(f); f.close()
        
        f=open(os.path.join(self.settings.outputFolder, self.plate, self.settings.traj_filename.format(self.well+'.ch5')), 'r')
        d=pickle.load(f); f.close()
#             
#         except:
#             print "Not able to load pre-existing tracks"
#             return -1
#         else:
        currSharp=sharp[self.plate][self.well]
            
        return currSharp, filter(lambda x: x.id in currSharp, d['tracklets dictionary'][self.plate][self.well].lstTraj)
        
#     
#     def findGaleries(self, tracklets):
#         file_=os.path.join(self.settings.allDataFolder, self.plate, 'hdf5', "{}.ch5".format(self.well))
#         path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary3".format(self.plate, self.well.split('_')[0])
#         path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary3/bounding_box".format(self.plate, self.well.split('_')[0])
# 
#             
#         objects = vi.readHDF5(file_, path_objects)
#         bounding_boxes = vi.readHDF5(file_, path_boundingBox)
#         boxes=defaultdict(list)
#         
#         for track in tracklets:
#             lstFrames=sorted(track.lstPoints.keys(), key=itemgetter(0))
#             k=0#permits to order the images
#             
#             for im, cell_id in lstFrames:
#                 where_=np.where((objects['time_idx']==im)&(objects['obj_label_id']==self.dict_corresp_nuclei_cyto[im][cell_id]))
#                 boxes[im].append((track.id, k, bounding_boxes[where_]))
#                 k+=1
# 
#         return boxes        

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
        currSharp, track_lstPoints = self.loadTracks()

        
        if outputFolder!=None:
            self.outputFolder= outputFolder

        if not os.path.isdir(self.outputFolder):
            os.mkdir(self.outputFolder)
        if not os.path.isdir(os.path.join(self.outputFolder, self.plate)):
            os.mkdir(os.path.join(self.outputFolder, self.plate))
        if not os.path.isdir(os.path.join(self.outputFolder, self.plate, 'galerie')):
            os.mkdir(os.path.join(self.outputFolder, self.plate, 'galerie'))
            
        self.crop_to_sequence(track_lstPoints, currSharp)
        
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




