import sys, os, pdb, vigra, getpass
import cPickle as pickle
import numpy as np
import matplotlib as mpl
from scipy.stats.stats import scoreatpercentile
from operator import itemgetter
if getpass.getuser()=='aschoenauer':
    mpl.use('Agg')
import matplotlib.pyplot as p
import brewer2mpl

import vigra.impex as vi
from collections import defaultdict
from util import settings
from util.plots import couleurs, markers
from util.listFileManagement import usable_MITO
from optparse import OptionParser
from itertools import product
    
#ax.imshow(arr, cmap=cm.RdBu_r, interpolation=None)

def comprehensivePlot(exp, inDir, inputFile="cell_cycle_cens_{}.pkl",
                      outputDir = '../resultData/cell_cycle/movies_median', 
                      features=['roisize','total intensity'], 
                      only_complete=False, complete_File = "cell_cycle_info_{}.pkl"):
    '''
    OK. So some trajectories are going to be complete and not be in the second _cens_ file because of their length, which we limited to min 5 frames
    
    '''
    folder = filter(lambda x: x[:9] == exp[:9], os.listdir(inDir))[0]
    try:
        f=open(os.path.join(inDir, folder, inputFile.format(exp[11:]+'_01')))
        d=pickle.load(f)
        f.close()
        
    except:
        print "Pas ", os.path.join(inDir, folder, inputFile.format(exp[11:]+'_01'))
        return
    else:
        filter_=lambda x:True
        if only_complete:
            try:
                f=open(os.path.join(inDir, folder, complete_File.format(exp[11:]+'_01')))
                d_complete=pickle.load(f)
                f.close()
            except:
                print "Not able to do it complete only"
            else:
                filter_=lambda x: x in d_complete['length']
            d={el: {x: d[el][x] for x in filter(filter_, d[el])} for el in d}
        
        num_track = len(d['length']); max_length = np.max(d['length'].values())
        arr=[np.zeros(shape=(num_track, max_length)) for k in range(len(features))]
        for k in range(len(arr)):
            arr[k].fill(-1)
        
        f, axes=p.subplots(10,(num_track/10+1),sharex=True, sharey=True, figsize=(12,12))
        
        for i,el in enumerate(d['roisize']):
            val=d['roisize'][el]
            arr[0][i,:val.shape[0]]=val[:,0,0]/float(val[0,0,0])
            arr[1][i,:val.shape[0]]=d['total intensity'][el][:,0,0]/float(d['total intensity'][el][0,0,0])
            
            if val.shape[0]>=10:
                for k in range(len(arr)):
                    der = [arr[k][i, u+5]- arr[k][i,u] for u in range(val.shape[0]-5)]
                    acc=[der[u+1]-der[u] for u in range(len(der)-1)]

                    try:
                        #Plotting feature
                        axes.flatten()[i].plot(range(val.shape[0]), arr[k][i, :val.shape[0]], label=features[k], color=couleurs[k])
                        axes.flatten()[i].plot(range(len(acc)), acc,color=couleurs[k], linewidth=3)
                        
                        axes.flatten()[i].set_title(el)
                    except:
                        pass
            
        axes.flatten()[i].legend()
                    
        f.savefig(os.path.join(outputDir, '{}_ranges-{}.png'.format(features[0], exp)))
        
        f=p.figure()
        ax=f.add_subplot(221)
        ax.imshow(arr[0], cmap=mpl.cm.RdBu_r, interpolation=None)
        ax.set_yticks(np.linspace(1, num_track, num_track))
        ax.set_yticklabels(d['roisize'].keys(), fontsize=5)
        ax.set_title(features[0])
        
        ax=f.add_subplot(222)
        ax.imshow(arr[1], cmap=mpl.cm.RdBu_r, interpolation=None)
        ax.set_yticks(np.linspace(1, num_track), num_track/5)
        ax.set_title(features[1])
        
        ax=f.add_subplot(212); l=list(d['roisize'].keys())
        ax.set_title(features[0])
        for i,el in enumerate(arr[0]):
            where_=np.where(el>0)[0]
            ax.plot(range(where_.shape[0]), el[where_])
            ax.text(where_[-1], el[where_][-1], l[i])
        
        f.savefig(os.path.join(outputDir, '{}-{}.png'.format(features[0], exp)))
        p.close('all')
        
    return
    
def _filtering_level(liste, possibility, level):
    r1=filter(lambda x: str(possibility)[1:-1]== x[1:5] and 'id0' in x, liste)
    if level=='simple':
        return r1
    increment=lambda x:str(int(x.split('_')[-2][1:])+1)
    r1=filter(lambda x: x[:-10]+increment(x)+'_id2.png' not in liste, r1)
    return r1

def scoreEvaluation(folder, no_folder="pasOK", yes_folder="cell_cycle", plot=True, level='simple'):
    
    possibilities = product(range(3), range(5))
    
    yes=[]; result={}
    no=os.listdir(os.path.join(folder, no_folder))
    for plate in filter(lambda x:'OK' in x, os.listdir(os.path.join(folder, yes_folder))):
        yes.extend(os.listdir(os.path.join(folder, yes_folder, plate)))
    
    for possibility in possibilities:
        currNo=_filtering_level(no, possibility, level)
        currYes = _filtering_level(yes, possibility, level)
        result[possibility]=(len(currNo)+0.0001, len(currYes)+0.0001)
        
    precision = np.array([[  float(100*result[(i,j)][1])/(result[(i,j)][0]+result[(i,j)][1]) for i in range(3)] for j in range(5)])
    count = np.array([[  result[(i,j)][0]+result[(i,j)][1] for i in range(3)] for j in range(5)])-0.0002
    
    if plot:
        print count
        print precision
        cmap = brewer2mpl.get_map('Blues', 'Sequential', 9).mpl_colormap
        f=p.figure()
        
        yticklabels=range(5); yticklabels.reverse()
        xticklabels=range(3)
        xticks = np.array(range(3))
        yticks = np.array(range(5))
        
        ax=f.add_subplot(121)
        ax.pcolormesh(np.flipud(precision), cmap=cmap, vmin=0, vmax=100,edgecolors = 'black')
        ax.set_xticks(xticks+0.5)
        ax.set_xticklabels(xticklabels)
        ax.set_yticks(yticks+0.5)
        ax.set_yticklabels(yticklabels)
        ax.set_xlabel('Score of the first split'); ax.set_ylabel('Score of the second split')
        ax.set_title('Percentage mitosis/splits according to their scores')
        
        ax=f.add_subplot(122)
        ax.pcolormesh(np.flipud(np.array(count, dtype=float)/np.sum(count)), cmap=cmap, vmin=0, vmax=0.4,edgecolors = 'black')
        ax.set_xticks(xticks+0.5)
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel('Score of the first split'); ax.set_ylabel('Score of the second split')
        ax.set_yticks(yticks+0.5)
        ax.set_yticklabels(yticklabels)
        ax.set_title('Repartition of splits according to their scores')
        p.show()
    
    return count,precision
    
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
        self.outputFolder=self.settings.outputFolder

        
        print "Taking tracks that end with movie into account", self.settings.not_ending_track
    
    def _usable(self):
        if not self.settings.new_h5:
            return usable_MITO(self.settings.trackingFolder, [(self.plate, self.well)], self.settings.qc_file, self.settings.mitocheck_file, self.settings.trackingFilename,
                               check_size=False)[0]
        else:
            if self.settings.exist_QC==True:
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
        t,c, m = d['tracklets dictionary'], d['connexions between tracklets'], d['movie_length']
    #NB this dictionary starts with image number 0
    
        self.movie_length = m[self.plate][self.well]
        return t[self.plate][self.well], c[self.plate][self.well]
    
    def _completeConnexions(self, complete_tracks, tracklets, connexions):
        '''
        Here we're interested in getting objects ids when they're complete tracks.
        
        In result we are going to have {object_id : [mom_id, me_id at first frame, me_id at last frame, children id (1,2 or 3) ]
        In siblings we are going to have {object_id : siblings id (1 or 2), appearing also at my first frame}.
        '''
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
        
    def _incompleteConnexions(self, incomplete_tracks, tracklets, connexions):
        '''
        Here we do the same as _completeConnexions except we're not bothering about my children since they don't exist in the video
        '''
        result = {el :[] for el in incomplete_tracks}; siblings={}
        
        c= {im:{el: connexions[im][el] for el in connexions[im] if el[0] in incomplete_tracks or 
                np.any([ww in incomplete_tracks for ww in connexions[im][el]])} for im in connexions}
        c= {im:c[im] for im in filter(lambda x: c[x]!={},c)}
        
        for im in sorted(c):
            for el in c[im]:
                try:
                    traj=filter(lambda x:x.id==el[0], tracklets.lstTraj)[0]
                except IndexError:
                    pdb.set_trace()
                else:
                    sorted_mom_points = sorted(traj.lstPoints.keys(), key=lambda tup:tup[0])
                
                sorted_child_points= {outEl: sorted(filter(lambda x:x.id==outEl, tracklets.lstTraj)[0].lstPoints.keys(), 
                                                    key=lambda tup:tup[0]) for outEl in c[im][el]}
                for outEl in c[im][el]:
                    if el[0] in incomplete_tracks:
                        #we should not have mothers in incomplete tracks as it is precisely the definition of
                        #incomplete tracks that they're not mothers (sorry dudes)
                        pdb.set_trace()
                    if outEl in incomplete_tracks:
                        result[outEl].append(sorted_mom_points[-1])
                        result[outEl].append(sorted_child_points[outEl][0])
                    #ici je mets l'id des siblings de outEl pour pouvoir faire le bon crop
                        siblings[outEl]=[sorted_child_points[sibling][0][1] for sibling in c[im][el] if sibling !=outEl]
        return result, siblings
    
    def findConnexions(self, tracklets, connexions):
        '''
        We need to find the cells for which we have the first appearance after a split, and the disappearance
        After a quick look at the crops, I find that we also need a crop of the cells in the next image
        
        Following Thomas' remark on the 5/29, adding the possibility to also take into account tracks that
        start with a mitosis and end with the end of the movie
        
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
        
    def _findCompleteObjects(self, splits, objects, classification, boxes, compute_boxes, bounding_boxes, siblings, score, frames):
        '''
        Modifies frames and boxes in place to add information about complete tracks
        '''
        for track_id in splits:
            
            local_box=None
            
            for k in range(len(splits[track_id])):
                im, cell_id= splits[track_id][k]
                where_=np.where((objects['time_idx']==im)&(objects['obj_label_id']==cell_id))
                classif = classification[where_]['label_idx'][0]
                #looking at mom or me
                if k==0 or k==2:
                    if compute_boxes:
                        boxes[im].append((track_id, k, bounding_boxes[where_]))
                    if k==0:
                        #looking at mom
                        score[track_id][0]+= int(classif in [6,8,9])
                        frames[track_id][0]=im+1
                    elif k==2:
                        #looking at me being a mom
                        score[track_id][1]+= int(classif in [6,8,9])
                        frames[track_id][1]=im
                elif k==1:
                    #looking at me, being born
                    score[track_id][0]+= int(classif ==7)
                    if compute_boxes:
                        currSibBox = bounding_boxes[where_]
                        for sibling in siblings[track_id]:
                            currSibBox=np.vstack((currSibBox, bounding_boxes[np.where((objects['time_idx']==im)&(objects['obj_label_id']==sibling))]))
                        
                        boxes[im].append((track_id, k, np.array([min(currSibBox['left']), max(currSibBox['right']), min(currSibBox['top']), max(currSibBox['bottom'])])))
                    
                #looking at daughters
                else:
                    if compute_boxes:
                        local_box=bounding_boxes[where_] if local_box is None else np.vstack((local_box, bounding_boxes[where_]))
                    score[track_id][1]+= int(classif ==7)
            if compute_boxes:
                boxes[im].append((track_id,3, np.array([min(local_box['left']), max(local_box['right']), min(local_box['top']), max(local_box['bottom'])]) ))
                
        return
    def _findIncompleteObjects(self, isplits, objects, classification, score, frames):
        '''
        Modifies frames and boxes in place to add information about complete tracks
        '''
        for track_id in isplits:
            for k in range(len(isplits[track_id])):
                im, cell_id= isplits[track_id][k]
                where_=np.where((objects['time_idx']==im)&(objects['obj_label_id']==cell_id))
                classif = classification[where_]['label_idx'][0]
                #looking at mom or me
                if k==0 or k==2:
#                     if compute_boxes:
#                         boxes[im].append((track_id, k, bounding_boxes[where_]))
                    if k==0:
                        #looking at mom
                        score[track_id][0]+= int(classif in [6,8,9])
                        frames[track_id][0]=im+1
                    elif k==2:
                        #looking at me being a mom - SHOULD NOT HAPPEN HERE
#                         score[track_id][1]+= int(classif in [6,8,9])
#                         frames[track_id][1]=im
                        pdb.set_trace()
                elif k==1:
                    #looking at me, being born
                    score[track_id][0]+= int(classif ==7)
#                     if compute_boxes:
#                         currSibBox = bounding_boxes[where_]
#                         for sibling in siblings[track_id]:
#                             currSibBox=np.vstack((currSibBox, bounding_boxes[np.where((objects['time_idx']==im)&(objects['obj_label_id']==sibling))]))
#                         
#                         boxes[im].append((track_id, k, np.array([min(currSibBox['left']), max(currSibBox['right']), min(currSibBox['top']), max(currSibBox['bottom'])])))
                    
                #looking at daughters
                else:#so k=3 4 5, should not happen here
                    pdb.set_trace()
                    
#                     if compute_boxes:
#                         local_box=bounding_boxes[where_] if local_box is None else np.vstack((local_box, bounding_boxes[where_]))
#                     score[track_id][1]+= int(classif ==7)
#             if compute_boxes:
#                 boxes[im].append((track_id,3, np.array([min(local_box['left']), max(local_box['right']), min(local_box['top']), max(local_box['bottom'])]) ))
        return
    
    def findGaleries(self, tracklets):
        if self.settings.new_h5:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.ch5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/{}/object/primary__primary3".format(self.plate, self.well.split('_')[0],
                                                                                                         self.well.split('_')[1])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/{}/feature/primary__primary3/bounding_box".format(self.plate, self.well.split('_')[0],
                                                                                                         self.well.split('_')[1])
        else:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.hdf5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary".format(self.plate, self.well.split('_')[0])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/bounding_box".format(self.plate, self.well.split('_')[0])
            
        objects = vi.readHDF5(file_, path_objects)
        bounding_boxes = vi.readHDF5(file_, path_boundingBox)
        boxes=defaultdict(list)
        
        for track in tracklets:
            lstFrames=sorted(track.lstPoints.keys(), key=itemgetter(0))
            k=0#permits to order the images
            
            for im, cell_id in lstFrames:
                where_=np.where((objects['time_idx']==im)&(objects['obj_label_id']==cell_id))
                boxes[im].append((track.id, k, bounding_boxes[where_]))
                k+=1

        return boxes

    
    def findObjects(self, splits, siblings,isplits, compute_boxes=False):
        if self.settings.new_h5:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.ch5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/{}/object/primary__primary3".format(self.plate, self.well.split('_')[0],
                                                                                            self.well.split('_')[1])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/{}/feature/primary__primary3/bounding_box".format(self.plate, self.well.split('_')[0],
                                                                                            self.well.split('_')[1])
            path_classif="/sample/0/plate/{}/experiment/{}/position/{}/feature/primary__primary3/object_classification/prediction".format(self.plate, self.well.split('_')[0],
                                                                                            self.well.split('_')[1])
            
        else:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.hdf5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary".format(self.plate, self.well.split('_')[0])
            path_boundingBox="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/bounding_box".format(self.plate, self.well.split('_')[0])
            path_classif="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/object_classification/prediction".format(self.plate, self.well.split('_')[0])
            
        objects = vi.readHDF5(file_, path_objects)
        bounding_boxes = vi.readHDF5(file_, path_boundingBox) if compute_boxes else None
        classification = vi.readHDF5(file_, path_classif)
        
        boxes=defaultdict(list)
        score={t:[0,0] for t in splits}
        frames={t:[0,0] for t in splits}

        
        self._findCompleteObjects(splits, objects, classification, boxes, compute_boxes, bounding_boxes, siblings, score, frames)
        
        if self.settings.not_ending_track:
            score.update({t:[0, None] for t in isplits})
            frames.update({t:[0, self.movie_length-1] for t in isplits})            
            self._findIncompleteObjects(isplits, objects, classification, score, frames)
        
        cell_cycle_lengths={el: frames[el][1]-frames[el][0]+1 for el in frames}

        return boxes, cell_cycle_lengths, score
    
    def getObjective(self,objective, lengths, scores, tracklets):
        '''
        This function extracts information dealing with complete tracks, ie starting with a mitosis and ending with a mitosis.
        
        One chooses which information by setting objective in the setting file. This should be under the form :
        objective = {'name':[name of the info to extract], 'function':[function taking as argument objective in case one wants to divide or do smt],
        'features':[list of feature names that are going to be used in objective['function']}
        
        '''
        if self.settings.new_h5:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.ch5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/{}/object/primary__primary3".format(self.plate, self.well.split('_')[0],
                                                                                            self.well.split('_')[1])
            path_features="/sample/0/plate/{}/experiment/{}/position/{}/feature/primary__primary3/object_features".format(self.plate, self.well.split('_')[0],
                                                                                            self.well.split('_')[1])
            path_feature_names = "definition/feature/primary__primary3/object_features"

        else:
            file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.hdf5".format(self.well))
            path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary".format(self.plate, self.well.split('_')[0])
            path_features="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/object_features".format(self.plate, self.well.split('_')[0])
            path_feature_names = "definition/feature/primary__primary/object_features"
            
        objects = vi.readHDF5(file_, path_objects)
        features =vi.readHDF5(file_, path_features)
        
        feature_names = vi.readHDF5(file_, path_feature_names)
        result = {}
        
        if objective in feature_names['name']:
            tab=features[:,np.where(feature_names['name']==objective)[0]]
        else:
            try:
                tabs={el: features[:,np.where(feature_names['name']==el)[0]] for el in objective['features']}
            except:
                raise ValueError
            else:
                tab = objective['function'](tabs)
        
        for track_id in scores:
            if ((scores[track_id][1] is not None and scores[track_id][0]>=1 and scores[track_id][1]>=1) or (scores[track_id][0]>=1 and scores[track_id][1]==None))\
                    and lengths[track_id]>=5:
                try:
                    traj=filter(lambda x:x.id==track_id, tracklets.lstTraj)[0]
                except IndexError:
                    pdb.set_trace()
                else:
                    result[track_id]=self._getObjective(objective, objects, tab, traj)
            else:
                del lengths[track_id]
        if type(objective)==dict:
            return {objective['name']: result, 'length':lengths}
        else:
            return {objective: result, 'length':lengths}
    
    def _getObjective(self, objective, objects, tab, traj):
        result=[]
        
        for im, cell_id in sorted(traj.lstPoints.keys(), key=lambda tup:tup[0]):
            result.append(tab[np.where((objects['time_idx']==im)&(objects['obj_label_id']==cell_id))])
            
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
        if self.settings.outputFile.format(self.well) in os.listdir(os.path.join(self.settings.trackingFolder, self.plate)):
            print "Re-opening existing file"
            f=open(os.path.join(self.settings.trackingFolder, self.plate, self.settings.outputFile.format(self.well)), 'r')
            d=pickle.load(f); f.close()
        else:
            d={}
        d.update(result)
        f=open(os.path.join(self.settings.trackingFolder, self.plate, self.settings.outputFile.format(self.well)), 'w')
        pickle.dump(d, f); f.close()
        
        return
    
    def loadTracks(self):
        try:
            f=open(os.path.join(self.settings.trackingFolder, self.plate, self.settings.outputFile.format(self.well)), 'r')
            d=pickle.load(f); f.close()
        except:
            print "Not able to load pre-existing tracks"
            return -1
        else:
            return d['length'].keys()
    
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
        boxes, _, scores=self.findObjects(splits, siblings, compute_boxes=True)

        if not os.path.isdir(self.outputFolder):
            os.mkdir(self.outputFolder)
        if not os.path.isdir(os.path.join(self.outputFolder, self.plate)):
            os.mkdir(os.path.join(self.outputFolder, self.plate))
        #iv. crop the boxes in the images
        self.crop_to_sequence(boxes, scores)
        self.saveBoxes(boxes)
        
        return
    
    def findObjective(self):
        #before anything checking that it passed the qc
        try:
            assert self._usable()
        except AssertionError:
            print "QC failed"
            return
        
        #i.load tracking info
        tracklets, connexions = self.load()
        
        #ii. find thrivisions
        splits, siblings, isplits, isiblings = self.findConnexions(tracklets, connexions)
        
        #iii. find their bounding boxes
        _, cell_cycle_lengths, scores=self.findObjects(splits, siblings, isplits, compute_boxes=False)
        
        noting_objective = self.getObjective(self.settings.objective, cell_cycle_lengths, scores, tracklets)
        
        self.saveResults(noting_objective)
        return noting_objective['length'].keys()
    
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
    
    parser.add_option("-f", "--settings_file", dest="settings_file", default='tracking/settings/settings_cell_cycle.py',
                      help="Settings_file")

    parser.add_option("-p", "--plate", dest="plate",
                      help="The plate which you are interested in")
    
    parser.add_option("-w", "--well", dest="well",
                      help="The well which you are interested in")
    
    parser.add_option("-o", dest="outputFolder",default=None,
                      help="Output folder")
    
    parser.add_option('--action', type=str, default=None,
                      help='Leave None for finding settings.objective, galerie for extracting images of complete tracks')

    (options, args) = parser.parse_args()
    
    hh=completeTrackExtraction(options.settings_file, options.plate, options.well)
    
    if options.action=='galerie':
        #This will compute galerie images for complete tracks of the plate,well given
        hh.exportGalerieImages(options.outputFolder)
        
    elif options.action=='objective':
        #This will compute length + intensity or roisize of the plate,well given
        hh.findObjective()
    else:
        raise AttributeError('This is not a valid action yet')
    print "Done"
    
    
    
    