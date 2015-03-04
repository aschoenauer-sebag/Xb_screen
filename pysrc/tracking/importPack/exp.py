'''
To export tracking into h5 files
'''
from util.settings import Settings
import os, pdb
import cPickle as pickle
import vigra.impex as vi
import numpy as np
from _collections import defaultdict

class TrackWriter(object):
    
    TRACK_DTYPE=np.dtype([('obj_idx1', '<u4'), ('obj_idx2', '<u4')])
    OBJECT_PATH_NAME="/sample/0/plate/{}/experiment/{}/position/1/object/{}"
    TRACKING_PATH_NAME="/sample/0/plate/{}/experiment/{}/position/1/object/tracking"

    def __init__(self, plate, well,  setting_file="tracking/settings/settings_track_writer.py", primary_channel_name = 'primary__primary3'):
        '''
        Well is expected as an int
        '''
        self.settings=Settings(setting_file)
        self.plate = plate
        self.well ="{:>05}".format( well)
        self.ch5_file=os.path.join(self.settings.CH5Folder, self.plate, 'hdf5', self.settings.CH5Filename.format(self.well))
        self.primary_channel_name= primary_channel_name
        
    def _getTrajectories(self):
        try:
            f=open(os.path.join(self.settings.trajFolder, self.plate, self.settings.trajFilename.format(self.well)))
            dataDict=pickle.load(f); f.close()
        except Exception as e:
            raise e
        else:
            return dataDict['tracklets dictionary'][self.plate]["{}_01".format(self.well)], dataDict['connexions between tracklets'][self.plate]["{}_01".format(self.well)]
        
    def _getObjectIds(self):
        path = self.OBJECT_PATH_NAME.format(self.plate, self.well, self.primary_channel_name)
        
        try:
            object_ids=vi.readHDF5(self.ch5_file, path)
        except Exception as e:
            raise e
        else:
            return object_ids
        
    def _quality_check(self):
        '''
        Warning:
            H5 file: image number starts at 0
            Browser and annotations: image number starts at 1
        '''
        try:
            f=open(self.settings.image_QC_file)
            l=pickle.load(f)
            f.close()
        except (OSError, IOError) as e:
            raise e
        else:
            try:
                return l[self.plate][int(self.well)]
            except KeyError:
                return []
            
    def _find_id(self, frame, object_frame_id, object_ids):
        try:
            qc_frame_correspondance= frame+len(np.where(self.skipped_frames<=frame)[0])
            assert (len(self.skipped_frames)<=2), "Need to rethink frame correspondences"
            if qc_frame_correspondance in self.skipped_frames:
                qc_frame_correspondance+=1
            
            object_movie_id=np.where((object_ids['time_idx']==qc_frame_correspondance)&(object_ids['obj_label_id']==object_frame_id))[0][0]
        except IndexError:
            print self.plate, self.well,frame, qc_frame_correspondance
            return None
        else:
            return object_movie_id
        
    def _getTracklets(self, tracklets, object_ids):
        '''
        track_id_memory will store for each trajectory with which object_movie_id it started and finished.
        As a consequence, it won't be necessary to go through object_ids again when performing the next step
        '''
        result=[]; track_id_memory=defaultdict(list)
        for track in tracklets.lstTraj:
            points = sorted(track.lstPoints.keys())
            nextId=None
            
            for i,point in enumerate(points):
                try:
                    nextPoint = points[i+1]
                except IndexError:
                    #Meaning we've reached the end of track
                    #Possible cases: either there is only one frame in this tracklet, either not
                    if len(points)==1:
                        currId = self._find_id(frame=point[0], object_frame_id=point[1], object_ids=object_ids)
                        track_id_memory[track.id]=[currId, currId]
                    else:
                        track_id_memory[track.id].append(nextId)
                    continue
                if nextId is None:
                    currId = self._find_id(frame=point[0], object_frame_id=point[1], object_ids=object_ids)
                    track_id_memory[track.id]=[currId]
                else:
                    currId = nextId
                nextId=self._find_id(frame=nextPoint[0], object_frame_id=nextPoint[1], object_ids=object_ids)
                
                result.append((currId, nextId))
                
        return result, track_id_memory
    
    def _getTrackConnexions(self, connexions):
        for frame in connexions:
            currConnexions = connexions[frame]
            for beginning in currConnexions:
                end = currConnexions[beginning]
                
                for track_id in beginning:
                    for track_id2 in end:
                        self.track_result.append((self.track_id_memory[track_id][-1], self.track_id_memory[track_id2][0]))
                    
        return
    
    def _saveResult(self):
        self.track_result=np.array(self.track_result, dtype=self.TRACK_DTYPE)
        path=self.TRACKING_PATH_NAME.format(self.plate, self.well)
        try:
            vi.writeHDF5(self.track_result, self.ch5_file, path)
        except Exception as e:
            raise e
        return 1
    
    def _alreadyDone(self):
        path = self.TRACKING_PATH_NAME.format(self.plate, self.well)
        
        try:
            tracking_table=vi.readHDF5(self.ch5_file, path)
        except:
            return False
        else:
            return True
        
    def __call__(self):
        #if not redo then check if it has already been done
        if not self.settings.redo:
            if self._alreadyDone():
                print "Already done"
                return 1
        
        #i. DATA LOADING
            # read trajectories
        tracklets, connexions = self._getTrajectories()
        
            #also get table of object_id in h5 file
        object_ids= self._getObjectIds()
            #see if some images were deleted during quality control
        self.skipped_frames = self._quality_check()
        
        #ii. get list of connexions for each track
        self.track_result, self.track_id_memory = self._getTracklets(tracklets, object_ids)
        #iii. Complete list of result by adding connexions between tracklets
        self._getTrackConnexions(connexions)
        
        #iv. Save this in the right h5 file
        try:
            self._saveResult()
        except Exception as e:
            print "##############", e.message, self.plate, self.well
        
        
        
        
        
        
        
        