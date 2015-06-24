import os, sys, pdb, getpass
import numpy as np
import cPickle as pickle
from collections import defaultdict
from optparse import OptionParser
from scipy.spatial.distance import squareform, pdist

from util.settings import Settings
from tracking.trajPack.thrivision import thrivisionExtraction
from vigra import impex as vi
from util.listFileManagement import correct_from_Nan, strToTuple

from tracking.histograms import transportation

if getpass.getuser()=='lalil0u':
    import matplotlib.pyplot as p
    from util.plots import couleurs

# 652 Loading error for  LT0159_17--ex2006_01_20--sp2006_01_10--tt17--c5 00134_01
# 654 Loading error for  LT0159_50--ex2006_02_01--sp2006_01_10--tt17--c5 00134_01
# 1150 Loading error for  LT0062_07--ex2006_03_08--sp2005_06_02--tt17--c3 00300_01
# 1151 Loading error for  LT0062_46--ex2005_07_08--sp2005_06_02--tt18--c4 00300_01
# 2675 Loading error for  LT0072_25--ex2005_07_15--sp2005_06_05--tt18--c2 00007_01
# 2676 Loading error for  LT0072_35--ex2007_10_24--sp2005_06_20--tt17--c5 00007_01
# 2677 Loading error for  LT0072_47--ex2005_09_30--sp2005_06_20--tt17--c5 00007_01
# 2678 Loading error for  LT0159_17--ex2006_01_20--sp2006_01_10--tt17--c5 00014_01
# 2679 Loading error for  LT0159_50--ex2006_02_01--sp2006_01_10--tt17--c5 00014_01
# 2698 Loading error for  LT0159_17--ex2006_01_20--sp2006_01_10--tt17--c5 00109_01
# 2699 Loading error for  LT0159_50--ex2006_02_01--sp2006_01_10--tt17--c5 00109_01
# 2865 Loading error for  LT0084_47--ex2005_08_03--sp2005_07_07--tt17--c5 00120_01

#to do jobs launch thrivision.scriptCommand(exp_list, baseName='pheno_seq', command="phenotypes/phenotype_seq.py")

def computingDistance(percentages, who, distance='transport', M=None):
    if len(set(who))!=len(who):
        unred_perc=[]; unred_who=[]
        for i,el in enumerate(who):
            if el not in unred_who:
                unred_perc.append(percentages[i])
                unred_who.append(el)
    
    if distance=='transport':
        result=np.zeros(shape=(percentages.shape[0], percentages.shape[0]))
        for i in range(percentages.shape[0]-1):
            print i,
            result[i, i+1:]=transportation.multSinkhorn(M, lamb=10, r=percentages[i], C=percentages[i+1:].T)
            result[i+1:, i]= result[i,i+1:].T
        return result
    try:
        result=squareform(pdist(percentages, metric=distance))
    except:
        print "Bad distance value"
        return None
    else:
        return result
    


def trajectory_phenotype_comparison(inputFolder, maskFile, inputData):
#     Ordre des phenotypes:
#     Grape
#     increased proliferation
#     Cell Death
#     Binuclear
#     Dynamic changes
#     Large
#     Polylobed
#     Mitotic delay/arrest

#     #i. How were the labels pre-computed?

    #Loading mask
    f=open(os.path.join(inputFolder, maskFile))
    _, mask=pickle.load(f); f.close()
    
    f=open(os.path.join('../data/mitocheck_exp_hitlist_perPheno.pkl'))
    hit_experiments=pickle.load(f); f.close()
    phenos=list(hit_experiments.keys())
    colors = dict(zip(hit_experiments.keys(), couleurs))
    
    new_hit_experiments=[]
    for el in hit_experiments:
        for exp in hit_experiments[el]:
            new_hit_experiments.append((exp, el))
            
    new_hit_experiments=np.array(new_hit_experiments)
    
    f,axes=p.subplots(4,2,sharex=True, sharey=True); k=0
    for i,el in enumerate(new_hit_experiments):
        if i not in mask:
            exp, pheno=el
            axes.flatten()[phenos.index(pheno)].scatter(inputData[k,0],inputData[k,1], color=colors[pheno])
            k+=1
        
    p.show()
    

class pheno_seq_extractor(thrivisionExtraction):
    def __init__(self, setting_file, plate, well):
        super(pheno_seq_extractor, self).__init__(setting_file, plate, well)
        return

    def loadResults(self,exp_list):
        '''
        Here we're loading results on a per experiment basis. This will be interesting to look at distances between experiments
        based on phenotypes, vs distances based on trajectory types.
        '''
        if len(exp_list[0])!=2:
            exp_list=strToTuple(exp_list, os.listdir(self.settings.outputFolder))
        result = None; i=0; missed=[]
        for pl,w in exp_list:
            print i,
            
            try:
                f=open(os.path.join(self.settings.outputFolder,pl, self.settings.outputFile.format(pl[:10], w)), 'r')
                pheno_seq_list, mask = pickle.load(f)
                f.close()
            except:
                print "Loading error for ", pl, w
                missed.append(i)
                continue
            else:
                pheno_seq_list = np.sum( np.array([np.bincount(pheno_seq_list[j], minlength=17) for j in range(len(pheno_seq_list)) if j not in mask]), 0)[:-2]
            #15 and 16 are respectively out of focus and artefact objects. We don't want them
                pheno_seq_list=pheno_seq_list/float(np.sum(pheno_seq_list))
                print pheno_seq_list.shape
                result = np.vstack((result, pheno_seq_list)) if result is not None else pheno_seq_list
            finally:
                i+=1
                
        print "Saving"
        
        f=open(os.path.join(self.settings.outputFolder,self.settings.outputFile.format("ALL", "hit_exp")), 'w')
        pickle.dump((result, missed),f); f.close()
        return
    
    def pheno_seq(self, tracklets,track_filter= (lambda x: x.fusion !=True and len(x.lstPoints)>11)):
        file_=os.path.join(self.settings.hdf5Folder, self.plate, 'hdf5', "{}.hdf5".format(self.well))
        path_objects="/sample/0/plate/{}/experiment/{}/position/1/object/primary__primary".format(self.plate, self.well.split('_')[0])
        path_classif="/sample/0/plate/{}/experiment/{}/position/1/feature/primary__primary/object_classification/prediction".format(self.plate, self.well.split('_')[0])
            
        objects = vi.readHDF5(file_, path_objects)
        classif = vi.readHDF5(file_, path_classif)
        result=[]
        
        for tracklet in filter(track_filter, tracklets):
            t=tracklet.lstPoints.keys(); t.sort(); 
            currR=[]
            for fr, cell_id in t:
                try:
                    where_=np.where((objects['time_idx']==fr)&(objects["obj_label_id"]==cell_id))[0]
                    currR.append(classif[where_])
                except IndexError:
                    currR.append(-1)
            currR = np.array(currR)['label_idx'][:,0]
            result.append(currR)
            
        return result
    
    def getMask(self):
        file_=os.path.join(self.settings.outputFolder, self.plate, self.settings.trajFeaturesFilename.format(self.well))
    
        f=open(file_, 'r')
        arr, _,__= pickle.load(f)
        f.close()
        
        _, toDel = correct_from_Nan(arr, perMovie=False)
        
        return toDel
    
    
    def __call__(self):
            #before anything checking that it passed the qc
        try:
            assert self._usable(check_size=False)
        except AssertionError:
            print "QC failed"
            return
        
        #i.load tracking info
        tracklets, _ = self.load()
        
        #ii. find their bounding boxes
        pheno_sequences=self.pheno_seq(tracklets)
        
        #iii. get the mask for those that are not considered in the feature array
        mask = self.getMask()

        if not os.path.isdir(self.settings.outputFolder):
            os.mkdir(self.settings.outputFolder)
        if not os.path.isdir(os.path.join(self.settings.outputFolder, self.plate)):
            os.mkdir(os.path.join(self.settings.outputFolder, self.plate))
        #iv. crop the boxes in the images
        self.save((pheno_sequences, mask))
        
        return
    
    
if __name__ == '__main__':
    verbose=0
    description =\
'''
%prog - Computing phenotype sequences given tracklets
Input:
- plate, well: experiment of interest. Well under the form 00315_01
- settings file

'''
    #Surreal: apparently I need to do this because I changed the architecture of my code between the computation of the trajectories and now that I need to reopen the files
    from tracking import PyPack
    sys.modules['PyPack.fHacktrack2']=PyPack.fHacktrack2
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    
    parser.add_option("-f", "--settings_file", dest="settings_file", default='phenotypes/settings/settings_pheno_seq.py',
                      help="Settings_file")

    parser.add_option("-p", "--plate", dest="plate",
                      help="The plate which you are interested in")
    
    parser.add_option("-w", "--well", dest="well",
                      help="The well which you are interested in")


    
    (options, args) = parser.parse_args()
    
    p=pheno_seq_extractor(options.settings_file, options.plate, options.well)
    p()
    print "Done"
        
        
