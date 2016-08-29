import pdb, os, csv
import numpy as np
import cPickle as pickle
from vigra import impex as vi
from scipy.stats import ranksums

from util import typeD, typeD2, ctrl_drugscreen
from _collections import defaultdict

raw_result_dir_Mitocheck= "/share/data40T/Thomas/mitocheck_full_hdf5/out_data"#"/share/data40T/aschoenauer/drug_screen/results_August_2016/mito_joint_classifier"
raw_result_dir_DS= "/share/data40T/aschoenauer/drug_screen/results_August_2016/joint_classifier"

qc_mitocheck_file = '../data/mapping_2014/qc_export.txt'
qc_drugscreen_file = '../data/qc_drugscreen.txt'

test_result_dir = "/share/data40T/aschoenauer/drug_screen/results/mitocheck_tests"

mitocheck_experimentFilename = '/cbio/donnees/aschoenauer/projects/drug_screen/MITO_experiments.pkl'
drugscreen_experimentFilename = '/cbio/donnees/aschoenauer/projects/drug_screen/DS_experiments.pkl'

mitocheck_primary_channel_name = 'primary__test'
DS_primary_channel_name = 'primary__primary3'

pathClassification = "/sample/0/plate/{}/experiment/{}/position/1/feature/%s/object_classification/prediction"

mitocheck_classes = ['Interphase', 'Large', 'Elongated', 'Shape1', 'Shape3', 'Grape',
       'Metaphase', 'Anaphase', 'MetaphaseAlignment', 'Prometaphase',
       'ADCCM', 'Apoptosis', 'Hole', 'Folded', 'SmallIrregular',
       'Artefact', 'UndefinedCondensed', 'OutOfFocus']

DS_classes = ['Interphase', 'Large', 'Elongated', 'Binucleated', 'Polylobed','Grape', 
        'Metaphase', 'Anaphase', 'MetaphaseAlignment','Prometaphase', 
        'ADCCM', 'Apoptosis', 'Hole', 'Folded', 'SmallIrregular', 
        'Artefact', 'Focus', 'Kidney']


class Wilcoxon_normalization(object):
    
    def __init__(self,goal="mitocheck"):
        self.goal = goal
        if self.goal == "mitocheck":
            self.plateList = np.array(os.listdir(raw_result_dir_Mitocheck))
            self.raw_result_dir = raw_result_dir_Mitocheck
            self.pathClassif = pathClassification%mitocheck_primary_channel_name
        else:
            self.plateList = np.array(['LT0900_0{}'.format(k) for k in range(1,4)])
            self.raw_result_dir=raw_result_dir_DS
            self.pathClassif = pathClassification%DS_primary_channel_name
        return
    
    def plateFinder(self):
        if self.goal =='mitocheck':
            plates = filter(lambda w: 'Valid' not in w, os.listdir(self.raw_result_dir))
            plateModels = list(set([el.split('_')[0] for el in plates]))
        else:
            plateModels=['LT0900']
    
        return sorted(plateModels)
    
    def separateControls(self, plates):    
        result={}
        if self.goal=="mitocheck":
            for plateModel in plates:
                if int(plateModel[2:])<50:
                    l = list(typeD["scrambled"])
                else:
                    l = list(typeD2["scrambled"])
                np.random.shuffle(l)
                result[plateModel] = (l[:6], l[6:])
        else:
            l=list(ctrl_drugscreen)
            np.random.shuffle(l)
            result[plates[0]] = (l[:6], l[6:])
            
        return result
    
    def loadAndTest(self, ctrls):
        
        for plateModel in ctrls:
            print 'Doing ', plateModel
            currCtrls, false_exp = ctrls[plateModel]
            
            ctrlData = self.loadData(plateModel, currCtrls)
            
            self.save(ctrlData,plateModel, 'CTRL')
            for i,well in enumerate(false_exp):
                false_expData = self.loadData(plateModel, [well])
                if false_expData.shape==():
                    continue
                statList = self.testRankSum(ctrlData, false_expData)
                
                self.save(statList, plateModel, 'scrambled{}'.format(i))
        return
    
    def testRankSum(self, ctrlData, expData):
    
        result=[]
        for k in range(ctrlData.shape[1]):
            result.append(ranksums(ctrlData[:,k], expData[:,k])[0])
            
        return result
        
    
    def loadData(self, plateModel, wellList):
        '''
       This should return a nb wells x nb replicates x 16 matrix of percentages of phenotypes for control wells
'''
        plates = self.plateList[np.where([plateModel in p for p in self.plateList])[0]]
        res = None
        for plate in plates:
            for well in wellList:
                if not well in self.QC[plate[:9]]:
                    continue
                if not 'hdf5' in os.listdir(os.path.join(self.raw_result_dir, plate))\
                            or not '00{}_01.ch5'.format(well) in os.listdir(os.path.join(self.raw_result_dir, plate, 'hdf5')):
                    print 'No H5 file ', plate, well
                    continue
                filename = os.path.join(self.raw_result_dir, plate, 'hdf5', '00{}_01.ch5'.format(well))
                    
                pathClassif = self.pathClassif.format(plate, '00{}'.format(well))
                tabClassification = np.array(vi.readHDF5(filename, pathClassif), dtype=int)
                
                r = np.bincount(tabClassification, minlength=18)/float(tabClassification.shape[0])
                res = r if res is None else np.vstack((res, r))
                
        #Deleting Out of focus and Artefact nuclei to do the tests
        res = np.delete(res, [15,16, 17], 1)
        return res
    
    def loadQC(self):
        self.QC=defaultdict(set)
        
        if self.goal == "mitocheck":
            f=open(qc_mitocheck_file, 'r')
            reader = csv.reader(f, delimiter='\t'); reader.next()
            
            for el in reader:
                if 'Valid' not in el[0] and int(el[0][2:6])<600 and el[-1]=="ok":
                    plate = el[0].split('--')[0]
                    well = el[0].split('--')[1]
                    self.QC[plate].add(well)
            
            
        else:
            f=open(qc_drugscreen_file, 'r')
            reader = csv.reader(f, delimiter='\t')
            for k in range(1,6):
                self.QC['LT0900_0{}'.format(k)]=['{:0>3}'.format(i) for i in range(1,385)]
            for el in reader:
                plate = el[0].split('--')[0]
                well = '{:0>3}'.format( el[0].split('--')[1])
                self.QC[plate].remove(well)
        f.close()
        print "Done loading QC"
        
    def save(self, data, plateModel, name):
        
        f=open(os.path.join(test_result_dir, '{}_{}.pkl'.format(plateModel, name)), 'w')
        pickle.dump(data, f)
        f.close()
        
        return
    
    def findExperimentNames(self):
        if self.goal=="mitocheck":
            f=open(mitocheck_experimentFilename,'r')
        else:
            f=open(drugscreen_experimentFilename, "r")
            
        l = pickle.load(f)
        f.close()
        
        well_types=set()
        for el in l:
            e=el.split('--')
            well_types.add((e[0][:6], e[1]))
            
        return well_types
    
    def loadCtrlData(self, plateModel):
        f=open(os.path.join(test_result_dir, '{}_CTRL.pkl'.format(plateModel)), 'r')
        d=pickle.load(f)
        f.close()
        
        return d
    
    def __call__(self, work_on_ctrl = True):
        '''
       Goal can be "mitocheck" or "drugscreen".
'''
        #load QC because it is going to be useful before loading data
        self.loadQC()
        
        if work_on_ctrl:
            #i. find plate models
            plates = self.plateFinder()
            
            #ii. for each model, separate into a group of six and one (or two)
            ctrls = self.separateControls(plates)

            #iii. for each plate model, load control data and other data, do tests for control and other and store control data and test result
            self.loadAndTest(ctrls)
            
        else:
            #for each siRNA, note what the plateModel is
            experiments= self.findExperimentNames()
            
            for experiment in experiments:
            #load the data
                try:
                    ctrlData = self.loadCtrlData(experiment[0])
                except:
                    print "No ctrl data ", experiment[0]
                    continue
                
                expData = self.loadData(experiment[0], [experiment[1]])
                if expData.shape==():
                    print "No data for ", experiment
                    continue
                
            #do the test
                stat = self.testRankSum(ctrlData, expData)
            #save the result
                self.save(stat, experiment[0], experiment[1])
                
                