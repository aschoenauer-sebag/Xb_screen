import pdb, os, csv
import numpy as np
import cPickle as pickle
from vigra import impex as vi
from scipy.stats import ranksums

from util import typeD, typeD2
from _collections import defaultdict

raw_result_dir_Mitocheck= "/share/data40T/aschoenauer/drug_screen/results_August_2016/mito_joint_classifier"
test_result_dir = "/share/data40T/aschoenauer/drug_screen/results/mitocheck_tests"

experimentFilename = '/cbio/donnees/aschoenauer/projects/drug_screen/MITO_experiments.pkl'

plateList = np.array(os.listdir(raw_result_dir_Mitocheck))
primary_channel_name = 'primary__primary3'
pathClassification = "/sample/0/plate/{}/experiment/{}/position/1/feature/%s/object_classification/prediction"%primary_channel_name

class Wilcoxon_normalization(object):
    
    def __init__(self):
        return
    
    def plateFinder(self):
        plates = filter(lambda w: 'Valid' not in w, os.listdir(raw_result_dir_Mitocheck))
        plateModels = list(set([el.split('_')[0] for el in plates]))
    
        return sorted(plateModels)
    
    def separateControls(self, plates):    
        result={}
        
        for plateModel in plates:
            if int(plateModel[2:])<50:
                l = list(typeD["scrambled"])
            else:
                l = list(typeD2["scrambled"])
            np.random.shuffle(l)
            result[plateModel] = (l[:6], l[6:])
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
        plates = plateList[np.where([plateModel in p for p in plateList])[0]]
        res = None
        for plate in plates:
            for well in wellList:
                if not well in self.QC[plate[:9]]:
                    continue
                if not 'hdf5' in os.listdir(os.path.join(raw_result_dir_Mitocheck, plate))\
                            or not '00{}_01.ch5'.format(well) in os.listdir(os.path.join(raw_result_dir_Mitocheck, plate, 'hdf5')):
                    print 'No H5 file ', plate, well
                    continue
                
                filename = os.path.join(raw_result_dir_Mitocheck, plate, 'hdf5', '00{}_01.ch5'.format(well))
                pathClassif = pathClassification.format(plate, '00{}'.format(well))
                tabClassification = np.array(vi.readHDF5(filename, pathClassif), dtype=int)
                
                r = np.bincount(tabClassification, minlength=18)/float(tabClassification.shape[0])
                res = r if res is None else np.vstack((res, r))
                
        #Deleting Out of focus and Artefact nuclei to do the tests
        res = np.delete(res, [15,16], 1)
        return res
    
    def loadQC(self):
        f=open('../data/mapping_2014/qc_export.txt', 'r')
        reader = csv.reader(f, delimiter='\t'); reader.next()
        
        self.QC=defaultdict(set)
        self.siRNAs = defaultdict(set)
        for el in reader:
            #if el[2]!='scrambled':
                
            #siRNA = el[]
            if 'Valid' not in el[0] and int(el[0][2:6])<600 and el[-1]=="ok":
                plate = el[0].split('--')[0]
                well = el[0].split('--')[1]
                self.QC[plate].add(well)
        print "Done loading QC"
        
    def save(self, data, plateModel, name):
        
        f=open(os.path.join(test_result_dir, '{}_{}.pkl'.format(plateModel, name)), 'w')
        pickle.dump(data, f)
        f.close()
        
        return
    
    def findExperimentNames(self):
        f=open(experimentFilename,'r')
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
                expData = self.loadData(experiment[0], [experiment[1]])
                ctrlData = self.loadCtrlData(experiment[0])
            #do the test
                stat = self.testRankSum(ctrlData, expData)
            #save the result
                self.save(stat, experiment[0], experiment[1])
                
                