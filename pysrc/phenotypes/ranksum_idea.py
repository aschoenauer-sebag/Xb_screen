import pdb, os, csv, pandas
import numpy as np
import cPickle as pickle
from vigra import impex as vi
from scipy.stats import ranksums

from util import typeD, typeD2, typeD3, ctrl_drugscreen
from _collections import defaultdict

raw_result_dir_Mitocheck= "/share/data40T/Thomas/mitocheck_full_hdf5/out_data"#"/share/data40T/aschoenauer/drug_screen/results_August_2016/mito_joint_classifier"
raw_result_dir_DS= "/share/data40T/aschoenauer/drug_screen/results_August_2016/joint_classifier"

qc_mitocheck_file = '../data/mapping_2014/qc_export.txt'
qc_mitocheckValid_file = '../data/mapping_2014/qc_export_validation.txt'
qc_drugscreen_file = '../data/qc_drugscreen.txt'

test_result_dir = "/share/data40T/aschoenauer/drug_screen/results/mitocheck_tests_pvals"

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

def qc_trsf(plate):
    return plate.split('--')[0]


class Wilcoxon_normalization(object):
    
    def __init__(self,goal="mitocheck", redo=False):
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
    
    @staticmethod
    def loadResults(goal='mitocheck'):
        cols=['Plate', 'Well']
        cols.extend(mitocheck_classes[:15])
        
        test_func = (lambda x: 'LT0900' in x)
        
        files = filter(lambda x: 'CTRL' not in x and not test_func(x), os.listdir(test_result_dir)) if goal=='mitocheck'\
                else filter(lambda x: 'CTRL' not in x and test_func(x), os.listdir(test_result_dir))
        output_=[]
        for el in files:
            f=open(os.path.join(test_result_dir, el))
            d=pickle.load(f); f.close()
            if 'scrambled' in el:
                r=[el.split('_')[0], 'CTRL']
            else:
                r=[el.split('_')[0], el.split('_')[1].split('.')[0]]
            r.extend(d)
            output_.append(r)
            
        return pandas.DataFrame.from_records(output_, columns=cols)
        
    def plateFinder(self):
        if self.goal =='mitocheck':
            plates = os.listdir(self.raw_result_dir)
            plateModels = list(set([el.split('_')[0] for el in plates]))
        else:
            plateModels=['LT0900']
    
        return sorted(plateModels)
    
    def separateControls(self, plates):    
        result={}
        if self.goal=="mitocheck":
            for plateModel in plates:
                if 'Valid' not in plateModel:
                    if int(plateModel[2:])<50:
                        l = list(typeD["scrambled"])
                    else:
                        l = list(typeD2["scrambled"])
                else:
                    l=list(typeD3["scrambled"])
                np.random.shuffle(l)
                result[plateModel] = (l[:6], l[6:])
        else:
            l=list(ctrl_drugscreen)
            np.random.shuffle(l)
            result[plates[0]] = (l[:6], l[6:])
            
        return result
    
    def loadAndTest(self, ctrls):
        
        for plateModel in ctrls:
            if '{}_{}.pkl'.format(plateModel, 'CTRL') in os.listdir(test_result_dir):
                continue
            print 'Doing ', plateModel
            currCtrls, false_exp = ctrls[plateModel]
            
            ctrlData = self.loadData(plateModel, currCtrls)
            
            self.save(ctrlData,plateModel, 'CTRL')
            for i,well in enumerate(false_exp):
                false_expData = self.loadData(plateModel, [well])
                print ctrlData.shape, false_expData.shape
                if false_expData is None or false_expData.shape==():
                    continue
                statList = self.testRankSum(ctrlData, false_expData)
                
                self.save(statList, plateModel, 'scrambled{}'.format(i))
        return
    
    def testRankSum(self, ctrlData, expData):
    
        result=[]
        for k in range(ctrlData.shape[1]):
            result.append(ranksums(ctrlData[:,k], expData[:,k])[1])
            
        return result
        
    
    def loadData(self, plateModel, wellList):
        '''
       This should return a nb wells x nb replicates x 16 matrix of percentages of phenotypes for control wells
'''
        plates = self.plateList[np.where([plateModel in p for p in self.plateList])[0]]
        res = None

        for plate in plates:
            for well in wellList:
                if not well in self.QC[qc_trsf(plate)]:
                    continue
                if not 'hdf5' in os.listdir(os.path.join(self.raw_result_dir, plate))\
                            or not '00{}_01.ch5'.format(well) in os.listdir(os.path.join(self.raw_result_dir, plate, 'hdf5')):
                    print 'No H5 file ', plate, well
                    continue
                filename = os.path.join(self.raw_result_dir, plate, 'hdf5', '00{}_01.ch5'.format(well))
                    
                pathClassif = self.pathClassif.format(plate, '00{}'.format(well))
                try:
                    tabClassification = np.array(vi.readHDF5(filename, pathClassif), dtype=int)
                except ValueError:
                    return None
                
                r = np.bincount(tabClassification, minlength=18)/float(tabClassification.shape[0])
                res = r if res is None else np.vstack((res, r))
                
        #Deleting Out of focus and Artefact nuclei to do the tests
        try:
            res = np.delete(res, [15,16, 17], 1)
        except IndexError:
            try:
                res=res[np.newaxis,:]
                res = np.delete(res, [15,16,17], 1)
            except:
                return None
                
        return res
    
    def loadQC(self):
        self.QC=defaultdict(set)
        
        if self.goal == "mitocheck":
            #QC for labteks 
            f=open(qc_mitocheck_file, 'r')
            reader = csv.reader(f, delimiter='\t'); reader.next()
            
            for el in reader:
                if el[-1]=="ok":
                    plate = el[0].split('--')[0]
                    well = el[0].split('--')[1]
                    self.QC[plate].add(well)
            f.close()
            #QC for validation labteks
            f=open(qc_mitocheckValid_file, 'r')
            reader = csv.reader(f, delimiter='\t'); reader.next()
            
            for el in reader:
                if el[-2]=="True":
                    plate = el[0].split('--')[0]
                    well = el[0].split('--')[1]
                    self.QC[plate].add(well)        
            f.close()
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
                if expData is None or expData.shape==():
                    print "No data for ", experiment
                    continue
                
            #do the test
                stat = self.testRankSum(ctrlData, expData)
            #save the result
                self.save(stat, experiment[0], experiment[1])
                
                