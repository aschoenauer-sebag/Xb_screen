import csv, os, pdb, getpass
import numpy as np
import cPickle as pickle

from _collections import defaultdict

#Looking at Mitocheck classification using Mitocheck classifier only
#raw_result_dir_Mitocheck= "/share/data40T/Thomas/mitocheck_full_hdf5/out_data"
#primary_channel_name = 'primary__test'
#Looking at Mitocheck classification using the joint classifier
raw_result_dir_Mitocheck= "/share/data40T/aschoenauer/drug_screen/results_August_2016/mito_joint_classifier"
primary_channel_name = 'primary__primary3'

if getpass.getuser()=='lalil0u':
    ds_result_dir = '/media/lalil0u/New/projects/drug_screen/results/'
else:
    ds_result_dir = '/cbio/donnees/aschoenauer/projects/drug_screen/results/'
    from analyzer import interface
    from tracking.importPack.imp import importTargetedFromHDF5

focus_classname = 'OutOfFocus'
artefact_classname = 'Artefact'

def phenotype_aggregated_test(folder='separated_classifier', phenotype="Interphase", choose_ctrls=True, mitocheck='old'):
    pheno_mito = []; pheno_ds =[]
    phenotype = '{}_ch1'.format(phenotype)
    types=[]
    mito_folder= os.path.join(ds_result_dir, 'plates')
    
    if mitocheck=='old':#this means that it is data from the original mitocheck classifier
        filename='mitocheck_ctrls_old'
    else:
        filename='mitocheck_ctrls_new'#data from the joint classifier
    
    if getpass.getuser()!='lalil0u':
        for file in filter(lambda x: filename in x, os.listdir(mito_folder)):
            f=open(os.path.join(mito_folder, file))
            d=pickle.load(f); f.close()
            
            for plate in sorted(d.keys()):
                for well in sorted(d[plate].keys()):
                    try:
                        s= np.sum( d[plate][well][phenotype]*d[plate][well]['object_count'])/float(np.sum(d[plate][well]['object_count']))
                        pheno_mito.append(s)
                        types.append('Mitocheck ctrl')
                    except KeyError:
                        pheno_mito.append(0)
                        print "No {} in mito file".format(phenotype)
                    
    ds_folder = os.path.join(ds_result_dir, folder)#this tells if we're looking at separated or joint classifier but for proliferation 
        #it should not change anything
    if choose_ctrls:
        test = (lambda x: x=='empty')
    else:
        test= lambda x: x!=''
        
    print "Skipping ",
    for el in sorted(os.listdir(ds_folder)):
        f=open(os.path.join(ds_folder, el))
        d=pickle.load(f); f.close()
        
        for well in sorted(filter(lambda x: x!='FAILED QC' and x not in d['FAILED QC'], d.keys())):
            if test(d[well]['Xenobiotic']):
                if phenotype in d[well]:
                    try:
                        s= np.sum(d[well][phenotype]*d[well]['object_count'])/float(np.sum(d[well]['object_count']))
                    except TypeError:
                        try:
                            arr1 = np.array(d[well][phenotype])[:,0]
                        except IndexError:
                            print well, 
                            continue
                        else:
                            s= np.sum(arr1*d[well]['object_count'])/float(np.sum(d[well]['object_count']))
                            pheno_ds.append(s)
                            types.append(d[well]['Xenobiotic'])
                            
                    else:
                        pheno_ds.append(s)
                        types.append(d[well]['Xenobiotic'])
                else:
                    pheno_ds.append(0)
                    types.append(d[well]['Xenobiotic'])
                    
    return np.array(pheno_ds), np.array(pheno_mito), np.array(types)

def phenotype_single_test(folder='separated_classifier', t=0, phenotype="Interphase"):
    pheno_mito = []; pheno_ds =[]
    phenotype = '{}_ch1'.format(phenotype)
    mito_folder= os.path.join(ds_result_dir, 'plates')
    for file in filter(lambda x: 'mitocheck_ctrls' in x, os.listdir(mito_folder)):
        f=open(os.path.join(mito_folder, file))
        d=pickle.load(f); f.close()
        
        for plate in d:
            for well in d[plate]:
                try:
                    pheno_mito.append(d[plate][well][phenotype][t])
                except KeyError:
                    pdb.set_trace()
                    
    ds_folder = os.path.join(ds_result_dir, folder)#this tells if we're looking at separated or joint classifier but for proliferation 
        #it should not change anything
        
    for el in os.listdir(ds_folder):
        f=open(os.path.join(ds_folder, el))
        d=pickle.load(f); f.close()
        
        for well in filter(lambda x: x!='FAILED QC', d.keys()):
            if d[well]['Xenobiotic']=="empty":
                try:
                    pheno_ds.append(d[well][phenotype][t])
                except KeyError:
                    pdb.set_trace()
                    
    return pheno_ds, pheno_mito

def proliferation_test(folder='separated_classifier'):
    prolif_mito = []; prolif_ds =[]
    
    mito_folder= os.path.join(ds_result_dir, 'plates')
    for file in filter(lambda x: 'mitocheck_ctrls' in x, os.listdir(mito_folder)):
        f=open(os.path.join(mito_folder, file))
        d=pickle.load(f); f.close()
        
        for plate in d:
            for well in d[plate]:
                try:
                    prolif_mito.append(d[plate][well]["proliferation"])
                except KeyError:
                    pdb.set_trace()
                    
    ds_folder = os.path.join(ds_result_dir, folder)#this tells if we're looking at separated or joint classifier but for proliferation 
        #it should not change anything
        
    for el in os.listdir(ds_folder):
        f=open(os.path.join(ds_folder, el))
        d=pickle.load(f); f.close()
        
        for well in filter(lambda x: x!='FAILED QC', d.keys()):
            if d[well]['Xenobiotic']=="empty":
                try:
                    prolif_ds.append(d[well]['proliferation'])
                except KeyError:
                    pdb.set_trace()
                    
    return prolif_ds, prolif_mito
    

def extractControlDataFromMitocheck(n=200):
    #i. get n control wells from Mitocheck
    f=open('../data/mapping_2014/qc_export.txt', 'r')
    reader = csv.reader(f, delimiter='\t'); reader.next()
    ctrls = []
    plateList = np.array(os.listdir(raw_result_dir_Mitocheck))
    for el in reader:
        if 'Valid' not in el[0] and int(el[0][2:6])<600 and el[2]=='scrambled' and el[-1]=="ok":
            plate = el[0].split('--')[0]
            well = el[0].split('--')[1]
            truePlate = plateList[np.where([plate in el for el in plateList])[0]][0]
            
            ctrls.append((truePlate, '00{}_01'.format(well)))
    print len(ctrls)
    np.random.shuffle(ctrls)
    
    #ii. load the data under the form {well: data}
    result = loadData(ctrls, n)

    #iii. save it
    f=open('/cbio/donnees/aschoenauer/projects/drug_screen/results/plates/mitocheck_ctrls.pkl', 'w')
    pickle.dump(result, f)
    f.close()
    
    return result
    
def loadData(ctrls, n):
    newFrameLot = None
    dataFolder = raw_result_dir_Mitocheck  
    featureL=['Interphase']
    found=0
    for plate, well in ctrls:
        if found==n:
            break
        filename = os.path.join(dataFolder, plate,"hdf5", '{}.ch5'.format(well))
        try:
            featureL,classes, frameLotC= importTargetedFromHDF5(filename, plate, well,featureL,primary_channel_name=primary_channel_name,
                                                                secondary=False)
        except:
            print "Error at loading from hdf5 ", plate, well
            continue
        else:
            found+=1
        if newFrameLot == None:
            newFrameLot = frameLotC 
        else: 
            newFrameLot.addFrameLot(frameLotC)
            
            featureL = list(featureL)
    if classes is not None:
        featureL.extend(filter(lambda x: x not in featureL, classes))
        
    print'--------------', featureL
    totalResult=defaultdict(dict)
    
    for plate in newFrameLot.lstFrames:
        for well in newFrameLot.lstFrames[plate]:
            result={'object_count':[], 'cell_count':[]}
            
            result.update({'{}_ch1'.format(el):[] for el in featureL})
            
            for frame_nb in newFrameLot.lstFrames[plate][well]:
                frame = newFrameLot.lstFrames[plate][well][frame_nb]
                
            #for each frame, the cell count = object count - [artefacts + out of focus objects]. Only possible if the classification was computed
                if classes is not None:
                    bincount = np.bincount(frame.classification, minlength=len(classes))
                    if np.sum(bincount)==0:
                        print "No info frame {} well {} plate {}".format(frame_nb, well, plate)
                        continue
                    out_of_focus = bincount[np.where(classes==focus_classname)] if focus_classname in classes else 0
#NB: no more smallUnidentified after classifier update on 2/20/15
                    artefacts = bincount[np.where(classes==artefact_classname)] if artefact_classname in classes else 0
                    result["cell_count"].append( frame.centers.shape[0] - (out_of_focus + artefacts))
                    
            #for each frame, the object count = all objects = cells+out of focus objects + artefacts
                result["object_count"].append(frame.centers.shape[0])
                
            #other parameters are specified in the setting file    
                for el in featureL:
                    result['{}_ch1'.format(el)].append(bincount[np.where(classes==el)]/float(np.sum(bincount)))
            
            result["initObjectCount"]=result["object_count"][0]
            result["initCellCount"]=result["cell_count"][0]
            result["proliferation"]=result["cell_count"][-1]/float(result["cell_count"][0])

            for el in classes:
                result['{}_ch1'.format(el)]=np.array(result['{}_ch1'.format(el)])[:,0]
                
            totalResult[plate][int(well[:-3])]=result
            
    return totalResult