import os, pdb, shutil
from collections import defaultdict
from operator import itemgetter
from collections import Counter
import cPickle as pickle
import numpy as np
import vigra.impex as vi
from PIL import Image

from util import typeD, typeD2
from tracking.trajPack import sdFeatures1, logTrsf

def is_ctrl(experiment):
    id_ = int(experiment[0].split('_')[0][2:])
    if id_<50:
        if np.any(np.array([experiment[1][2:5] in l for l in typeD.values()])):
            return True
    else:
        if np.any(np.array([experiment[1][2:5] in l for l in typeD2.values()])):
            return True
    return False


def renameFromZeiss(inputFolder, plate, outputFolder=None):        
#ici c'est pour renommer des images dans le cas ou pour une raison connue de Dieu seul
#le puits n'est pas passe dans le nom du fichier
#inputfolder doit etre le nom du dossier avec la plaque
#    for well in os.listdir(inputFolder):
#        for file_ in os.listdir(os.path.join(inputFolder,well)):
#            os.rename(os.path.join(inputFolder,well, file_), os.path.join(inputFolder,well, file_[:-18]+'_--'+well+'--P00001'+file_[-18:]))

    #SI ON A DES IMAGES ZEISS PUREMENT ET SIMPLEMENT, noms des dossiers deja changes
     for well in range(1,97):
        w = "W%05i"%well
        for el in os.listdir(os.path.join(inputFolder, w)):
            timep = int(el.split('_')[2][1:4])
            ch = int(el.split('_')[2][-1])
            os.rename(os.path.join(inputFolder, w, el), os.path.join(inputFolder,w, '%s--%s--P0001_t%05i_c%05i.tif'%(plate, w, timep, ch)))


'''
File containing scripts to deal with files: copying from a place to another, counting how many there are of each type etc.
'''

def deletingAllForPlateWell(todel, folder='/share/data20T/mitocheck/tracking_results'):
    '''
    To delete all files for a list of (plate, well) elements in the tracking results folder.
    
    Input:
    - todel: list of (plate, well) elements
    - folder: folder where to delete files
    '''
    
    baseNames2 = ['centers_P{}_w{}.pkl', "avcum_P{}_{}.pkl"]
    base1=['traj_noF_densities_w{}.hdf5.pkl', 'TERd_tabFeatures_{}.pkl','d_tabFeatures_{}.pkl']
    for el in todel:
        print el
        ll=el.split('/')
        plate = ll[-3]; well=ll[-1][:-5]
        for name in baseNames2:
            try:
                os.remove(os.path.join(folder, plate, name.format(plate, well)))
            except OSError:
                print "This file {} does not exist".format(os.path.join(folder, plate, name.format(plate, well)))
        for name in base1:
            try:
                os.remove(os.path.join(folder, plate, name.format(well)))
            except OSError:
                print "This file {} does not exist".format(os.path.join(folder, plate, name.format(well)))
        
    return

def countingHDF5(filename, plate, well):
    '''
    Counting the number of images contained in an hdf5 file as a proxy for the number of images of the original experiment
    
    Input:
    - filename: complete path to the hdf5 file as produced by CellCognition when extracting features from the raw data
    - 
    '''

    pathObjects = "/sample/0/plate/"+plate+"/experiment/"+well[:-3]+"/position/"+well[-1]+"/object/primary__primary"
    try:
        tabObjects = vi.readHDF5(filename, pathObjects)
    except:
        print 'No file'
        return 1000
    else:
        #changed to take into account that for some movies frames are not in the chronological order in the hdf5 file
        return np.max(np.array(tabObjects, dtype=int))+1#otherwise we forget the last frame

def checkingAllHDF5(folder="/share/data20T/mitocheck/Alice/results",
    folderRaw='/share/data20T/mitocheck/compressed_data', liste=None):
    '''
    Checking the number of images for all experiments in a certain folder.
    The function checks the number of images in the raw data folder as well as the number of images
    in the hdf5 file. If any of those is under 90, it is recorded. A table summing up the results is saved
    
    Input:
    - folder: location of hdf5 files 
    - folderRaw: location of the raw data
    
    Output:
    BAD: experiments where nbImages<90 in hdf5 file but not in raw data
    veryBAD: experiments where nbImages<90 in both hdf5 file and raw data
    + saves a file hdf5ToDel.pkl with both outputs
    '''    
    
    result_arr=[]
    if liste == None:
        liste=os.listdir(folderRaw)
    liste.sort()

    for plate in liste:
        print plate
        fileList = os.listdir(os.path.join(folderRaw, plate))
        fileList.sort()
        for fichier in fileList:
            well=fichier[:3]

            nbImages = len(os.listdir(os.path.join(folderRaw, plate, fichier)))

            filename = os.path.join(folder, plate, 'hdf5', '00{}_01.hdf5'.format(well))
            number = countingHDF5(filename, plate, '00{}_01'.format(well))
                        
            result_arr.append((plate, well, nbImages, number))
    f=open('nb_Images2_raw_hdf5.pkl', 'w'); pickle.dump(result_arr, f); f.close()
    return 

def countingPNG(folder='/share/data20T/mitocheck/compressed_data'):
    '''
    Counts the number of PNG images of the raw data in the given folder
    
    Input:
    -folder: where the raw data is located
    
    Output:
    saves a file PNG.pkl with the experiments with png images as opposed to tiff, in the raw data
    '''
    liste=os.listdir(folder); liste.sort()
    print len(liste)
    small=[]
    for plate in liste:
        print "'",plate, "',"
        ll=os.listdir(os.path.join(folder, plate))
        for well in ll:
            images = os.listdir(os.path.join(folder, plate, well))
            zou = np.any(np.array([image[-3:]=="png" for image in images]))
            if zou:
                print '******************',plate, well
                small.append(os.path.join(folder, plate, well))
                
    f=open('PNG.pkl', 'w'); pickle.dump(small, f); f.close()
    return
                    

def copying(folder, bigfolder):
    '''
    Copies all hdf5 files in folder to bigfolder, keeping the same folder architecture
    as produced by CellCognition
    
    Input:
    -folder: where to copy from
    - bigfolder: where to copy to
    '''
    depart = os.listdir(folder)
    arrivee=os.listdir(bigfolder)
    for pl in depart:
        print pl, "au depart"
        if pl not in arrivee:
            shutil.copytree(os.path.join(folder,pl), os.path.join(bigfolder, pl))
            print 'copying all', pl
        elif 'hdf5' not in os.listdir(os.path.join(bigfolder, pl)) and 'hdf5' in os.listdir(os.path.join(folder, pl)):
            shutil.copytree(os.path.join(folder, pl, 'hdf5'), os.path.join(bigfolder, pl, 'hdf5'))
            print 'copying hdf5', pl
        elif 'hdf5' in os.listdir(os.path.join(folder, pl)):
            for w in os.listdir(os.path.join(folder, pl, 'hdf5')):
                if w not in os.listdir(os.path.join(bigfolder, pl, 'hdf5')):
                    shutil.copy(os.path.join(folder, pl, 'hdf5', w), os.path.join(bigfolder, pl, 'hdf5'))
                    print 'copying ', pl, w
        else:
            print "hdf5 not in ", os.path.join(folder, pl)
    return

def deletingDoublets(folderDoub, folderToKeep):
    '''
    Deletes files that are both in folderDoub and folderToKeep from folderDoub
    
    Input:
    -folderDoub: where the doublet files will be deleted
    -folderToKeep: where the files are kept
    
    Output:
    nb of deletions
    '''
    countr=0
    for pl in os.listdir(folderDoub):
        try:
            wellCour = os.listdir(os.path.join(folderDoub, pl, 'hdf5'))
        except:
            pass
        else:
            if pl in os.listdir(folderToKeep):
                try:
                    lp = os.listdir(os.path.join(folderToKeep, pl, 'hdf5'))
                except:
                    pass
                else:
                    for well in lp:
                        if well in wellCour:
                            print pl, well
                            countr+=1
                            os.remove(os.path.join(folderDoub, pl, 'hdf5', well))
    return countr

def txtToList(fichier):
    '''
    Making an array from a text file that contains \r and \t and \n, without giving the first line
    
    Input: filename with whole path
    Ouput: array with text that was in the file
    '''
    f=open(fichier, 'r');
    lines = [x.rstrip().strip('\n').split('\t') for x in f.readlines()]; f.close()
    #print lines[0]
    lines=np.array(lines[1:])
    return lines

def geneListToFile(genes, name):
    '''
    Writing a gene list into a text file
    
    Input: 
    - genes: gene list
    - name: filename, will be written in working directory
    '''
    f=open(name, 'w')
    for k in range(len(genes)):
        f.write('{}\t \n'.format(genes[k]))
    f.close()    
    print 'fichier ecrit ', os.path.join(os.getcwd(), name)
    return 1

def writeGSEARankingFile(genes, name):
    if name[-3:]!='rnk':
        print 'Warning you need to put rkn extension for GSEA'
        return
    #genes=sorted(genes, key=itemgetter(1))
    f=open(name, 'w')
    f.write('#Name \t Mitocheck \n')
    for k in range(len(genes)):
        f.write('{}\t {} \n'.format(genes[k][0], genes[k][1]))
    f.close()    
    print 'fichier ecrit ', os.path.join(os.getcwd(), name)
    return 1

def strToTuple(strList, platesList):
    '''
    Goes from an experiment list in the form LT00*****--00*** to 
    the whole name as recorded in the mitocheck raw data folder
    
    Input:
    - strList: initial experiment list
    - platesList: list of plates in the raw data
    
    Output: experiment list in the form usable to find experiments in the data
    '''
    result = []
#    f=open(fichier, "r")
#    platesList = pickle.load(f); f.close()
    
    for el in strList:
        try:
            pl = filter(lambda x: el[:9] in x, platesList)[0]
        except IndexError:
            print el, 'not found'
        else:
            result.append((pl, '{:>05}_01'.format(el[-3:])))
    return result
    
def noRaw(arrExp, moviesD='/share/data20T/mitocheck/compressed_data'):
    '''
    Checks the existence of raw data for a given experiment list
    
    Input:
    - arrExp: experiment list
    - moviesD: raw data folder
    
    Output:
    -bou: number of experiments that are validation experiments (hence not yet taken into account)
    - resultNoRaw: experiments for which raw data is not in raw data folder
    - resultToDo: experiments for which raw data exists
    '''
    
    resultNoRaw=[]; resultToDo=[]; bou=0
    arrExp=np.array(arrExp)
    for k in range(arrExp.shape[0]):
        if 'LTValid' in arrExp[k,0]:
            bou+=1
        else:
            l=os.listdir(os.path.join(moviesD, arrExp[k,0]))
            zz=filter(lambda x: arrExp[k,1][2:5]==x[:3], l)
            if len(zz)==0:resultNoRaw.append(arrExp[k])
            elif len(zz)>1: raise
            else:
                resultToDo.append(arrExp[k])
    return bou, resultNoRaw, resultToDo    

def countingDone(experiments,featlistonly=True, name=None,rawD='/share/data20T/mitocheck/Alice/results',\
                 trajD = "/share/data20T/mitocheck/tracking_results", featureTabName='hist2_tabFeatures_{}.pkl'):
    '''
    Counts of a list of experiments how many have been dealt with with regard to object feature extraction,
    cell trajectory extraction and finally trajectory feature extraction.
    
    Input:
    - experiments: list of experiments
    - featlistonly: if the output should only be list of experiments for which trajectory features have been extracted
    - name: in case the output should be saved in a file
    - rawD: folder where the raw data is
    - trajD: folder where the extracted trajectories are
    
    Output:
    if featlistonly: list of experiments for which trajectories have been extracted
    else: lists of experiments for which object features, cell trajectories and trajectory features
    should respectively be extracted
    '''


    baseNames = ['{}.hdf5', 'traj_noF_densities_w{}.hdf5.pkl', featureTabName]
    print "experiments : ", len(experiments)
    no_hdf5=[]
    no_tracking=[]; no_trajfeat =[]
    counthdf5=0; counttracking = 0; countfeat=[]
    
    for line in experiments:
        for baseName in baseNames:
            if baseName == '{}.hdf5':

                if line[0] not in os.listdir(rawD) or 'hdf5' not in os.listdir(os.path.join(rawD, line[0])) or baseName.format(line[1]) not in os.listdir(os.path.join(rawD, line[0], 'hdf5')): 
                    no_hdf5.append(line)
                else:
                    counthdf5+=1
                    
            elif baseName == 'traj_noF_densities_w{}.hdf5.pkl':
                if line[0] not in os.listdir(os.path.join(trajD)) or baseName.format(line[1]) not in os.listdir(os.path.join(trajD, line[0])):
                    no_tracking.append(line)
                else:
                    counttracking+=1
            elif baseName == featureTabName:
                if line[0] not in os.listdir(os.path.join(trajD)) or baseName.format(line[1]) not in os.listdir(os.path.join(trajD, line[0])):
                    no_trajfeat.append(line)
                else:
                    countfeat.append(line)
                    
    no_hdf5=np.array(no_hdf5); no_tracking=np.array(no_tracking)
    print len(no_hdf5)
    lt, noraw, todoraw = noRaw(no_hdf5)
#    
    print noraw
#    f=open("/cbio/donnees/aschoenauer/workspace2/Tracking/{}.pkl".format('rawToDo'), 'w'); pickle.dump(todo,f); f.close()
    
    lt2, noraw2, todotracking =noRaw(no_tracking)
    lt3, noraw3, todofeat = noRaw(no_trajfeat)
    if lt!=lt2 or not np.all(np.array(noraw)==np.array(noraw2)): raise
    print 'LTValid', lt
    print 'no raw', len(noraw)
    print 'raw to do', len(todoraw), 'and raw done', counthdf5
    print 'tracking to do', len(todotracking), 'and tracking done', counttracking
    todotracking=[(line[0], line[1]) for line in todotracking]
    todofeat = [(line[0], line[1]) for line in todofeat]
    print "traj features to do", len(todofeat), 'and traj features done', len(countfeat)
    if name is not None:
        print 'saving features traj exp list in ', "/cbio/donnees/aschoenauer/workspace2/Tracking/{}.pkl".format(name)
        f=open("/cbio/donnees/aschoenauer/workspace2/Tracking/{}.pkl".format(name), 'w'); pickle.dump(countfeat,f); f.close()
    if not featlistonly: return todoraw, todotracking, todofeat
    else:return countfeat


def fromShareToCBIO(expdict,wellFolder=False, dst='/cbio/donnees/aschoenauer/cluster',\
                    folderIn="/share/data20T/mitocheck/compressed_data"):
    '''
    Copying raw data of a dictionary of experiments from /share to /cbio/donnees/aschoenauer
    
    Input:
    - expdict: dictionary of experiments of the type {gene:experiment list}
    - wellFolder: bool, if there should be a folder per plate + inside, a folder per well in the destination folder
    - dst: destination folder
    - folderIn: input folder
    
    '''
    if not wellFolder:
        for gene in expdict:
            for pl, w in expdict[gene]:
                ww=filter(lambda x: w[2:5] == x[:3],os.listdir(os.path.join(folderIn, pl)))
                try:
                    shutil.copytree(os.path.join(folderIn, pl, ww[0]), os.path.join(dst, gene, pl+'_'+ww[0][:3]))
                except:
                    pass
            
    else:
        for gene in expdict:
            for pl, w in expdict[gene]:
                ww=filter(lambda x: w[2:5] == x[:3],os.listdir(os.path.join(folderIn, pl)))
                try:
                    shutil.copytree(os.path.join(folderIn, pl, ww[0]), os.path.join(dst, gene, pl, ww[0][:3]))
                except:
                    pass

def appendingControl(plateL):
    todo=[]
    if type(plateL[0])==tuple:
        plateL = Counter(np.array(plateL)[:,0]).keys()

    for el in plateL:
        id_ = int(el.split('_')[0][2:])
        if id_<50:
            for si in typeD:
                for val in typeD[si]:
                    todo.append((el, '{:>05}_01'.format(int(val))))
        else:
            for si in typeD2:
                for val in typeD2[si]:
                    todo.append((el, '{:>05}_01'.format(int(val))))

    return todo


def gettingSiRNA(labels, who, ctrlStatus, length, genes, qc):
    '''
    This function gives the siRNA list of a set of Mitocheck experiments, also removing experiments
    that did not go through the quality control
    '''
    sirna=[]; unfound=[]
    yqualDict=expSi(qc)
    for k in range(len(who)):
        try:
            sirna.append(yqualDict[who[k][0][:9]+'--'+who[k][1][2:5]])
        except KeyError:
            unfound.append( k)
            sirna.append('unfound')
            print who[k][0][:9]+'--'+who[k][1][2:5]
    todel=[]
    for indice in unfound:
        todel.extend(range(np.sum(length[:indice]), np.sum(length[:indice+1])))
        
    todel.sort(); unfound.sort(reverse=True)
    nlabels = np.delete(labels, todel, 0)
    for indice in unfound:
        del genes[indice]
        del sirna[indice]
        del length[indice]
        del ctrlStatus[indice]
        del who[indice]
    if np.sum(length)!=len(nlabels) or len(length)!=len(who) or len(length)!=len(ctrlStatus) or len(length)!=len(genes) or len(length)!=len(sirna) \
        or 'unfound' in sirna:
        raise
    else: return nlabels, who, ctrlStatus, length, genes, sirna

def EntrezToExp(corresEntrezEnsemble='NatureEntrezEnsemble.txt', mitocheck='mitocheck_siRNAs_target_genes_Ens72.txt', qc='qc_export.txt', sub_EntrezList=None):
    EE=txtToList(corresEntrezEnsemble)
    EElist=[]
    for k in range(len(EE)):
        EElist.append((EE[k][0], EE[k][1]))
    EElist=np.array(EElist) #entrez ID of the Nature genes to ensembl
    
    entrrezs=defaultdict(list)
    if sub_EntrezList == None:
        for k in range(len(EElist[:,0])):
            entrrezs[EElist[k][0]].append(EElist[k][1])
    else:
        for k in range(len(EElist[:,0])):
            if EElist[k][0] in sub_EntrezList:
                entrrezs[EElist[k][0]].append(EElist[k][1])

    mitocheck_siRNAs=txtToList(mitocheck)    
    siEnsembl=[]
    for k in range(len(mitocheck_siRNAs)):
        siEnsembl.append((mitocheck_siRNAs[k][1], mitocheck_siRNAs[k][2]))
        
    siEns=defaultdict(list)
    for k in range(len(siEnsembl)):
        siEns[siEnsembl[k][1]].append(siEnsembl[k][0]) #ensembl ID of the Mitocheck genes to siRNA
    siRNAresult = []
    yeSiExp=expSi(qc, 0)
    expDict=defaultdict(list)
    for entrez in entrrezs:
        for ens in entrrezs[entrez]:
            for sirna in siEns[ens]:
                siRNAresult.append(sirna)
                expDict[entrez].extend(yeSiExp[sirna])

    missingEntrez=filter(lambda x: x not in expDict.keys(), entrrezs.keys())
    return expDict, missingEntrez, siRNAresult

def expSi(qc, sens=1):
    qc=txtToList(qc)
    yes=qc[np.where(qc[:,-2]=='True')]
    
    if sens==0:
        yeSiExp=defaultdict(list)
    #siRNA to plate,well for experiments that have passed the quality control    
        for k in range(yes.shape[0]):
            yeSiExp[yes[k,1]].append(yes[k,0])
    else:
        yeSiExp={}
        for k in range(yes.shape[0]):
            yeSiExp[yes[k,0]]=yes[k,1]
    return yeSiExp

def siEntrez(mitocheck):
    mitocheck=txtToList(mitocheck)
    result = {}
    for k in range(len(mitocheck)):
        if len(mitocheck[k])==4:
            result[mitocheck[k][1]]=mitocheck[k][-1]

    return result

class ArffReader(object):
    '''
    modified (from Cecog + me) ARFF reader from Pradeep Kishore Gowda
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/440533
    '''

    DCT_DATA_TYPES = {'numeric': float,
                      'string' : str,
                     }

    def __init__(self, strFilename):
        self.oFile = file(strFilename, "r")
        self.strRelationName = ""
        self.dctClassNames = {}
        self.dctClassLabels = {}
        self.lstFeatureNames = []
        self.dctDataTypes = {}
        self.dctFeatureData = {}
        self.dctHexColors = {}
        self.hasZeroInsert = True
        self._read()
        self._close()

    def _convert_string(self, string):
        return string.replace('\'','').replace('"','')

    def _read(self):
        in_data = False

        lstClassNames = []
        lstClassLabels = []
        lstHexColors = []

        for line in self.oFile.readlines():
            line = line.rstrip()

            if len(line) > 0:

                if line[0] in ['@', '%']:
                    lstToken = line.split(' ')
                    first = lstToken[0].lower()

                    if first == '@relation':
                        self.strRelationName = lstToken[1]

                    elif first == '@attribute':
                        if lstToken[2][0] == '{':
                            i1 = line.find('{')
                            i2 = line.find('}')
                            assert i1 != i2 != -1
                            lstClassNames = [self._convert_string(x).strip(" ")
                                             for x in line[i1+1:i2].split(',')]
                            setClassNames = set(lstClassNames)
                        else:
                            strFeatureName = lstToken[1]
                            strDataType = lstToken[2].lower()
                            self.lstFeatureNames.append(strFeatureName)
                            self.dctDataTypes[strFeatureName] = \
                                self.DCT_DATA_TYPES[strDataType]

                    elif first == '%class-labels':
                        if lstToken[1][0] == '{':
                            i1 = line.find('{')
                            i2 = line.find('}')
                            assert i1 != i2 != -1
                            lstClassLabels = [int(self._convert_string(x))
                                              for x in
                                              line[i1+1:i2].split(',')]

                    elif first == '%class-colors':
                        if lstToken[1][0] == '{':
                            i1 = line.find('{')
                            i2 = line.find('}')
                            assert i1 != i2 != -1
                            lstHexColors = [self._convert_string(x).upper()
                                            for x in
                                            line[i1+1:i2].split(',')]

                    elif first == '%has-zero-inserted-in-feature-vector':
                        self.hasZeroInsert = int(lstToken[1]) != 0

                    elif first == '@data':
                        in_data = True

                elif in_data:
                    lstItems = line.split(',')
                    strClassName = self._convert_string(lstItems[-1])
                    lstItems = lstItems[:-1]
                    #print strClassName, lstClassNames
                    assert strClassName in setClassNames
                    assert len(lstItems) == len(self.lstFeatureNames)
                    lstFeatureData = []
                    for strFeatureName, item in zip(self.lstFeatureNames,
                                                    lstItems):
                        if self.dctDataTypes[strFeatureName] == str:
                            value = self._convert_string(item)
                        else:
                            value = self.dctDataTypes[strFeatureName](item)
                        lstFeatureData.append(value)
                    if strClassName not in self.dctFeatureData:
                        self.dctFeatureData[strClassName] = []
                    self.dctFeatureData[strClassName].append(lstFeatureData)

        for strClassName in self.dctFeatureData:
            self.dctFeatureData[strClassName] = np.array(self.dctFeatureData[strClassName], np.float)

        for iClassLabel, strClassName in zip(lstClassLabels, lstClassNames):
            self.dctClassLabels[strClassName] = iClassLabel
            self.dctClassNames[iClassLabel] = strClassName

        for hexColor, strClassName in zip(lstHexColors, lstClassNames):
            self.dctHexColors[strClassName] = hexColor


    def _close(self):
        self.oFile.close()


    
if __name__ == '__main__':
    checkingAllHDF5()
