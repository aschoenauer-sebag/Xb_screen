import os, pdb, shutil
from collections import defaultdict
from operator import itemgetter
import cPickle as pickle
import numpy as np

from fileManagement import txtToList, strToTuple
#ctrl wells for labteks, id<50
typeD={}
typeD["scrambled"] = ["015", "026", "063", "074", "304", "315", "352"]
typeD["bcop"] = ["363"]
typeD["eg5"] = ["290", "301", "338", "349"]
typeD["incenp"] = ["004", "049", "052"]
typeD["empty"] = []
typeD["marker"] = ["001"]

typeD2={}
typeD2["scrambled"] = ["015", "026", "063", "074", "304", "315", "352", "363"]
typeD2["bcop"] = ["333", "336", "381", "384"]
typeD2["eg5"] = ["290", "301", "338", "349"]
typeD2["incenp"] = ["004", "049", "052"]
typeD2["empty"] = ["311", "322", "359", "370"]
typeD2["marker"] = ["001"]

def appendingControl(plateL):
    todo=[]

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

def EntrezToExp(corresEntrezEnsemble='NatureEntrezEnsemble.txt', mitocheck='mitocheck_siRNAs_target_genes_Ens72.txt', qc='qc_export.txt'):
    EE=txtToList(corresEntrezEnsemble)
    EElist=[]
    for k in range(len(EE)):
        EElist.append((EE[k][0], EE[k][1]))
    EElist=np.array(EElist) #entrez ID of the Nature genes to ensembl
    
    entrrezs=defaultdict(list)
    for k in range(len(EElist[:,0])):
        entrrezs[EElist[k][0]].append(EElist[k][1])

    mitocheck_siRNAs=txtToList(mitocheck)    
    siEnsembl=[]
    for k in range(len(mitocheck_siRNAs)):
        siEnsembl.append((mitocheck_siRNAs[k][1], mitocheck_siRNAs[k][2]))
        
    siEns=defaultdict(list)
    for k in range(len(siEnsembl)):
        siEns[siEnsembl[k][1]].append(siEnsembl[k][0]) #ensembl ID of the Mitocheck genes to siRNA
    yeSiExp=expSi(qc, 0)
    expDict=defaultdict(list)
    for entrez in entrrezs:
        for ens in entrrezs[entrez]:
            for sirna in siEns[ens]:
                for exp in yeSiExp[sirna]:
                    expDict[entrez].append(exp)
    missingEntrez=filter(lambda x: x not in expDict.keys(), entrrezs.keys())
    return expDict, missingEntrez

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

