import numpy as np
import os
from tracking.importPack import imp, EVENTS

#returns two dictionaries : one with the cells that have NaN entries in their features, with the count of the number of NaN, and one with
#the features that have NaN entries

def whichAreNan(X):
    nanB = np.isnan(X)
    nan = np.where(nanB==True)
    cellsWithNan = nan[0]
    featuresWithNan = nan[1]
    c={}
    f={}
    c = nan#np.where(nan[1] in [77,78])
#    for i in cellsWithNan:
#        if i not in c:
#            c[i]=1
#        else:
#            c[i]+=1
    for i in featuresWithNan:
        if i not in f:
            f[i]=1
        else:
            f[i]+=1
            
    return c,f

def returnBadFeatures():
    print os.getcwd()
    tabF= None
#    small = [128, 129, 130, 131, 132, 133, 134, 169, 170, 173, 174, 177, 178, 181, 182, 189, 197, 198, 201, 202, 109, 190, 105, 106, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 123, 124, 125, 126, 127]
#small pr <500
    small = [128, 129, 130, 131, 132, 133, 134, 169, 170, 173, 174, 177, 178, 181, 182, 189, 197, 198, 201, 202, 109, 16, 17, 190, 105, 106, 107, 108, 18, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 123, 124, 125, 126, 127]
    smallR = {}
    listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
    for plate in listD:
        print plate
        listW = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate)
        for well in listW:
            well=well[:-5]
            print well
            filename = '/media/lalil0u/New/workspace2/Tracking/data/raw/'+plate+"/"+well+".hdf5"
#        filenameT = '/media/lalil0u/New/workspace2/Tracking/data/trainingset/'+plate+"_P"+well+".hdf5"

            print "features loading : plate = "+plate+",well = "+well
            result = imp.importFeaturesOnly(filename, plate, well)
#            i=0
#            for c in result:
#                for f in small:
#                    if np.isnan(c[f]):
#                        if f not in smallR:
#                            smallR[f] = []
#                        smallR[f].append((plate, well, i))
#                i+=1
                                                          
            if tabF==None:
                tabF = result
            else:
                tabF=np.vstack((tabF, result))
            
    featuresNames = imp.importFeaturesNames(filename)
#    print featuresNames, type(featuresNames), featuresNames.shape
    c, f = whichAreNan(tabF)
    print "nb de features pblmatiques :", len(f.keys())
#    print "nom des small features", featuresNames[small]
#    print "lignes des cellules avec les small NaN :" 
#    for feat in smallR:
#        print "\n", featuresNames[feat][0], "\n", smallR[feat][0], smallR[feat][1] 
#    for feat in f.keys():
#        print ("{}{:> 5d} ".format(featuresNames[feat][0], f[feat]))  
    return f, featuresNames

def clean(X,f):
    return np.delete(X, f, 1)

def minMax(tab):
    result = {"min":[], "max":[]}
    for feat in range(tab.shape[1]): 
        result["min"].append(min(np.transpose(tab)[feat]))
        result["max"].append(max(np.transpose(tab)[feat]))
        
    return result
    

def normalize(Mat1, Mat2=None):
    if Mat2==None:
        Mat2 = Mat1
    average = np.mean(Mat2, 0)      
    std=np.std(Mat2, 0)
    return (Mat1 - average)/std
#la forme attendue est que les features a normaliser soient en colonnes et les lignes les objets

def f(sol):
    r = [0 for x in range(len(EVENTS))]
    i=0; v = sol.truthVec()
    for t in v:
        r[i]=sum(t)
        i+=1
#    total = sum(r)
#    r=list(np.array(r)/float(total)*100)
#    i=0
#    for k in r:
#        reste = k%5
#        plus = 5 if 5-reste<reste else 0
#        r[i]=k-reste + plus
#        i+=1
    return r
    