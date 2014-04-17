import os
import numpy as np
from importPack import imp
from dataPack import treatments
import pylab as p

print os.getcwd()
tabF= None
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
        if tabF==None:
            tabF = result
        else:
            tabF=np.vstack((tabF, result))
            
featuresNames = imp.importFeaturesNames(filename)
#print featuresNames, type(featuresNames), featuresNames.shape,type(tabF), type(tabF[0][0])
c, f = treatments.whichAreNan(tabF)
print "nb de features pblmatiques :", len(f.keys())

small = []
for feat in f.keys():
    print ("{}{:> 5d} ".format(featuresNames[feat][0], f[feat]))
    if f[feat]<500:
        small.append(feat)
print small    
names=featuresNames[f.keys()]
val = f.values()
pos = p.arange(len(f.keys()))+0.5
p.barh(pos, val, align='center')
p.yticks(pos, names)
p.xlabel('Number of objects for which this feature cannot be defined')
#p.title('Training set size :'+str(total)+" objects")
p.grid(True)
p.show()

#    if f[feat]<11:
#        small.append(feat)
        
#return f
#print "features pblmatiques :", f