import os, csv, getpass, pandas
import numpy as np
import cPickle as pickle
from sklearn.manifold import TSNE
import matplotlib.pyplot as p

primary_channel_name = 'primary__primary3'
pathClassification = "/sample/0/plate/{}/experiment/{}/position/1/feature/%s/object_classification/prediction"%primary_channel_name

if getpass.getuser()!='lalil0u':
    import vigra.impex as vi
    
def plot_results(norm=True):
    #Opening rank sum statistics
    f=open('/mnt/projects/drug_screen/results/ranksum1.pkl')
    mitocheck, ds = pickle.load(f)
    f.close()
    
    f=open('/mnt/projects/drug_screen/results/DS_hits_well_drugs.pkl')
    wells, drugs = pickle.load(f); f.close()
    
    ds_ctrl = ds[ds.Well=='CTRL'].iloc[:,2:]
    mitocheck_ctrl = mitocheck[mitocheck.Well=='CTRL'].iloc[:,2:]
    if norm:
        ds_m, mito_m = np.mean(ds_ctrl, 0), np.mean(mitocheck_ctrl, 0)
        ds_std, mito_std = np.std(ds_ctrl, 0), np.std(mitocheck_ctrl, 0)
        
        nds = (ds.iloc[:,2:]-ds_m)/ds_std
        nmito_c = (mitocheck_ctrl - mito_m)/mito_std
    else:
        nds = ds.iloc[:,2:]
        nmito_c = mitocheck_ctrl
    
    tsne = TSNE(n_components = 2)
    a = tsne.fit_transform(np.vstack((nds, nmito_c)))
    
    colors = []
    for index,el in ds.iterrows():
        if el.Well == "CTRL":
            colors.append('g')
        elif int(el.Well) in wells:
            colors.append('r')
        else:
            colors.append('black')
            
    colors.extend(['c' for k in range(mitocheck_ctrl.shape[0])])
    p.scatter(a[:,0], a[:,1], color=colors, alpha=0.7)
    p.scatter(0,0, color='g', label='DS ctrl')
    p.scatter(0,0, color='r', label='DS hit')
    p.scatter(0,0, color='black', label='DS non-hit')
    p.scatter(0,0, color='c', label='Mitocheck ctrl')
    p.legend()
    p.title('Rank-sum statistics, t-SNE projection')
    p.show()
    
    
    #now looking at raw phenotypic profiles
    f=open('/mnt/projects/drug_screen/results/all_joint_aggregated_values_oldmito.pkl')
    ds, mitocheck, types, phenotypes=pickle.load(f)
    f.close()
    
    
    
    return
    

def check_h5_files_exist(rawDataFolder = '/share/data20T/mitocheck/compressed_data', 
                         h5Folder = '/share/data40T/aschoenauer/drug_screen/results_August_2016/mito_joint_classifier'):
    '''
    Goal: see how many H5 files are missing with respect to the raw data
'''
    plates = filter(lambda x: 'Valid' not in x, os.listdir(rawDataFolder))
    missing_h5=[]
    past_len=0
    
    for plate in plates:
        well_h5files = ['00{}_01'.format(el.split('--')[0]) for el in os.listdir(os.path.join(rawDataFolder, plate)) if el.split('.')[-1]!='xml']
        try:
            existing_h5files=set(os.listdir(os.path.join(h5Folder, plate, 'hdf5')))
        except:
            missing_h5.extend([(plate, el) for el in well_h5files])
        else:
            missing_h5.extend([(plate, el) for el in well_h5files if '{}.ch5'.format(el) not in existing_h5files])
        if len(missing_h5)>past_len:
            print plate, ' missing ', len(missing_h5)-past_len
            
        past_len = len(missing_h5)
        
    return missing_h5 

def kidney_looker(h5Folder = '/share/data40T/aschoenauer/drug_screen/results_August_2016/mito_joint_classifier'):
    '''
   Goal: identify if there exists certain siRNAs which are characterized by a high number of "Kidney" nuclei 
   Kidney corresponds to label 17
'''
    
    f=open('../data/mapping_2014/qc_export.txt', 'r')
    reader = csv.reader(f, delimiter='\t'); reader.next()
    whatList = []
    plateList = np.array(os.listdir(h5Folder))
    i=0
    for el in reader:
        if 'Valid' not in el[0] and int(el[0][2:6])<600 and el[-1]=="ok":
            plate = el[0].split('--')[0]
            well = el[0].split('--')[1]
            
            truePlate = plateList[np.where([plate in p for p in plateList])[0]][0]
            
            if 'hdf5' in os.listdir(os.path.join(h5Folder, truePlate)) and '00{}_01.ch5'.format(well) in os.listdir(os.path.join(h5Folder, truePlate, 'hdf5')):
                pathClassif = pathClassification.format(truePlate, '00{}'.format(well))
                try:
                    tabClassification = np.array(vi.readHDF5(os.path.join(h5Folder, truePlate, 'hdf5', '00{}_01.ch5'.format(well)), 
                                                         pathClassif), dtype=int)
                except IOError:
                    print "Error ", truePlate, well
                    continue
                
                kidneys = np.bincount(tabClassification, minlength=18)[-1]/float(tabClassification.shape[0])
                
                
                if el[2]=='scrambled':
                    whatList.append((truePlate, well, 'CT', kidneys))
                else:
                    whatList.append((truePlate, well, el[1], kidneys))
                i+=1
                if i%1000==0:
                    print "Done ", i
                    
    return pandas.DataFrame.from_records(whatList, columns = ('Plate', 'Well', 'siRNA', 'KidneyPerc'))
                    
            