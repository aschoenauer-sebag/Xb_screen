import os, pdb
import cPickle as pickle
import numpy as np
import vigra.impex as vi
import matplotlib.pyplot as p
from collections import defaultdict, Counter
from scipy.stats import scoreatpercentile

from analyzer import compoundL, plates, total_experiment_number
from plateSetting import fromXBToWells
from util.plots import couleurs
from scipy.stats.stats import spearmanr

def normalizationCorrelationMeasures(arr_plate, arr_negctrl, who_plate, who_negctrl):
    
    who_plate2=['{}--{}'.format(*el) for el in who_plate]
    who_negctrl2=['{}--{}'.format(*el) for el in who_negctrl]
    
    comm=filter(lambda x: x in who_plate2, who_negctrl2)
    
    return arr_plate[np.where(np.array([el in comm for el in who_plate2]))], arr_negctrl[np.where(np.array([el in comm for el in who_negctrl2]))]

def conditionCorrelationMeasures(arr, who,conditions, condition_list):
    result=np.zeros(shape=(len(plates), len(plates)))
    if type(arr)==list:
        arr=np.array(arr)
    if type(conditions)==list:
        conditions=np.array(conditions)
    if type(who)==list:
        who=np.array(list)
    
    for i in range(len(plates)):
        result[i,i]=1
#         for cond in condition_list:
#             print conditions[np.where((conditions==cond)&(who[:,0]==plates[i]))],
#             print who[np.where((conditions==cond)&(who[:,0]==plates[i]))]
        list_val1=[np.median(arr[np.where((conditions==cond)&(who[:,0]==plates[i]))]) for cond in condition_list]
        #print list_val1
        for j in range(i+1, len(plates)):
            list_val2=[np.median(arr[np.where((conditions==cond)&(who[:,0]==plates[j]))]) for cond in condition_list]
            aa=np.vstack((list_val1, list_val2))
            aa=aa[:,np.where((np.isfinite(aa[1])&(np.isfinite(aa[0]))))]
            result[i,j]=spearmanr(aa[0,0], aa[1,0])[0]
            result[j,i]=result[i,j]
    
    return result

def replicateCorrelationPlots(arr, who, conditions, ax=None, refPlate=None):
    '''
    This is to look at replicability from a well perspective (and not condition perspective).
    If the 
    Arguments:
    - conditions: list of conditions ordered as the data in arr
    - who: list or array of (plate, well) info
    '''
    if refPlate==None:
        if type(who)==list:
            who=np.array(who)
        if type(conditions)==list:
            conditions=np.array(conditions)

        f,axes=p.subplots(3,4,sharex=True, sharey=True)
        for i,plate in enumerate(['201114', '121214', '271214']):
            replicateCorrelationPlots(arr, who, conditions, ax=axes[i], refPlate=plate)
        p.show()
        return
    else:
        if ax==None:
            raise ValueError
        
        for i,currPlate in enumerate(filter(lambda x: x!=refPlate, plates)):
            currConditions = conditions[np.where(who[:,0]==currPlate)]
            currArr = arr[np.where((who[:,0]==currPlate))]
            for k,cond in enumerate(currConditions):
                refCond = np.median(arr[np.where((conditions==cond)&(who[:,0]==refPlate))])
                ax[i].scatter(refCond, currArr[k], color=couleurs[plates.index(currPlate)])
                ax[i].plot(np.linspace(np.min(currArr), np.max(currArr), 10), np.linspace(np.min(currArr), np.max(currArr), 10))
                ax[i].set_xlabel(refPlate); ax[i].set_ylabel(currPlate); ax[i].grid(True)
                
        return
            
        

def conditionReplicabilityBoxplots(arr, condition_list, conditions, key_list=None):
    '''
    This function looks at replicability from a condition perspective.
    It plots boxplots of the data that is in arr for each condition given in condition_list. Conditions should contain
    once all possible conditions.
    '''
    res=defaultdict(list)
    data=[]
    
    if key_list==None:
        for condition in conditions:
            data.append(np.median(arr,0)[np.where(condition_list==condition)]/np.median(np.median(arr,0)[np.where(condition_list==condition)]))
            data.append(np.max(arr,0)[np.where(condition_list==condition)]/np.median(np.max(arr,0)[np.where(condition_list==condition)]))
        for el in data[::2]:
            res['Median'].append(scoreatpercentile(a=el, per=75,axis=0)-scoreatpercentile(a=el, per=25,axis=0))
        for el in data[1::2]:
            res['Max'].append(scoreatpercentile(a=el, per=75,axis=0)-scoreatpercentile(a=el, per=25,axis=0))
        title="Comparison of replicability for trajectory statistic for all conditions, between med(iterations) and max(iterations)"
    else:
        for condition in conditions:
            for key in key_list:
                el = arr[key][np.where(condition_list==condition)]/np.median(arr[key][np.where(condition_list==condition)])
                data.append(el)
                res[key].append(scoreatpercentile(a=el, per=75,axis=0)-scoreatpercentile(a=el, per=25,axis=0))
        title="Replicability for {} distance for all conditions".format(key)

    f=p.figure()
    ax=f.add_subplot(121)
    ax.boxplot(data)
    if key is None:
        xtickNames = p.setp(ax, xticklabels=np.repeat(conditions, 2))
    else:
        xtickNames = p.setp(ax, xticklabels=np.repeat(conditions, len(key_list)))
    ax.grid(True)
    p.setp(xtickNames, rotation=90, fontsize=8)
    p.title(title)
    ax=f.add_subplot(122)
    for el in res:
        ax.hist(res[el], label=el, bins=50, range=(0,1), alpha=0.5)
    ax.legend(); ax.grid(True)
    p.title('IQR distributions')
    
    p.show()
            


def zFactor(data, ctrlStatus):
    '''
    Function to compute z-factor to evaluate screen quality
    '''
    
    if type(data)==list or len(data.shape)==1 or data.shape[1]==1:
        ctrlStatus=np.array(ctrlStatus)
        ctrl=np.array(data)[np.where(ctrlStatus==0)]
        experiments = np.array(data)[np.where(ctrlStatus==1)]
        
        return 1-3*(np.std(ctrl)+np.std(experiments))/np.absolute(np.mean(ctrl)-np.mean(experiments))
    else:
        result = np.zeros(shape=(data.shape[1],), dtype=float)
        for k in range(data.shape[1]):
            result[k]=zFactor(data[:,k], ctrlStatus)
        return result

def intensity_qc(input_folder, output_folder):
    print 'Doing quality control text file for plates ', plates
    
    print "Loading automatic intensity quality control results"
    f=open(os.path.join(input_folder, 'xb_intensity_qc.pkl'), 'r')
    d_intensity=pickle.load(f)
    f.close()
    print "Loading automatic out of focus and cell count quality control results"
    f=open(os.path.join(input_folder, 'xb_focus_qc.pkl'), 'r')
    flou_qc=pickle.load(f)
    f.close()
   
    print "Loading manual quality control results"
    f=open(os.path.join(input_folder, 'xb_manual_qc.pkl'), 'r')
    d_manual=pickle.load(f)
    f.close()
    
    f=open(os.path.join(output_folder, 'qc_xbscreen.txt'), 'w')
    f.write('{}\t{}\t{}\t{}\t{}\n'.format('Plate', 'Well', 'WellQC', 'ImageQC', 'ManualQC'))
    line_type='{}\t{}\t{}\t{}\t{}\n'
    for plate in plates:
        for well in d_intensity[plate]:
            visual_qc=(plate not in d_manual or well not in d_manual[plate])
            flou=(plate not in flou_qc or well not in flou_qc[plate])
            image_qc=list(np.array(d_intensity[plate][well])+1)
            f.write(line_type.format(plate, well, flou, image_qc, visual_qc ))
    f.close()
    return

def usable_XBSC(compound, dose, plate=None, input_folder='../data', confDir='/media/lalil0u/New/projects/Xb_screen/protocols_etal/plate_setups'):
    _, info_dict=fromXBToWells([compound], dose_filter=dose, plate=plate,confDir=confDir)
    
    print "Loading manual quality control results"
    f=open(os.path.join(input_folder, 'xb_manual_qc.pkl'), 'r')
    d_manual=pickle.load(f)
    f.close()
    print "Loading automatic out of focus and cell count quality control results"
    f=open(os.path.join(input_folder, 'xb_focus_qc.pkl'), 'r')
    flou_qc=pickle.load(f)
    f.close()
    
    expL=[]
    if plate is not None:
        wells=info_dict[compound][dose][plate]
        for well in wells:
            if ((plate not in flou_qc) or (plate in flou_qc and well not in flou_qc[plate])):
                if ((plate not in d_manual) or (plate in d_manual and well not in d_manual[plate])):
                    expL.append((plate, '{:>05}'.format(well)))
    else:
        for plate in info_dict[compound][dose]:
            wells=info_dict[compound][dose][plate]
            for well in wells:
                if ((plate not in flou_qc) or (plate in flou_qc and well not in flou_qc[plate])):
                    if ((plate not in d_manual) or (plate in d_manual and well not in d_manual[plate])):
                        expL.append((plate, '{:>05}'.format(well)))
                
    return expL         
    
    

def computingToRedo(threshold_flou=0.37, threshold_init_cell_count=23, threshold_control_count=3,
                    input_folder='../data', 
                     hdf5Folder = "/media/lalil0u/New/projects/Xb_screen/plates__all_features_2bis",
                     savingFolder = "/media/lalil0u/New/projects/Xb_screen/dry_lab_results",
                     verbose=False):
    '''
    What should the goal of this function be?
    I have a list of xenobiotics and doses and I want to know if I have three experiments on a different day for each
    
    ALSO
    the number of experiments discarded by this qc rule
    
    ALSO
    normally there are five controls/solvent/plate. There should not be less than three otherwise the plate is to be discarded (oh god)
    '''
    print 'Counting experiments starting with plates ', plates
    
    failed=defaultdict(dict); passed=defaultdict(dict)
    number_failed=0; total_number=0
    flou_qc_dict=defaultdict(list)
    _,well_groups=fromXBToWells(compoundL)
    
    who=[]; qc_indicators = np.empty(shape=(total_experiment_number,3), dtype=float); qc_indicators.fill(-1)
    
    print "Loading manual quality control results"
    f=open(os.path.join(input_folder, 'xb_manual_qc.pkl'), 'r')
    d_manual=pickle.load(f)
    f.close()
    result=defaultdict(dict)
    for plate in plates:
        f=open(os.path.join(savingFolder, 'processedDictResult_P{}.pkl'.format(plate)))
        resCour=pickle.load(f)
        f.close()
        result[plate]=resCour
    for compound in compoundL:
        print "################################# Loading wells for xenobiotic ", compound
        #by default this function returns all wells with this xenobiotic starting on the 20th of Nov
        curr_well_groups=well_groups[compound]
        failed[compound]=defaultdict(list)
        passed[compound]=defaultdict(list)
        for dose in curr_well_groups:
            print "----DOSE ", dose
            total_bio_replicates=0; total_num_wells=0
            for plate in curr_well_groups[dose]:
                usable_plate=False

                for well in curr_well_groups[dose][plate]:
                    total_number+=1
                    print plate, well; who.append((plate, well, compound, dose))
                    #i.checking if the well is in the manual qc failed list
                    if plate in d_manual and well in d_manual[plate]:
                        print 'Manual QC failed', plate, well
                        failed[compound][dose].append((plate, well))
                        number_failed+=1
                    else:
                        #ii. checking if the initial number of objects is above the threshold
                        try:
                            avg_init_cell_count=np.mean(np.array(result[plate][well]['cell_count'])[:10])
                        except ValueError:
                            print "                                                                        MISSING INFO "
                            continue
                        if verbose:
                            print "Avg init cell count ", avg_init_cell_count
                        qc_indicators[total_number-1, 0]=avg_init_cell_count
                        
                        avg_init_cell_perc=np.mean((np.array(result[plate][well]['cell_count'], dtype=float)/np.array(result[plate][well]['object_count']))[:10])
                        if verbose:
                            print "Avg init cell perc ", avg_init_cell_perc
                        qc_indicators[total_number-1, 1]=avg_init_cell_perc
                        
                        flou_arr=np.array(result[plate][well]['Flou_ch1'])
                        end_flou_perc=np.mean(flou_arr[min(190, flou_arr.shape[0]-10):min(200, flou_arr.shape[0])])
                        if verbose:
                            print "Avg end out of focus perc ", end_flou_perc
                        qc_indicators[total_number-1, 2]=end_flou_perc
                        
                        if avg_init_cell_count<threshold_init_cell_count:
                            if verbose:
                                print 'Cell count failed', plate, well
                            failed[compound][dose].append((plate, well))
                            number_failed+=1
                            
                            flou_qc_dict[plate].append(well)
                            
                        #iii. checking if the percentage of out of focus objects in the last frames is ok
                        elif end_flou_perc>threshold_flou:
                            if verbose:
                                print 'Out of focus count failed', plate, well
                            failed[compound][dose].append((plate, well))
                            number_failed+=1
                            
                            flou_qc_dict[plate].append(well)
                        else:
                            passed[compound][dose].append((plate,well))
                            usable_plate=True
                            total_num_wells+=1
            
                if usable_plate:
                    total_bio_replicates+=1
            print '****number of wells ', total_num_wells
            print '*****************************biological rep ', total_bio_replicates
    print 'Saving list of failed wells in ', os.path.join(input_folder, 'xb_focus_qc.pkl')
    f=open(os.path.join(input_folder, 'xb_focus_qc.pkl'), 'w')
    pickle.dump(flou_qc_dict, f); f.close()
    print "Percentage of experiments failed ", number_failed/float(total_number)
    
    return passed, failed, number_failed, np.hstack((np.array(who), qc_indicators))
            















