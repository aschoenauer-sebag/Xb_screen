import os, pdb
import cPickle as pickle
import numpy as np
import vigra.impex as vi
from collections import defaultdict

from plateSetting import fromXBToWells

def intensity_qc(input_folder, output_folder, plateL=['201114', '271114', '121214']):
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
    pdb.set_trace()
    line_type='{}\t{}\t{}\t{}\t{}\n'
    for plate in plateL:
        for well in d_intensity[plate]:
            visual_qc=(plate not in d_manual or well not in d_manual[plate])
            flou=(plate not in flou_qc or well not in flou_qc[plate])
            image_qc=list(np.array(d_intensity[plate][well])+1)
            f.write(line_type.format(plate, well, flou, image_qc, visual_qc ))
    f.close()
    return

def computingToRedo(threshold_flou=0.4, threshold_init_cell_count=20, threshold_control_count=3,
                    input_folder='../data', 
                    compoundL=['TCDD', 'MeHg', 'BPA', 'PCB', 'Endo','DMSO', 'Nonane', 'Rien', 'TGF'],
                     plateL=['201114', '271114', '121214'],
                     hdf5Folder = "/media/lalil0u/New/projects/Xb_screen/plates__all_features_2bis",
                     savingFolder = "/media/lalil0u/New/projects/Xb_screen/dry_lab_results"):
    '''
    What should the goal of this function be?
    I have a list of xenobiotics and doses and I want to know if I have three experiments on a different day for each
    
    ALSO
    the number of experiments discarded by this qc rule
    
    ALSO
    normally there are five controls/solvent/plate. There should not be less than three otherwise the plate is to be discarded (oh god)
    '''
    failed=defaultdict(dict)
    number_failed=0; total_number=0
    flou_qc_dict=defaultdict(list)
    _,well_groups=fromXBToWells(compoundL)
    print "Loading manual quality control results"
    f=open(os.path.join(input_folder, 'xb_manual_qc.pkl'), 'r')
    d_manual=pickle.load(f)
    f.close()
    
    result=defaultdict(dict)
    for plate in plateL:
        f=open(os.path.join(savingFolder, 'processedDictResult_P{}.pkl'.format(plate)))
        resCour=pickle.load(f)
        f.close()
        result[plate]=resCour
    
    for compound in compoundL:
        print "################### Loading wells for xenobiotic ", compound
        #by default this function returns all wells with this xenobiotic starting on the 20th of Nov
        curr_well_groups=well_groups[compound]
        failed[compound]=defaultdict(list)
        for dose in curr_well_groups:
            print "----DOSE ", dose
            
            for plate in curr_well_groups[dose]:
                for well in curr_well_groups[dose][plate]:
                    total_number+=1
                    #i.checking if the well is in the manual qc failed list
                    if plate in d_manual and well in d_manual[plate]:
                        print 'Manual QC failed', plate, well
                        failed[compound][dose].append((plate, well))
                        number_failed+=1
                    else:
                        #ii. checking if the initial number of objects is above the threshold
                        avg_init_cell_count=np.mean(np.array(result[plate][well]['cell_count'])[:10])
                        if avg_init_cell_count<threshold_init_cell_count:
                            print 'Cell count failed', plate, well
                            failed[compound][dose].append((plate, well))
                            number_failed+=1
                            
                            flou_qc_dict[plate].append(well)
                            
                        #iii. checking if the percentage of out of focus objects in the last frames is ok
                        else:
                            end_flou_perc=np.mean(np.array(result[plate][well]['Flou_ch1'])[190:200])
                            if end_flou_perc>threshold_flou:
                                print 'Out of focus count failed', plate, well
                                failed[compound][dose].append((plate, well))
                                number_failed+=1
                                
                                flou_qc_dict[plate].append(well)
    print 'Saving list of failed wells in ', os.path.join(input_folder, 'xb_focus_qc.pkl')
    f=open(os.path.join(input_folder, 'xb_focus_qc.pkl'), 'w')
    pickle.dump(flou_qc_dict, f); f.close()
    print "Percentage of experiments failed ", number_failed/float(total_number)
    return failed, number_failed, total_number
            















