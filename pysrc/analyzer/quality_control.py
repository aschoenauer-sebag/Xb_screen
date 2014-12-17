import os, pdb
import cPickle as pickle
import numpy as np
import vigra.impex as vi
from collections import defaultdict

from plateSetting import fromXBToWells

def intensity_qc(input_folder, output_folder, flou_qc={}):
    print "Loading automatic quality control results"
    f=open(os.path.join(input_folder, 'xb_intensity_qc.pkl'), 'r')
    d_intensity=pickle.load(f)
    f.close()
    print "Loading manual quality control results"
    f=open(os.path.join(input_folder, 'xb_manual_qc.pkl'), 'r')
    d_manual=pickle.load(f)
    f.close()
    
    f=open(os.path.join(output_folder, 'qc_xbscreen.txt'), 'w')
    f.write('{}\t{}\t{}\t{}\t{}\n'.format('Plate', 'Well', 'WellQC', 'ImageQC', 'ManualQC'))
    
    line_type='{}\t{}\t{}\t{}\t{}\n'
    for plate in d_intensity:
        for well in d_intensity[plate]:
            visual_qc=(plate not in d_manual or well not in d_manual[plate])
            flou=(plate not in flou_qc or well not in flou_qc[plate])
            image_qc=list(np.array(d_intensity[plate][well])+1)
            f.write(line_type.format(plate, well, flou, image_qc, visual_qc ))
    f.close()
    return

def computingToRedo(threshold_flou=0.5, threshold_init_cell_count=30, threshold_control_count=3,
                    input_folder='../data', 
                    xbL=['TCDD', 'MeHg', 'BPA', 'PCB', 'Endo'],
                    others=['DMSO', 'Nonane', 'Rien', 'TGF'],
                     plateL=['201114', '271114', '121214'],
                     hdf5Folder = "/media/lalil0u/New/projects/Xb_screen/plates__all_features_2bis",
                     savingFolder = "/media/lalil0u/New/projects/Xb_screen/dry_lab_results/track_predictions__settings2"):
    '''
    What should the goal of this function be?
    I have a list of xenobiotics and doses and I want to know if I have three experiments on a different day for each
    
    ALSO
    the number of experiments discarded by this qc rule
    
    ALSO
    normally there are five controls/solvent/plate. There should not be less than three otherwise the plate is to be discarded (oh god)
    '''
    
    print "Loading manual quality control results"
    f=open(os.path.join(input_folder, 'xb_manual_qc.pkl'), 'r')
    d_manual=pickle.load(f)
    f.close()
    for xb in xbL:
        for plate in plateL:
            well_groups=fromXBToWells(xb, plate=plate)
