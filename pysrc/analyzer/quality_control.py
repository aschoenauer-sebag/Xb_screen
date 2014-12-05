import os, pdb
import cPickle as pickle
import numpy as np
import vigra.impex as vi
from collections import defaultdict

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
