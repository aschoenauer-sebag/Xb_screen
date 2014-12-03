import os, pdb
import cPickle as pickle
import numpy as np
import vigra.impex as vi
from collections import defaultdict

def intensity_qc(input_folder, output_folder, visual_qc={}, flou_qc={}):
    f=open(os.path.join(input_folder, 'intensity_qc.pkl'))
    d=pickle.load(f)
    f.close()
    
    f=open(os.path.join(output_folder, 'qc_xbscreen.txt'), 'w')
    f.write('{}\t{}\t{}\t{}\t{}\n'.format('Plate', 'Well', 'WellQC', 'ImageQC', 'ManualQC'))
    
    line_type='{}\t{}\t{}\t{}\t{}\n'
    for plate in d:
        for well in d[plate]:
            manual_qc=(plate not in visual_qc or well not in visual_qc[plate])
            flou=(plate not in flou_qc or well not in flou_qc[plate])
            image_qc=list(np.array(d[plate][well])+1)
            f.write(line_type.format(plate, well, flou, image_qc, manual_qc ))
    f.close()
    return
