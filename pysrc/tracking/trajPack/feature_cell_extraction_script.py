import os, pdb, time
from optparse import OptionParser
from collections import Counter
import numpy as np
import cPickle as pickle

from feature_cell_extraction import nb_exp_list
from util.listFileManagement import expSi, appendingControl, strToTuple, countingDone

jobSize = 10
progFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc'
scriptFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/scripts'
path_command = """setenv PATH /cbio/donnees/nvaroquaux/.local/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/cbio/donnees/nvaroquaux/.local/lib
setenv LIBRARY_PATH /cbio/donnees/nvaroquaux/.local/lib
setenv PYTHONPATH /cbio/donnees/aschoenauer/workspace2/cecog/pysrc:/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc
setenv DRMAA_LIBRARY_PATH /opt/gridengine/lib/lx26-amd64/libdrmaa.so
"""
pbsOutDir = '/cbio/donnees/aschoenauer/PBS/OUT'
pbsErrDir = '/cbio/donnees/aschoenauer/PBS/ERR'
pbsArrayEnvVar = 'SGE_TASK_ID'

def globalSummaryScript(baseName, expFile, nb_exp, iter_,bins_type):
    
    jobCount = 0
    i=0
    head = """#!/bin/sh
cd %s""" %progFolder
    baseName = baseName+'{}'.format(bins_type)
#A. DEALING WITH EXPERIMENTS
    for nb in nb_exp:
        for iteration in range(iter_):
            i+=1; jobCount +=1
            cmd ="""
python tracking/trajPack/feature_cell_extraction.py --experimentFile %s --nb_exp %i --iter %i --bins_type %s
"""
            cmd%=(
                  expFile, nb, iteration, bins_type
                  )
            
            # this is now written to a script file (simple text file)
            # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
            script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i))
            script_file = file(script_name, "w")
            script_file.write(head + cmd)
            script_file.close()
    
            # make the script executable (without this, the cluster node cannot call it)
            os.system('chmod a+x %s' % script_name)
    

            # write the main script
    array_script_name = '%s.sh' % os.path.join(scriptFolder, baseName)
    main_script_file = file(array_script_name, 'w')
    main_content = """#!/bin/sh
%s
#$ -o %s
#$ -e %s
%s$%s.sh
""" % (path_command,
       pbsOutDir,  
       pbsErrDir, 
       os.path.join(scriptFolder, baseName),
       pbsArrayEnvVar)

    main_script_file.write(main_content)
    main_script_file.close()
    os.system('chmod a+x %s' % array_script_name)
    sub_cmd = 'qsub -t 1-%i %s' % (jobCount, array_script_name)

    print sub_cmd
    
    return 1

if __name__ == '__main__':
    description =\
'''
%prog - Script generating for experiment summarization, summary clustering and hit finder
'''

    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    parser.add_option("-b", "--base_name", dest="baseName",
                      help="Base name for script")    
#    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
#    parser.add_option('--bin_size', type=int, dest="bin_size", default=50)    
    
    iter_=1
    expFile = '../data/exp_123etctrl.pkl'

    (options, args) = parser.parse_args()
    
    globalSummaryScript(options.baseName,expFile,
                        [12050], iter_, options.bins_type
                      )
    