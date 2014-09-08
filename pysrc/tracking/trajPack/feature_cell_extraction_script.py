import os, pdb, time
from optparse import OptionParser
import cPickle as pickle


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

def globalSummaryScript(baseName, siRNAFile,div_name, bins_type, bin_size, testCtrl=False):
    
    jobCount = 0
    i=0
    if not testCtrl:
        #A. DEALING WITH EXPERIMENTS
        f=open(siRNAFile, 'r')
        siRNAList = pickle.load(f); f.close()
        cmd ="""
python tracking/trajPack/feature_cell_extraction.py --siRNA %s --div_name %s --bins_type %s --bin_size %s
"""

    else:
        #A. DEALING WITH CONTROLS
        baseName+='CTRL_'
        siRNAList = os.listdir('/share/data20T/mitocheck/compressed_data/')
        cmd ="""
python tracking/trajPack/feature_cell_extraction.py --testCtrl %s --div_name %s --bins_type %s --bin_size %s
"""
    
    head = """#!/bin/sh
cd %s""" %progFolder
    baseName = baseName+'{}{}{}'.format(div_name, bins_type, bin_size)

    for siRNA in siRNAList:
        i+=1; jobCount +=1
        
        cour_cmd= cmd%(
              siRNA, div_name, bins_type, bin_size
              )
        
        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i))
        script_file = file(script_name, "w")
        script_file.write(head + cour_cmd)
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
%prog - Script generating for calculating distances between films
'''

    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    parser.add_option("-b", "--base_name", dest="baseName",
                      help="Base name for script")    
    #possible values: transportation (for Sinkhorn divergence), etransportation (for exact transportation distances in 1D), total_variation or hellinger
    parser.add_option('--div_name', type=str, dest='div_name', default='etransportation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)    
    parser.add_option('--testCtrl', type=str, dest='testCtrl', default=0)

    parser.add_option('--siRNA', type=str, dest='siRNAFile', default='../data/siRNA_Simpson.pkl')

    (options, args) = parser.parse_args()
    
    globalSummaryScript(options.baseName,options.siRNAFile, options.div_name,
                        options.bins_type, options.bin_size, options.testCtrl
                      )
    