import os, pdb, time
from optparse import OptionParser
from collections import Counter
import numpy as np
import cPickle as pickle

from util.listFileManagement import expSi, appendingControl, strToTuple, countingHDF5Done

jobSize = 10
progFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc'
scriptFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/scripts'
path_command = """setenv PATH /cbio/donnees/nvaroquaux/.local/bin:/cbio/donnees/twalter/software/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/cbio/donnees/nvaroquaux/.local/lib:/cbio/donnees/twalter/software/lib64/R/lib:/cbio/donnees/twalter/software/lib
setenv LIBRARY_PATH /cbio/donnees/nvaroquaux/.local/lib:/cbio/donnees/twalter/software/lib:/cbio/donnees/twalter/software/lib64/R/lib
setenv PYTHONPATH /cbio/donnees/aschoenauer/workspace2/cecog/pysrc:/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc:/cbio/donnees/aschoenauer/workspace2/interface_screen
setenv DRMAA_LIBRARY_PATH /opt/gridengine/lib/lx26-amd64/libdrmaa.so
setenv DJANGO_SETTINGS_MODULE interface_screen.settings
"""
pbsOutDir = '/cbio/donnees/aschoenauer/PBS/OUT'
pbsErrDir = '/cbio/donnees/aschoenauer/PBS/ERR'
pbsArrayEnvVar = 'SGE_TASK_ID'

quality_control_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/qc_export.txt'
data_folder = "/share/data20T/mitocheck/tracking_results"

def globalSummaryScript(baseName,  siRNAFile,
                        n_clusters_min, n_clusters_max,
                       div_name,  lambda_,  weights, 
                       bins_type,  bin_size,  cost_type,
                       batch_size,  n_init,  init, 
                       ddim, iter_=0):
    
    f=open(siRNAFile, 'r')
    siRNAList = pickle.load(f); f.close()
    
    siExpDict = expSi(qc = quality_control_file , sens=0)
    jobCount = 0
    i=0
    total_expList = []
    head = """#!/bin/sh
cd %s""" %progFolder
    baseName = baseName+'{}{}_w{}_{}_{}_{}'.format(iter_,div_name[:5], weights, bins_type, bin_size, cost_type)
#A. DEALING WITH EXPERIMENTS
    for siRNA in siRNAList:
        try:
            expList = siExpDict[siRNA]
        except KeyError:
            print "siRNA not in siRNA-experiment dictionary"
        else:
            expList = strToTuple(expList, os.listdir(data_folder))
            total_expList.extend(expList)
            for plate, well in expList:        
                jobCount += 1; i+=1
                cmd = plateWellSummaryScript(plate, well, div_name, lambda_, weights, bins_type, bin_size, cost_type, batch_size, n_init, init, ddim, iter_)

                # this is now written to a script file (simple text file)
                # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
                script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i))
                script_file = file(script_name, "w")
                script_file.write(head + cmd)
                script_file.close()
        
                # make the script executable (without this, the cluster node cannot call it)
                os.system('chmod a+x %s' % script_name)
    
#B. DEALING WITH CONTROLS
    ctrlExp = appendingControl(total_expList)
    ctrlExp = countingDone(ctrlExp)
    np.random.shuffle(ctrlExp)
    ctrlExp=ctrlExp[:int(0.2*len(total_expList))]
    for plate, well in ctrlExp:
        jobCount += 1; i+=1
        cmd = plateWellSummaryScript(plate, well, div_name, lambda_, weights, bins_type, bin_size, cost_type, batch_size, n_init, init, ddim, iter_)

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
    
#C. DOING EXPERIMENT CLUSTERING STEP
    expFilename = 'exp_Simpson_{}.pkl'.format(int(time.time()))
    total_expList.extend(ctrlExp)
    f=open(expFilename, 'w')
    pickle.dump(total_expList, f)
    f.close()
    baseName = baseName+'_clustering'
    for n_clusters in range(n_clusters_min, n_clusters_max):
        script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(n_clusters-n_clusters_min))
        script_file = file(script_name, "w")
        cmd="""
    python tracking/histograms/summarization_clustering.py -a clustering --experimentFile %s -k %i --ddimensional %i --bins_type %s --cost_type %s --bin_size %i --div_name %s -w %i --init %s --batch_size %i --iter %i
    """
        cmd %= (
                expFilename,
                n_clusters,
                 ddim,
                 bins_type,
                 cost_type,
                 bin_size,
                 div_name,
                 weights,
                 init,
                 batch_size,
                 iter_
            )
        script_file.write(head + cmd)
        script_file.close()
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
    sub_cmd = 'qsub -hold_jid  -t 1-%i %s' % (n_clusters_max - n_clusters_min, array_script_name)

    print sub_cmd
    
#D. GOING BACK TO EXPERIMENTS AND TESTING IF DIFFERENT FROM CONTROLS
    
    return 1

def plateWellSummaryScript(plate, well,
                       div_name,  lambda_,  weights, 
                       bins_type,  bin_size,  cost_type,
                       batch_size,  n_init,  init, 
                       ddim, iter_=0):

    # command to be executed on the cluster
    temp_cmd = """
python tracking/histograms/summarization_clustering.py -a summary --plate %s --well %s --ddimensional %i --bins_type %s --cost_type %s --bin_size %i --div_name %s -w %i --init %s --iter %i --batch_size %i
"""
    temp_cmd %= (
                 plate, 
                 well,
                 ddim,
                 bins_type,
                 cost_type,
                 bin_size,
                 div_name,
                 weights,
                 init,
                 iter_,
                 batch_size
            )


    return temp_cmd

def hitFinderScript(baseName, siRNAFile, testCtrl=False, iter_=0):
    f=open(siRNAFile, 'r')
    siRNAList = pickle.load(f); f.close()
    jobCount = 0
    head = """#!/bin/sh
cd %s""" %progFolder
    baseName+='{}'.format(iter_)
    if testCtrl:
        baseName = baseName+'CTRL'
        expList = []
        yqualDict=expSi(quality_control_file, sens=0)
        for siRNA in siRNAList:
            expList.extend(strToTuple(yqualDict[siRNA], os.listdir(data_folder)))
        plates = Counter(np.array(expList)[:,0]).keys()
        for i,plate in enumerate(plates):
            jobCount+=1; i+=1
            cmd = '''
python tracking/histograms/summarization_clustering.py -a hitFinder --verbose 0 --testCtrl %s --iter %i
    '''
            cmd%=(plate, iter_)
            # this is now written to a script file (simple text file)
            # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
            script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i))
            script_file = file(script_name, "w")
            script_file.write(head + cmd)
            script_file.close()
    
            # make the script executable (without this, the cluster node cannot call it)
            os.system('chmod a+x %s' % script_name)
    
    else:
        for i,siRNA in enumerate(siRNAList):
            jobCount+=1; i+=1
            cmd = '''
python tracking/histograms/summarization_clustering.py -a hitFinder --siRNA %s --verbose 0 --iter %i
    '''
            cmd%=(siRNA, iter_)
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
        
    return

if __name__ == '__main__':
    description =\
'''
%prog - Script generating for experiment summarization, summary clustering and hit finder
'''

    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    parser.add_option("-b", "--base_name", dest="baseName",
                      help="Base name for script")
    parser.add_option('-a', dest='action', help='Telling which action to perform', default="summary")
    parser.add_option('--siRNA', type=str, dest='siRNAListFile', default=None)
    parser.add_option('--iter', type=int, dest='iter_', default=0)
    parser.add_option('--testCtrl', type=int, dest='testCtrl', default = 0)
    
    parser.add_option('--nmin', type=int, dest='n_clusters_min',default=4, help="Min number of clusters to use in the experiment clustering step")
    parser.add_option('--nmax', type=int, dest='n_clusters_max',default=13, help="Max number of clusters to use in the experiment clustering step")
    
    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')#possible values: number or value
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)
    parser.add_option('--ddimensional', type=int, dest='ddim', default=0)
    
    parser.add_option("-w", type=int, dest="weights", default=0)
    parser.add_option("-l",type=int, dest="lambda_", default=10)
    parser.add_option("--batch_size", dest="batch_size", type=int,default=1000)
    parser.add_option("--init", dest="init", type=str,default='k-means++')
    parser.add_option("--n_init", dest="n_init", type=int,default=5)
    parser.add_option("--verbose", dest="verbose", type=int,default=0)
    
    
    (options, args) = parser.parse_args()
    if options.action == 'hitFinder':
        if int(options.testCtrl)==0:
            testCtrl=0
        else:
            testCtrl=options.testCtrl
        hitFinderScript(options.baseName, options.siRNAListFile, iter_=options.iter_, testCtrl = testCtrl)
    else:
        globalSummaryScript(options.baseName, options.siRNAListFile,
                        options.n_clusters_min, options.n_clusters_max,
                      options.div_name, options.lambda_, options.weights, 
                      options.bins_type, options.bin_size, options.cost_type,
                      options.batch_size, options.n_init, options.init, 
                      options.ddim, options.iter_
                      )
    