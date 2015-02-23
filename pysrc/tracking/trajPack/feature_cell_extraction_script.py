import os, pdb, time
from optparse import OptionParser
import cPickle as pickle
from analyzer import plates
from tracking.histograms.summaries_script import jobSize, progFolder, scriptFolder, path_command, pbsArrayEnvVar, pbsErrDir, pbsOutDir
# path_command = """setenv PATH /cbio/donnees/nvaroquaux/.local/bin:/cbio/donnees/twalter/software/bin:${PATH}
# setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/cbio/donnees/nvaroquaux/.local/lib:/cbio/donnees/twalter/software/lib64/R/lib:/cbio/donnees/twalter/software/lib
# setenv LIBRARY_PATH /cbio/donnees/nvaroquaux/.local/lib:/cbio/donnees/twalter/software/lib:/cbio/donnees/twalter/software/lib64/R/lib
# setenv PYTHONPATH /cbio/donnees/aschoenauer/workspace2/cecog/pysrc:/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc
# setenv DRMAA_LIBRARY_PATH /opt/gridengine/lib/lx26-amd64/libdrmaa.so
# """
# 
# pbsOutDir = '/cbio/donnees/aschoenauer/PBS/OUT'
# pbsErrDir = '/cbio/donnees/aschoenauer/PBS/ERR'
# pbsArrayEnvVar = 'SGE_TASK_ID'

def script_hierarchical_clustering(outFolder='../scripts', baseName='hierclust', outputname='halfM_max_05'):
    cmd ="""
python util/hierarchical_clustering.py --level %f --outputname %s
"""
    head = """#!/bin/sh
cd %s""" %progFolder
    levels = [0.4, 0.5,0.6,0.7,0.8]
    for k, level in enumerate(levels):
        cour_cmd= cmd%(level,outputname)
        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(k+1))
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
    sub_cmd = 'qsub -t 1-%i %s' % (len(levels), array_script_name)

    print sub_cmd
    
    return 1

def script_usable(outFolder='../scripts', baseName='usable', plates=True):
    if plates:
        range_=os.listdir('/share/data20T/mitocheck/tracking_results')
        cmd="""
python util/listFileManagement.py -p %s
"""
    else:
        range_=range(1000)
        cmd ="""
python util/listFileManagement.py -s %i
"""

    head = """#!/bin/sh
cd %s""" %progFolder
    i=1
    for k in range_:
        cour_cmd= cmd%k        
        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i))
        i+=1
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
    sub_cmd = 'qsub -t 1-%i %s' % (len(range_), array_script_name)

    print sub_cmd
    
    return 1
    

def globalSummaryScript(baseName, siRNAFile,div_name,iter, bins_type, bin_size, testCtrl=False, type="simulated"):
    
    jobCount = 0
    i=0
    head = """#!/bin/sh
cd %s""" %progFolder
    baseName = baseName+'{}{}'.format(div_name, iter)
    
    if type=="XBSC":
        settings_file = 'tracking/settings/settings_feature_extraction_XBSC.py'
        
        #Dealing with controls for the first iteration
        if iter==0:
            cmd ="""
python tracking/trajPack/feature_cell_extraction.py --testCtrl %s --solvent plate --div_name %s --settings_file %s --iter 0 --time_window 0
"""         
            for plate in plates:
                print i,
                i+=1; jobCount +=1
                cour_cmd=cmd%(plate, div_name, settings_file)
                script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i))
                script_file = file(script_name, "w")
                script_file.write(head + cour_cmd)
                script_file.close()
        
                # make the script executable (without this, the cluster node cannot call it)
                os.system('chmod a+x %s' % script_name)
        num_ctrl=i
        #Then doing the others
        f=open(siRNAFile, 'r')
        siRNAList = pickle.load(f); f.close()
        cmd ="""
python tracking/trajPack/feature_cell_extraction.py --siRNA %s --div_name %s --solvent plate --settings_file %s --iter %i
"""
        for siRNA in siRNAList:
            print i,
            i+=1; jobCount +=1
            
            cour_cmd= cmd%(
                  siRNA, div_name, settings_file, iter
                  )
            script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i))
            script_file = file(script_name, "w")
            script_file.write(head + cour_cmd)
            script_file.close()
    
            # make the script executable (without this, the cluster node cannot call it)
            os.system('chmod a+x %s' % script_name)
        
    else:
        if type=="simulated":
            settings_file = 'tracking/settings/settings_feature_extraction_SIMULATED_DATA.py'
            folder = '../resultData/simulated_traj/simres/plates'
        else:
            settings_file = 'tracking/settings/settings_feature_extraction.py'
            folder = '/share/data20T/mitocheck/tracking_results'
        
        if not testCtrl:
            #A. DEALING WITH EXPERIMENTS
            
            f=open(siRNAFile, 'r')
            siRNAList = pickle.load(f); f.close()
            cmd ="""
python tracking/trajPack/feature_cell_extraction.py --siRNA %s --div_name %s --bins_type %s --bin_size %s --settings_file %s --iter %i
"""
    
        else:
            #A. DEALING WITH CONTROLS, for a specific list of plates
            baseName+='CTRL_'
            siRNAList = os.listdir(folder)
    
            cmd ="""
python tracking/trajPack/feature_cell_extraction.py --testCtrl %s --div_name %s --bins_type %s --bin_size %s --settings_file %s --iter %i
"""
    
        for siRNA in siRNAList:
            print i,
            i+=1; jobCount +=1
            
            cour_cmd= cmd%(
                  siRNA, div_name, bins_type, bin_size, settings_file, iter
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

Options:
- baseName : script filename
- div_name : possible values: KS or MW
- bins_type : quantile or minmax, in case of transportation or etransportation distances
- bin_size : a number between 5 and 100 (up to you), in case of transportation or etransportation distances

-testCtrl : precise 1 if producing control pvalues. In this case, the iteration variable is not used (only for experiments)
-iter : 0 to 4 at least, to test an experiment against a randomly chosen set of controls of the same plate
-simulated: if working with simulated data
- siRNAFile : the file indicating the correct siRNAs to study (still targeted to exactly one gene according to the state of knowledge)
'''

    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    parser.add_option("-b", "--base_name", dest="baseName",help="Base name for script")
    parser.add_option('--div_name', type=str, dest='div_name', default='KS')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')#possible values: quantile or minmax
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)    
    parser.add_option('--testCtrl', type=str, dest='testCtrl', default=0)   
    parser.add_option('--iter', type=int, dest='iter_num', default=5)
    parser.add_option("-t", "--type", dest="type", default = "simulated", type=str, 
                      help="Use of simulated trajectories or XB SC trajectories or Mitocheck trajectories")
    
    parser.add_option('--siRNA', type=str, dest='siRNAFile', default='../data/siRNA_Simpson.pkl')

    (options, args) = parser.parse_args()
    for k in range(options.iter_num):
        globalSummaryScript(options.baseName,options.siRNAFile, options.div_name, k,
                        options.bins_type, options.bin_size, options.testCtrl, options.type
                      )
    
#NB before Bioimage Informatics, siRNAs impacted by the bug siRNAList = ['123438', '123438', '148427', '148427', '122325', '122325', '28902', '28902']
#and the controls that were modified 
#["LT0061_05--ex2006_03_08--sp2005_06_02--tt17--c5",
#"LT0061_03--ex2005_07_15--sp2005_06_02--tt173--c4",
#"LT0004_47--ex2005_04_08--sp2005_3_18--tt18--c2",
#"LT0004_06--ex2005_11_18--sp2005_03_18-tt173--c3",
#"LT0079_13--ex2005_07_20--sp2005_07_04--tt17--c4",
#"LT0079_05--ex2006_03_10--sp2005_07_14--tt17--c5"]
    