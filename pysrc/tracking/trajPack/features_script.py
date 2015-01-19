import os, pdb
from optparse import OptionParser
from math import ceil
import cPickle as pickle
import numpy as np
from collections import Counter, defaultdict
from tracking.trajPack import time_windows

jobSize = 2
progFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc'
scriptFolder = '/cbio/donnees/aschoenauer/data/features_scripts'
path_command = """setenv PATH /cbio/donnees/nvaroquaux/.local/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/cbio/donnees/nvaroquaux/.local/lib
setenv LIBRARY_PATH /cbio/donnees/nvaroquaux/.local/lib
setenv PYTHONPATH /cbio/donnees/aschoenauer/workspace2/cecog/pysrc:/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc
setenv DRMAA_LIBRARY_PATH /opt/gridengine/lib/lx26-amd64/libdrmaa.so
setenv DJANGO_SETTINGS_MODULE interface_screen.settings
"""
pbsOutDir = '/cbio/donnees/aschoenauer/PBS/OUT'
pbsErrDir = '/cbio/donnees/aschoenauer/PBS/ERR'
pbsArrayEnvVar = 'SGE_TASK_ID'

def scriptTrajFeat(trajDict, featOnly, baseName):
    jobCount = 0
    fileNumber = int(ceil(len(trajDict)/float(jobSize)))
    head = """#!/bin/sh
cd %s""" %progFolder

    for plate in trajDict:
        lstTraj = [(plate, w) for w in trajDict[plate]]
        lstFeat=[(plate, w) for w in trajDict[plate]]
        lstFeat.extend([(plate, w) for w in featOnly[plate]])
        jobCount += 1
        cmd = ''
        for pl, w in lstTraj:
    #FIRST if trajectories do not exist we do them
            temp_cmd = """
    python tracking/trajPack/parallel_trajFeatures.py -p %s -w %s -c %i"""
            temp_cmd %= (
                    pl,
                    w+'.hdf5',
                    0
                    )
    #        print temp_cmd
            cmd += temp_cmd
    #THEN we compute trajectories features anyway
        for pl,w in lstFeat: 
            temp_cmd = """
    python tracking/trajPack/parallel_trajFeatures.py -p %s -w %s -c %i """
            temp_cmd %= (
                    pl,
                    w+'.hdf5',
                    1
                    )
    
            cmd += temp_cmd

        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, '%s%i.sh' % (baseName, jobCount))
        script_file = file(script_name, "w")
        script_file.write(head + cmd)
        script_file.close()

        # make the script executable (without this, the cluster node cannot call it)
        os.system('chmod a+x %s' % script_name)


    for plate in filter(lambda x: x not in trajDict, featOnly):
        lstFeat=[(plate, w) for w in featOnly[plate]]
        jobCount += 1
        cmd = ''
        for pl,w in lstFeat: 
            temp_cmd = """
    python tracking/trajPack/parallel_trajFeatures.py -p %s -w %s -c %i """
            temp_cmd %= (
                    pl,
                    w+'.hdf5',
                    1
                    )
    
            cmd += temp_cmd

        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, '%s%i.sh' % (baseName, jobCount))
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
    os.system('chmod a+x %s' % array_script_name)

    # the submission commando is:
    #sub_cmd = 'qsub -o %s -e %s -t 1-%i %s' % (self.oBatchSettings.pbsOutDir,  
    #                                           self.oBatchSettings.pbsErrDir, 
    #                                           jobCount, array_script_name)
    sub_cmd = 'qsub -t 1-%i -l mem_free=1G %s' % (jobCount, array_script_name)

    print 'array containing %i jobs' % jobCount
    print sub_cmd
    return 1

def scriptCible(param_list, baseName, simulated=False, 
                settings_file='tracking/settings/settings_trajFeatures_XBSC.py', time_window_range=[None]):
    
    jobCount = 0
    fileNumber = int(ceil(len(param_list)/float(jobSize)))
    head = "#!/bin/sh \ncd {}\n" .format(progFolder)

    for i in range(fileNumber):
        lstJobPositions = [param_list[k] for k in range(i*jobSize, min(len(param_list), (i+1)*jobSize))]
        jobCount += 1
        cmd = ''
        for plate, w in lstJobPositions:
#FIRST if trajectories do not exist we do them, if it's not simulated data of course
            if not simulated:
                w+='.hdf5'
                temp_cmd = "python tracking/trajPack/parallel_trajFeatures.py -p {} -w {} -f {} -c {} -s 0\n".format(
                        plate,
                        w,
                        settings_file,
                        0
                        )
        #        print temp_cmd
                cmd += temp_cmd
    #THEN we compute trajectories features anyway
            for time_window in time_window_range:
                temp_cmd = "python tracking/trajPack/parallel_trajFeatures.py -p {} -w {} -f {} -t {} -c {} -s {}\n".format(
                        plate,
                        w,
                        settings_file,
                        time_window,
                        1,
                        int(simulated)
                        )
    
                cmd += temp_cmd

        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, '%s%i.sh' % (baseName, jobCount))
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
    os.system('chmod a+x %s' % array_script_name)

    # the submission commando is:
    #sub_cmd = 'qsub -o %s -e %s -t 1-%i %s' % (self.oBatchSettings.pbsOutDir,  
    #                                           self.oBatchSettings.pbsErrDir, 
    #                                           jobCount, array_script_name)
    sub_cmd = 'qsub -t 1-%i -l mem_free=10G %s' % (jobCount, array_script_name)

    print 'array containing %i jobs' % jobCount
    print sub_cmd
    return 1

if __name__ == '__main__':
    description =\
'''
%prog - Script generating for tracking
'''

    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    parser.add_option("-b", "--base_name", dest="baseName",
                      help="Base name for script")
    parser.add_option("-f", "--settings_file", dest="settings_file", type=str, default='tracking/settings/settings_trajFeatures_XBSC.py',
                      help="Settings_file")
    parser.add_option("-t", "--use_time_windows", dest="use_time_windows", default=0, type=int,
                      help="Choose to use time windows index in settings file or all frames")
    parser.add_option("--traj", dest="trajToDo", type=str, default= None,
                      help="Give absolute path to the file where you have listed the experiments you're interested in")
    parser.add_option("--feat", dest="featToDo", type=str, default= None,
                      help="Give absolute path to the file where you have listed the experiments you're interested in")
    parser.add_option("-s", "--simulated", dest="simulated", default = 0, type=int, 
                      help="Use of simulated trajectories or no")
    parser.add_option("--dataFolder", type=str, dest="dataFolder", 
                       help="In case of using simulated data, indicate where it is")

    (options, args) = parser.parse_args()

    if options.simulated:
        processedPlates = os.listdir(options.dataFolder)
        toProcess=[]
        for plate in processedPlates:
            toProcess.extend([(plate, '{:>03}'.format(k)) for k in range(1,385)])
        scriptCible(toProcess, options.baseName, simulated=True)
        
    elif options.trajToDo is not None:
        f=open(options.trajToDo, 'r')
        trajToDo=pickle.load(f); f.close()
        f=open(options.featToDo, 'r')
        featToDo=filter(lambda x: x not in trajToDo, pickle.load(f)); f.close()
        
        trajToDoDict=defaultdict(list)
        featOnlyDict=defaultdict(list)
        for pl,w in trajToDo:
            trajToDoDict[pl].append(w)
        for pl, w in featToDo:
            featOnlyDict[pl].append(w)
            
        scriptTrajFeat(trajToDoDict, featOnlyDict, options.baseName)
    else:
        f=open(options.featToDo, 'r')
        featToDo=pickle.load(f); f.close()
        if not options.use_time_windows:
            scriptCible(featToDo, options.baseName, settings_file=options.settings_file)
        else:
            scriptCible(featToDo, options.baseName, settings_file=options.settings_file, time_window_range=range(len(time_windows)))
    
    
    
    
    
    
