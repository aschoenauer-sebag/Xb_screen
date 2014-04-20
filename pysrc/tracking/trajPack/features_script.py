import os, pdb
from optparse import OptionParser
from math import ceil
import cPickle as pickle
import numpy as np
from collections import Counter, defaultdict

jobSize = 10
progFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc'
scriptFolder = '/cbio/donnees/aschoenauer/data/features_scripts'
path_command = """setenv PATH /cbio/donnees/nvaroquaux/.local/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/cbio/donnees/nvaroquaux/.local/lib
setenv LIBRARY_PATH /cbio/donnees/nvaroquaux/.local/lib
setenv PYTHONPATH /cbio/donnees/aschoenauer/workspace2/cecog/pysrc:/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc
setenv DRMAA_LIBRARY_PATH /opt/gridengine/lib/lx26-amd64/libdrmaa.so
"""
pbsOutDir = '/cbio/donnees/aschoenauer/PBS/OUT'
pbsErrDir = '/cbio/donnees/aschoenauer/PBS/ERR'
pbsArrayEnvVar = 'SGE_TASK_ID'

def scriptTrajFeat(trajDict, featOnly, dataFolder, baseName):
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
    python tracking/trajPack/parallel_trajFeatures.py -p %s -w %s -c %i -d %s"""
            temp_cmd %= (
                    pl,
                    w+'.hdf5',
                    0,
                    dataFolder
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
    sub_cmd = 'qsub -t 1-%i %s' % (jobCount, array_script_name)

    print 'array containing %i jobs' % jobCount
    print sub_cmd
    return 1

def scriptCible(param_list, dataFolder, baseName):
    jobCount = 0
    fileNumber = int(ceil(len(param_list)/float(jobSize)))
    head = """#!/bin/sh
cd %s""" %progFolder

    for i in range(fileNumber):
        lstJobPositions = [param_list[k] for k in range(i*jobSize, min(len(param_list), (i+1)*jobSize))]
        jobCount += 1
        cmd = ''
        for plate, w in lstJobPositions:
    #FIRST if trajectories do not exist we do them
            temp_cmd = """
    python tracking/trajPack/parallel_trajFeatures.py -p %s -w %s -c %i -d %s"""
            temp_cmd %= (
                    plate,
                    w+'.hdf5',
                    0,
                    dataFolder
                    )
#    #        print temp_cmd
#            cmd += temp_cmd
    #THEN we compute trajectories features anyway
            temp_cmd = """
    python tracking/trajPack/parallel_trajFeatures.py -p %s -w %s -c %i """
            temp_cmd %= (
                    plate,
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
    sub_cmd = 'qsub -t 1-%i %s' % (jobCount, array_script_name)

    print 'array containing %i jobs' % jobCount
    print sub_cmd
    return 1

#def generationScript(platesList, dataFolder, baseName):
#    jobCount = 0
#    nb_jobs = int(len(platesList)/float(jobSize))
#    
#    head = """#!/bin/sh
#cd %s""" %progFolder
#    for i in range(nb_jobs):
#            jobCount += 1
#            cmd = ''
#            lstJobPositions = [platesList[k] for k in range(i*jobSize, min(len(platesList), (i+1)*jobSize))]
#
#            for plate in lstJobPositions:
##                #check if plate has already been computed
##                dataFolder =self.oBatchSettings.baseOutDir#"/cbio/donnees/aschoenauer/data/mitocheck/results/"
##                if 'hdf5' in os.listdir(os.path.join(dataFolder, plate)) and pos+'.hdf5' in os.listdir(os.path.join(dataFolder, plate, 'hdf5')):
##                    continue
#
#                # command to be executed on the cluster
#                temp_cmd = """
#python trajPack/parallel_trajFeatures.py -p %s -w %s -d %s -c %i """
#
#                temp_cmd %= (
#                        plate,
##                        dataFolder, pas besoin pour le moment
#                        1
#                        )
#
#                cmd += temp_cmd
#
#            # this is now written to a script file (simple text file)
#            # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
#            script_name = os.path.join(scriptFolder, '%s%i.sh' % (baseName, jobCount))
#            script_file = file(script_name, "w")
#            script_file.write(head + cmd)
#            script_file.close()
#
#            # make the script executable (without this, the cluster node cannot call it)
#            os.system('chmod a+x %s' % script_name)
#
#        # write the main script
#    array_script_name = '%s.sh' % os.path.join(scriptFolder, baseName)
#    main_script_file = file(array_script_name, 'w')
#    main_content = """#!/bin/sh
#%s
##$ -o %s
##$ -e %s
#%s$%s.sh
#""" % (path_command,
#       pbsOutDir,  
#       pbsErrDir, 
#       os.path.join(scriptFolder, baseName),
#       pbsArrayEnvVar)
#
#    main_script_file.write(main_content)
#    os.system('chmod a+x %s' % array_script_name)
#
#    # the submission commando is:
#    #sub_cmd = 'qsub -o %s -e %s -t 1-%i %s' % (self.oBatchSettings.pbsOutDir,  
#    #                                           self.oBatchSettings.pbsErrDir, 
#    #                                           jobCount, array_script_name)
#    sub_cmd = 'qsub -t 1-%i %s' % (jobCount, array_script_name)
#
#    print 'array containing %i jobs' % jobCount
#    print sub_cmd
#    return 1

if __name__ == '__main__':
    description =\
'''
%prog - Script generating for tracking
'''

    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    parser.add_option("-b", "--base_name", dest="baseName",
                      help="Base name for script")
#    parser.add_option("-t", "--three_d", dest="TroisD",
#                      help="If you want plots in 3D True, else False")
    parser.add_option("-d", "--data_folder", dest="dataFolder", type=str, default='/share/data20T/mitocheck/Alice/results',
                      help="Give absolute path to tracking data")
    parser.add_option("--traj", dest="trajToDo", type=str, default= None,
                      help="Give absolute path to the file where you have listed the experiments you're interested in")
    parser.add_option("--feat", dest="featToDo", type=str, default= None,
                      help="Give absolute path to the file where you have listed the experiments you're interested in")

    (options, args) = parser.parse_args()
    dataFolder = options.dataFolder#"/cbio/donnees/aschoenauer/data/tracking/migration"
    
    
##here we are going to compute again same plates, subensemble of migration hits to get to know the new features
#    
#    f=open("/cbio/donnees/aschoenauer/workspace2/Tracking/src/trsf_maxd50inf3_data.pkl", 'r'); d=pickle.load(f); f.close()
#    toProcess = d[1]
#VERIFIER PUITS PAR PUITS QUI EST LA ET QUI N'EST PAS LA...
    if options.trajToDo is None and options.featToDo is None:
        processedPlates = os.listdir(dataFolder)#"/cbio/donnees/aschoenauer/data/tracking/migration")
        toProcess=[]
        for plate in processedPlates:
            try:
                liste = os.listdir(os.path.join(dataFolder, plate, 'hdf5'))
            except:
                pass
            else:
                toProcess.extend([(plate, w) for w in liste])
        scriptCible(toProcess, dataFolder, options.baseName)
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
            
        scriptTrajFeat(trajToDoDict, featOnlyDict, dataFolder, options.baseName)
    else:
        f=open(options.featToDo, 'r')
        featToDo=pickle.load(f); f.close()
        scriptCible(featToDo, dataFolder, options.baseName)
    
    
    
    
    
    
