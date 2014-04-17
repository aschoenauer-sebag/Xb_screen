import os, pdb
from optparse import OptionParser

jobSize = 10
progFolder = '/cbio/donnees/aschoenauer/workspace2/Tracking/src'
scriptFolder = '/cbio/donnees/aschoenauer/data/tracking/scripts'
path_command = """setenv PATH /cbio/donnees/nvaroquaux/.local/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/cbio/donnees/nvaroquaux/.local/lib
setenv LIBRARY_PATH /cbio/donnees/nvaroquaux/.local/lib
setenv PYTHONPATH /cbio/donnees/aschoenauer/workspace2/cecog/pysrc:/cbio/donnees/aschoenauer/workspace2/Tracking/src
setenv DRMAA_LIBRARY_PATH /opt/gridengine/lib/lx26-amd64/libdrmaa.so
"""
pbsOutDir = '/cbio/donnees/aschoenauer/PBS/OUT'
pbsErrDir = '/cbio/donnees/aschoenauer/PBS/ERR'
pbsArrayEnvVar = 'SGE_TASK_ID'

def scriptCible(param_list, dataFolder, baseName):
    jobCount = 0
    nb_jobs = len(param_list)
    
    head = """#!/bin/sh
cd %s""" %progFolder
    for plate, well in param_list:
        jobCount += 1
        cmd = ''
        # command to be executed on the cluster
        temp_cmd = """
python trajPack/parallel_trajFeatures.py -p %s -w %s -d %s """

        temp_cmd %= (
                plate,
                well,
                dataFolder
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

def generationScript(platesList, dataFolder, baseName, density):
    jobCount = 0
    nb_jobs = int(len(platesList)/float(jobSize))+1
    
    head = """#!/bin/sh
cd %s""" %progFolder
    for i in range(nb_jobs):
            jobCount += 1
            cmd = ''
            lstJobPositions = [platesList[k] for k in range(i*jobSize, min(len(platesList), (i+1)*jobSize))]

            for plate in lstJobPositions:
                try:
                    wells = os.listdir(os.path.join(dataFolder, plate, 'hdf5'))
                except OSError:
                    print 'no hdf5 files for plate {}'.format(plate)
                    continue
                else:
                    print "plate ", plate
#                #check if plate has already been computed
#                dataFolder =self.oBatchSettings.baseOutDir#"/cbio/donnees/aschoenauer/data/mitocheck/results/"
#                if 'hdf5' in os.listdir(os.path.join(dataFolder, plate)) and pos+'.hdf5' in os.listdir(os.path.join(dataFolder, plate, 'hdf5')):
#                    continue
                for w in wells:
                    # command to be executed on the cluster
                    temp_cmd = """
    python trajPack/parallel_trajFeatures.py -p %s -w %s -d %s -s %f """
    
                    temp_cmd %= (
                            plate,
                            w,
                            dataFolder,
                            density
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

if __name__ == '__main__':
    description =\
'''
%prog - Script generating for tracking
'''

    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)

    parser.add_option("-n", "--nb_plates", type=int, dest="N",
                      help="Number of plates to treat")
    parser.add_option("-b", "--base_name", dest="baseName",
                      help="Base name for script")
#    parser.add_option("-t", "--three_d", dest="TroisD",
#                      help="If you want plots in 3D True, else False")
    parser.add_option("-d", "--data_folder", dest="dataFolder",
                      help="Give absolute path to data")
#    parser.add_option("-s", "--density", dest="density", default=40, type=float,
#                      help="Give absolute path to data")


    (options, args) = parser.parse_args()

    allPlates = os.listdir(options.dataFolder)#pr l'instant "/cbio/donnees/aschoenauer/data/mitocheck/migration"
    processedPlates = os.listdir("/cbio/donnees/aschoenauer/data/tracking/migration")
    #len(allPlates): 497 pour /cbio/donnees/aschoenauer/data/mitocheck/results/migration
    #len(processedPlates) : 257 au 1/8/13
    #je decide de refaire l'extraction des trajectoires AVEC les evaluation de densite pour voir
    
    for i, density in enumerate([50, 60, 80, 100]):
        baseName=options.baseName+'_d{}'.format(density)      
        toProcessed = filter(lambda x: x in processedPlates, allPlates)[:100]

        generationScript(toProcessed, options.dataFolder, baseName, density)
    
    
    
    
    
