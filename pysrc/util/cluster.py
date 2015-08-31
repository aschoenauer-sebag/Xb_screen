import os

from util.listFileManagement import strToTuple
from util import progFolder, scriptFolder, pbsArrayEnvVar, pbsErrDir, pbsOutDir
from tracking.trajPack.tracking_script import path_command


def scriptCommand(exp_list, baseName='comp_track', command="tracking/trajPack/cell_cycle.py",jobSize=10,
                  h5_result_dir="/share/data20T/mitocheck/Alice/results",
                  max_nodes=500,
                   **kwargs):
    perExperiment=False
    if type(exp_list[0])!=int:
        perExperiment=True
        
        if len(exp_list[0])!=2:
            exp_list=strToTuple(exp_list, os.listdir(h5_result_dir))
    
    fileNumber = int(len(exp_list)/float(jobSize))+1
    
    head = """#!/bin/sh
cd %s""" %progFolder
    for keyword in kwargs:
        command = command +" -{} {}".format(keyword, kwargs[keyword])
    
    for k in range(fileNumber):
        cmd = ''
        for exp in exp_list[jobSize*k:jobSize*(k+1)]:
            if perExperiment:
                pl,w=exp
                temp_cmd = """
python %s -p %s -w %s"""
                temp_cmd %= (
                        command,
                        pl,
                        w
                        )
            else:
                temp_cmd = """
python %s -i %i"""
                temp_cmd %= (
                        command,
                        exp
                        )
    
            cmd += temp_cmd
    
        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, '%s%i.sh' % (baseName, k+1))
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

    sub_cmd = 'qsub -tc %i -t 1-%i %s' % (max_nodes, fileNumber, array_script_name)

    print sub_cmd
    return 1