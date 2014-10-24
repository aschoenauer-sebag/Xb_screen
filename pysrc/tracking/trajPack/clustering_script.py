import os, pdb
from optparse import OptionParser
import numpy as np

#env variables without the need for R
from util import progFolder, jobSize, scriptFolder, path_command, pbsOutDir, pbsArrayEnvVar, pbsErrDir

def generationTrsport(baseName,algo,simulated, debut, fin, only_dataprep,
            div_name, lambda_, datasize, weight, iter_, 
            batch_size, n_init, init, bins_type, bin_size, cost_type, ddim, preprocessed):
    #options.baseName,options.algo, options.data, int(options.n_debut), int(options.n_fin), 
    #options.div_name, options.lambda_, options.size,options.weights, options.iter,\
    #options.batch_size, options.n_init, options.init, options.bins_type, options.bin_size
    jobCount = 0
    
    if debut>0:
        head = """#!/bin/sh
cd %s""" %progFolder
        nb_jobs = fin-debut
        baseName = baseName+'{}_w{}_f{}_{}_{}_{}_a{}'.format(div_name[:5], weight, iter_, init, bins_type, bin_size, algo)
        for i in range(nb_jobs):
            jobCount += 1
            cmd = ''
            # command to be executed on the cluster
            temp_cmd = """
python tracking/histograms/k_means_transportation.py --only_dataprep %i --sim %i --ddimensional %i --preprocessed %i --iter %i -a %i --bins_type %s --cost_type %s --bin_size %i --div_name %s -k %i -w %i --init %s -s %i --batch_size %i
"""
            temp_cmd %= (
                         only_dataprep,
                         simulated,
                         ddim,
                         preprocessed,
                         iter_,
                         algo,
                         bins_type,
                         cost_type,
                         bin_size,
                         div_name,
                         debut+i,
                         weight,
                         init,
                         datasize,
                         batch_size
                    )

            cmd += temp_cmd

            # this is now written to a script file (simple text file)
            # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
            script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i+1))
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
    
        # the submission commando is:
        #sub_cmd = 'qsub -o %s -e %s -t 1-%i %s' % (self.oBatchSettings.pbsOutDir,  
        #                                           self.oBatchSettings.pbsErrDir, 
        #                                           jobCount, array_script_name)
        sub_cmd = 'qsub -t 1-%i %s' % (jobCount, array_script_name)

    else:
        baseName = baseName+'{}_w{}_f{}_{}_{}_{}_a{}'.format(div_name[:5],weight, iter_, init, bins_type, bin_size, algo)
        array_script_name = '%s.sh' % os.path.join(scriptFolder, baseName)
        main_script_file = file(array_script_name, 'w')
        main_content = """#!/bin/sh
%s
#$ -o %s
#$ -e %s

cd %s
""" % (path_command,
           pbsOutDir,  
           pbsErrDir,
           progFolder 
           )
    
        main_script_file.write(main_content)
        temp_cmd = """
python tracking/histograms/k_means_transportation.py --only_dataprep %i --sim %i --ddimensional %i --preprocessed %i --iter %i -a %i --bins_type %s --cost_type %s --bin_size %i --div_name %s -k %i -w %i --init %s -s %i --batch_size %i
"""
        temp_cmd %= (
                     only_dataprep,
                     simulated,
                     ddim,
                     preprocessed,
                     iter_,
                     algo,
                     bins_type,
                     cost_type,
                     bin_size,
                     div_name,
                     debut,
                     weight,
                     init,
                     datasize,
                     batch_size
                    )
        main_script_file.write(temp_cmd)
        main_script_file.close()
        os.system('chmod a+x %s' % array_script_name)
        sub_cmd = 'qsub %s' % (array_script_name)
    #print 'array containing %i jobs' % jobCount
    print sub_cmd
    return 1


def generationScript(baseName,outputname='halfM_max_05', simulated=0, algo=None, data=None, debut=None, fin=None, neighbours=None, sigma=None, density=None, covar=None, fuzzifier=None, num_samp=None):
    jobCount = 10
#    nb_jobs = fin-debut+1
#    if algo==0:
#        baseName = baseName+'_n{}_'.format(neighbours)
#    if algo==4:
#        baseName = baseName+'_n{}_s{}'.format(neighbours, sigma)
    head = """#!/bin/sh
cd %s""" %progFolder
    for i in range(jobCount):
        cmd = ''
        # command to be executed on the cluster
        temp_cmd = """
python tracking/trajPack/clustering.py --action clustering --outputname %s -n %i --simulated %i
"""
        temp_cmd %= (
                     outputname,
                     i,
                     simulated
                )

        cmd += temp_cmd

        # this is now written to a script file (simple text file)
        # the script file is called ltarray<x>.sh, where x is 1, 2, 3, 4, ... and corresponds to the job index.
        script_name = os.path.join(scriptFolder, '%s%i.sh' % (baseName, i+1))
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
    parser.add_option("-b", "--base_name", dest="baseName",
                      help="Base name for script")
    parser.add_option("-i", "--nb_debut", type=int, dest="n_debut", default = 0,
                      help="First index for naming files")
    parser.add_option("-l", "--nb_fin", type=int, dest="n_fin", default = 50,
                      help="Last index for naming files")
    parser.add_option('--only_dataprep', type=int, dest='only_dataprep', default=0)

    parser.add_option("-a", "--algo", dest="algo", default=1,type=int,
                      help="Chosen algo - 1 for mini-batch KMeans")
    parser.add_option("-d", "--data", dest="data", default = 0, type=int,
                      help="Quelles donnees, 0 for real data, 1 for simulated data")
    parser.add_option("-w", type=int, dest="weights", default=0,
                      help="Parameter to choose weights for Sinkhorn distance")
    parser.add_option('--lambda', dest='lambda_', type=int, default=10,
                      help ='Parameter for Sinkhorn distance')
    parser.add_option('--iter', dest='iter', type=int, default=0,
                      help="To choose a random sample of trajectories from workspace2/Tracking/histogramsNtotfirst.pkl")
    parser.add_option('--div_name', type=str, dest='div_name', default='transportation')
    parser.add_option('--ddimensional', type=int, dest='ddim', default=0)
    parser.add_option('--preprocessed', type=int, dest='preprocessed', default=0)
    parser.add_option("-s", dest="size", type=int,default=1000)
    parser.add_option("--batch_size", dest="batch_size", type=int,default=1000)
    parser.add_option("--n_init", dest="n_init", type=int,default=5)
    parser.add_option("--init", dest="init", type=str,default='k-means++')
    parser.add_option('--bins_type', type=str, dest="bins_type", default='quantile')
    parser.add_option('--cost_type', type=str, dest="cost_type", default='number')
    parser.add_option('--bin_size', type=int, dest="bin_size", default=10)
    #dataFolder = "/cbio/donnees/aschoenauer/data/tracking/migration"
    (options, args) = parser.parse_args()
##VERIFIER PUITS PAR PUITS QUI EST LA ET QUI N'EST PAS LA...
#    processedPlates = os.listdir("/cbio/donnees/aschoenauer/data/tracking/migration")
#    toProcessed=[]; p=processedPlates[0]
#    i=1;
#    while i<len(processedPlates):
#        if not np.any(['tabFeatures' in fi for fi in os.listdir(os.path.join(dataFolder, p))]):
#            toProcessed.append(p)
#        if len(toProcessed)==options.N:
#            break
#        p=processedPlates[i]; i+=1
    #toProcessed = processedPlates[:options.N]

    generationTrsport(options.baseName,options.algo, options.data, int(options.n_debut), int(options.n_fin), options.only_dataprep, 
                      options.div_name, options.lambda_, options.size,options.weights, options.iter,\
                      options.batch_size, options.n_init, options.init, options.bins_type, options.bin_size, options.cost_type,
                      options.ddim, options.preprocessed)
    
#SCRIPT POUR LANCER LA GENERATION DES SCRIPTS DES JOBS 
#   for bin_type in ['quantile', 'minmax']:
#        for bin_size in [10,50]:
#            for w in weights:
#                %run trajPack/clustering_script -b totvarNS -d 0 --iter 0 -a 1 --bins_type $bin_type --bin_size $bin_size -i -1 --div_name total_variation -w $w -s 40000 --init k-means++ --batch_size 1000
#
#    
    
    
    
    
#    else:
#        generationScript(options.baseName, int(options.algo), int(options.data), int(options.n_debut), int(options.n_fin), int(options.neighbours), options.sigma, 
#                     int(options.dens), options.covar, int(options.fuzzifier), int(options.num_samp))
    
    
#    parser.add_option("-g", "--neighbours", dest="neighbours", default = 15,
#                      help="Number of neighbours")
#    parser.add_option("-s", "--sigma", dest="sigma", default=0.1, type=float,
#                      help="Parameter sigma for spectral clustering")
#    parser.add_option("--density", dest="dens",type=int, default = 5, 
#                      help="Max nb of neighbours during track")
#    parser.add_option("--covariance", dest="covar",type=str, default = 'full', 
#                      help="Covariance type for GMM algo")
#    parser.add_option("--fuzzy", dest="fuzzifier",type=int, default = 2, 
#                      help="Fuzzifier parameter for Fuzzy CMeans algo")
#    
#    parser.add_option("--numsampling", dest="num_samp",type=int, default = 10, 
#                      help="Number of iterations for Gibbs sampling to estimate parameters for IBP")
    
    
