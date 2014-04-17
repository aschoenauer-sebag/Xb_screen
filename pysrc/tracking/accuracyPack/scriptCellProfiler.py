import os, subprocess, pdb
from accuracyCellProfiler import gettingDataIndexes, index_file



if __name__=='__main__':
    path = "/media/lalil0u/New/projects/groundtruth/results/all2"
    #pour enregistrer la segmentation
    pathSave = "/media/lalil0u/New/data/ground_truth_tracking/CecogSegmentation"
    pathData = '/media/lalil0u/New/workspace2/Tracking/data/CeproPredictions'
    nom='accuracy_Cepro'; param =['split_alternative_cost', 'max_split_score', 'merge_alternative_cost', 'max_merge_score']
    

    first_pass=[(40,80),(60,100), (100, 140)]
    second_pass = [(100,300),(150,300), (200, 300)]
#ATTENTION QUAND ON AJOUTE DES PARAMETRES A TESTER IL FAUT AUSSI CREER LA PIPELINE CELL PROFILER CORRESPONDANTE
    third_pass = [(40, 140),(60, 140), (80, 140)]
    first_pass_m = [(40,100),(40,120), (20, 80), (60,80)]
    param_list=[] 
    metadata=gettingDataIndexes(index_file)
    #for j in range(len(first_pass)):
    parameter='m{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (100,120))
    param_list.append(parameter)
    parameter='m{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((60,80), (80,100))
    param_list.append(parameter)
    parameter='m{0[0]}_{0[1]}_s{1[0]}_{1[1]}'.format((70,80), (100,140))
    param_list.append(parameter)
    print "PARAMETER ", param_list
    for parameter in param_list:
        try:
            os.mkdir(os.path.join(pathData, parameter))
        except OSError: 
            print "Folder %s already created"%os.path.join(pathData, parameter)
        finally:
            for p,w in metadata:
                print "plate ", p, "well ", w
                try:
                    l=os.listdir(os.path.join(pathData, parameter, p))
                except OSError:
                    pass
                else:
                    if os.path.join(pathData, parameter, p, w+'.h5') in os.listdir(os.path.join(pathData, parameter, p)):
                        continue
                try:
                    os.mkdir(os.path.join(pathData, parameter, p))
                except OSError: 
                    print "Folder %s already created"%os.path.join(pathData, parameter, p)
                finally:
                    CellProOutput = open(os.path.join(pathData, "output{0}{1}{2}.txt".format(parameter, p, w)), 'w')
                    
                    output_file= os.path.join(pathData, parameter, p, w+'.h5') #ajouter les param et le dossier de la cross validation
                    output_folder = os.path.join(pathData,parameter, p)
                    first, last = metadata[(p,w)]
                    subprocess.call(["/usr/bin/python", "/media/lalil0u/New/workspace2/CellProfiler/CellProfiler.py",\
                                     "-c", "-r", "-i", "/media/lalil0u/New/data/ground_truth_tracking",\
                                     "-o", output_folder, "-p", "/home/lalil0u/pipeline_finale_%s.cp"%parameter, "-f", str(first),\
                                     "-l", str(last), output_file],stdout=CellProOutput)
                    CellProOutput.close()
    print "FINI FINI FINI"
    
    
        #CODE TO MAKE BINARY MASK CORRESPONDING TO CECOG SEGMENTATION
#    for plate in os.listdir(pathSave):
#        print plate
#        for well in os.listdir(os.path.join(pathSave, plate)):
#            print well
#            if os.listdir(os.path.join(pathSave, plate, well)) ==[]:
#                print "vide"
#                pathSearch = path+'/'+plate+'/analyzed/00'+well+'_01'
#                images=[]
#                for dirpath, dirnames, filenames in os.walk(pathSearch):
#                    print dirpath
#                    for filename in filter(lambda x: os.path.splitext(x)[-1]=='.tif', filenames):
#                        images.append((dirpath, filename))
#                         
#                pathC=os.path.join(pathSave, plate, well)
#                for cPath, image in images:
#                    if not makeTif(cPath, pathC, image):
#                        print cPath, pathC, image
#                        raise