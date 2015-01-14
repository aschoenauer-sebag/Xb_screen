#Imports
import os, pdb, time, sys, argparse, subprocess, gc
import cPickle as pickle
import numpy as np

from sklearn.cross_validation import KFold
from sklearn.externals.joblib import Parallel
from sklearn.externals.joblib import delayed

import svmapi

from dataPack import treatments, joining, cplexSolving, classify
from importPack import imp, nb_folds, EVENTS, FEATURE_NUMBER, c_list, l_indexes, nb_big_folds

count = 0

def gettingSolu(loadingFolder,allDataFolder):
    global FEATURE_NUMBER
    #var pour indiquer si on prend un toy example ou pas
    puppet = False
    #tableau qui contiendra toutes les features de tlm pr voir lesquelles contiennent des NaN
    tabF= None
    ctrl=7; pheno = 6

    newFrameLot = None
    dataFolder=os.path.join(allDataFolder, 'raw')
    listD = os.listdir(dataFolder)
    for plate in listD:
        print plate
        listW = os.listdir(os.path.join(dataFolder, plate))
        for well in listW:
            well=well[:-5]
            print well
#            from trajPack import d_ctrl
#            if well not in d_ctrl[plate]:
#                if pheno<=0:
#                    print "------------------------PAS PRIS------------------------"
#                    continue
#                else:
#                    print "PHENO"
#                    pheno-=1
#            else:
#                if ctrl<=0:
#                    print "------------------------PAS PRIS------------------------"
#                    continue
#                else:
#                    print "CONTROL"
#                    ctrl-=1
                
            filename = os.path.join(dataFolder, plate,well+".hdf5")
            if puppet:
                filenameT = os.path.join(allDataFolder, 'puppet_trainingset', 'PL'+plate+"___P"+well+"___T00000.xml")
            else:
                filenameT =os.path.join(allDataFolder,'trainingset', 'PL'+plate+"___P"+well+"___T00000.xml")
            
            #ajout du frameLot et du tabF
            frameLotC, tabFC = gettingRaw(filename, filenameT, plate, well)
            if frameLotC==None:
                sys.stdout.write("File {} containing data for plate {}, well {} does not contain all necessary data".format(filename, plate, well))
                continue
            if newFrameLot == None:
                newFrameLot = frameLotC 
            else: newFrameLot.addFrameLot(frameLotC)
            tabF = tabFC if tabF == None else np.vstack((tabF, tabFC))
    
    print "final training set content :"
    count, total= newFrameLot.statisticsTraining2()
    print count, total
    
    #en ce qui concerne le nettoyage des NaN
    c, f = treatments.whichAreNan(tabF)
    print tabF.shape
    featuresToDelete = f.keys()
    
    newFrameLot.clean(featuresToDelete) ##np.delete(X, f, 1)
    FEATURE_NUMBER -=len(featuresToDelete)

    fichier = open(os.path.join(loadingFolder, "featuresToDelete.pkl"), 'w')
    pickle.dump(featuresToDelete, fichier)
    fichier.close()

    print FEATURE_NUMBER
    print "uplets now"
    #ICI ON RECUPERE DONC LES SINGLETS ET DOUBLETS AVEC LA VALEUR DU TRAINING DANS CELL.TO SI ILS Y SONT, NONE SINON
    #POUR LE CENTRE ET LES FEATURES C'EST LA MOYENNE DES OBJETS DU SINGLET
    singlets, doublets = newFrameLot.getTrainingUplets()

    return j(singlets, doublets, FEATURE_NUMBER)

def gettingRaw(filename, filenameT, plate, well, secondary=False,name_primary_channel='primary__primary3', frames_to_skip=None):
    global FEATURE_NUMBER
    print "images loading : plate = "+plate+",well = "+well    
    tabF = None
    try:
        frameLotC, tabF = imp.importRawSegFromHDF5(filename, plate, well, secondary=secondary, name_primary_channel=name_primary_channel, frames_to_skip=frames_to_skip)
    except ValueError:
        sys.stderr.write( sys.exc_info()[1])
        sys.stderr.write("File {} containing data for plate {}, well {} does not contain all necessary data".format(filename, plate, well))
        return None, None

    if filenameT is not None:
        print "training set loading : filename ="+filenameT
        frameLotC.addTraining2(filename, filenameT)
        print "current training set content :"
        print frameLotC.statisticsTraining2()
    
    return frameLotC, tabF

def j(singlets, doublets, featSize, training= True):
    solutions= None
    for plate in singlets:
            print plate
            for well in singlets[plate]:
                print well
                sys.stderr.write("\n plate {}, well {}\n".format(plate, well))
                for index in singlets[plate][well]:
                    print '-- ',
                    if index+1 not in singlets[plate][well] or singlets[plate][well][index+1]==[]:
                        continue
                    singletsL = singlets[plate][well][index]
                    nextSinglets = singlets[plate][well][index+1]
                    doubletsL = doublets[plate][well][index]
                    nextDoublets = doublets[plate][well][index+1]
                    if len(nextSinglets)==1:
                        continue
                    solution = joining.Solution(plate, well, index, singletsL, nextSinglets, doubletsL, nextDoublets, featSize, training)
                    if solutions == None:
                        solutions= joining.Solutions(solution, lstSolutions = None)
                    else:
                        solutions.append(solution)    
    return solutions

def read_examples(filename, sparm):
    """Reads and returns x,y example pairs from a file.
    
    This reads the examples contained at the file at path filename and
    returns them as a sequence.  Each element of the sequence should
    be an object 'e' where e[0] and e[1] is the pattern (x) and label
    (y) respectively.  Specifically, the intention is that the element
    be a two-element tuple containing an x-y pair."""
    #solutions = gettingSolu()#on est obliges, ca marche pas de les passer dans filename pcq il dit qqpart que c est un string
    pkl_file = open(filename, "r")
    solutions = pickle.load(pkl_file)
    pkl_file.close()
    
    result = []
    for sol in solutions.lstSolutions:
        #result.append((sol.featuresMatrix(), sol.truthVec()))
        result.append((sol, sol.truthVec()))
        
    return result

def init_model(sample, sm, sparm):
    """Initializes the learning model.
    
    Initialize the structure model sm.  The sm.size_psi must be set to
    the number of features.  The ancillary purpose is to add any
    information to sm that is necessary from the user code
    perspective.  This function returns nothing."""
    sparm.loss_type = 2
    sparm.slack_norm = 1
    #pdb.set_trace()
    if sparm.c==0:
        sparm.c = 0.01 #SET in the command line
    sparm.epsilon = 0.001 #?? TO SET
    sparm.newconstretrain = 50
    sm.size_psi = sample[0][0].length_feat()*sample[0][0].length_events()
    
def loss(y, ybar, sparm):
    """Return the loss of ybar relative to the true labeling y.
    
    Returns the loss for the correct label y and the predicted label
    ybar.  In the event that y and ybar are identical loss must be 0.
    Presumably as y and ybar grow more and more dissimilar the
    returned value will increase from that point.  sparm.loss_function
    holds the loss function option specified on the command line via
    the -l option.

    The default behavior is to perform 0/1 loss based on the truth of
    y==ybar."""
    #HAMMING LOSS
    number_eventsC = len(y)
    sigma = 0
    prod = 0
    class_weights = [0.96, 0.005, 0.005, 0.01, 0.02]
    for i in range(number_eventsC):
        sigma+=np.sum(y[i])
        #FINAL CHOICE : we are not using the class weights because on average it is less good (accuracy + precision)
        
        #prod +=np.sum( np.array(y[i])*np.array(ybar[i])/float(class_weights[i]))
        prod +=np.sum( np.array(y[i])*np.array(ybar[i]))
    
    return (1 - prod/float(sigma)) 
    
def psi(x, y, sm, sparm):
    """Return a feature vector representing pattern x and label y.

    This is the combined feature function, which this returns either a
    svmapi.Sparse object, or sequence of svmapi.Sparse objects (useful
    during kernel evaluations, as all components undergo kernel
    evaluation separately).  There is no default behavior."""
    #x is a list of the features of all hypotheses of different events, in the order of EVENTS (importPack)
    #y is also a list, of the truth of all hypotheses of different events, same order
    #m in x is np.ndarray of length length_hyp(e), length_feat
    xTilde = x.featuresMatrix()
    number_events = len(xTilde)
    psi = np.array([])
    #pdb.set_trace()
    for m in range(number_events):
        fz = np.dot(np.transpose(xTilde[m]), y[m])
        psi = np.hstack((psi, fz))

    return svmapi.Sparse(psi)    

def find_most_violated_constraint_margin(x, y, sm, sparm):
    """Return ybar associated with x's most violated constraint.

    The find most violated constraint function for margin rescaling.
    The default behavior is that this returns the value from the
    general find_most_violated_constraint function."""
    #dans l'exemple donne il y a une iteration sur toutes les classifications possibles pour trouver la contrainte la plus violee
    #nous on ne peut pas faire ca donc on doit utiliser cplex

    global count
    number_events = len(y)
    length_feat= x.length_feat()
    sigma = 0
    Y = [] ; L = []
    
    #weights :
    w = sm.w
    class_weights = [0.96, 0.005, 0.005, 0.01, 0.02] #[ 0.9597998 ,  0.00417081,  0.00468414,  0.01087619,  0.02046906]
    #FINAL CHOICE : we are not using the class weights because on average it is less good (accuracy + precision)
   
    longueurs = []
    for k in range(number_events):
        sigma+=np.sum(y[k])
        longueurs.append(len(y[k]))
        #Y=np.hstack((Y, np.array(y[k])/float(class_weights[k])))
        Y=np.hstack((Y, np.array(y[k])))
        
    
    Y=Y/sigma#type np.array, taille length_hyp,
    xTilde = x.featuresMatrix()
#SET UP COST FUNCTION
    for k in range(number_events):
        Mf = np.transpose(xTilde[k])
        #print Mf.shape #type :np.ndarray, taille: length_feat, length_hyp(k)
        if Mf.shape==(0,):
            continue
        wC = w[length_feat*k:length_feat*(k+1)]
        #print len(wC) #type : list longueur: length_feat
        L = np.hstack((L, np.dot(wC, Mf)))
    #print L.shape#normalement c'est length_hyp
    costFunc = L-Y
    constraints = x.constraintsMatrix()
    #print constraints.shape#(t1+t2, length_hyp)
    
#FIND CUTTING PLANE
    output2, msg = cplexSolving.solve(costFunc, constraints, x.singletsSize)#c'est ybar qu'il faut retourner
    
    print "----------------", sparm.c, "--------------------------####################MSG CPLEX"
    print msg
    output = [int(j) for j in output2]
    ybar=[];last=0
    for k in range(number_events):
        ybar.append(list(output[last:last+longueurs[k]]))
        last +=longueurs[k]
    count+=1; print "nb de fois que l'on a appele find_most_violated_constraint", count
#    if count>= 745:
#        pdb.set_trace()
#    print "ybar : "
#    for k in range(number_events):
#        print ybar[k], "et la verite :", y[k]
    print "loss :", loss(y, ybar, sparm)
    return ybar

    
def print_iteration_stats(ceps, cached_constraint, sample, sm,
                          cset, alpha, sparm):
    """Called just before the end of each cutting plane iteration.

    This is called just before the end of each cutting plane
    iteration, primarily to print statistics.  The 'ceps' argument is
    how much the most violated constraint was violated by.  The
    'cached_constraint' argument is true if this constraint was
    constructed from the cache.
    
    The default behavior is that nothing is printed."""
    print "###############################################################################################################################BOULOUM BOULOUM", ceps
    print "nouvelle valeur de w apres cette cutting plane iteration :"
    for i in sm.w: print i
    #pdb.set_trace()
    
def write_model(filename, sm, sparm):
    global monOutput
    global count
    """Dump the structmodel sm to a file.
    
    Write the structmodel sm to a file at path filename.

    The default behavior is equivalent to
    'cPickle.dump(sm,bz2.BZ2File(filename,'w'))'."""
#    f = bz2.BZ2File(filename, 'w')
#MARCHE PAS :
    f = open(filename, 'w')
    pickle.dump(sm.w, f)
    f.close()

def sousProcessLearn(c, n_big_fold, n_fold, outputLearn=None, loadingFolder='../results/'):
    #pdb.set_trace()
    #juste pour preciser, ici je peux dire stdout = existing file object => tout le log sera ecrit dans un fichier
    if n_fold ==-1:
        if n_big_fold ==-1:
            subprocess.call(["/media/lalil0u/New/software2/downloads/unpacked/svm-python-v204/svm_python_learn", "--m", "test", "-c", c, "-v", "0", "-y", "0", 
                     os.path.join(loadingFolder,"normalized_data_all.pkl"), os.path.join(loadingFolder,"modelfile_all"+c+".pkl")], 
                    stdout=outputLearn)
        else:
            subprocess.call(["/media/lalil0u/New/software2/downloads/unpacked/svm-python-v204/svm_python_learn", "--m", "test", "-c", c, "-v", "0", "-y", "0", 
                     os.path.join(loadingFolder,"data_TRAIN_fold_BIG"+str(n_big_fold)+".pkl"), 
                     os.path.join(loadingFolder,"modelfile_all"+c+"_"+str(n_big_fold)+".pkl")], stdout=outputLearn)
    else:
        subprocess.call(["/media/lalil0u/New/software2/downloads/unpacked/svm-python-v204/svm_python_learn", "--m", "test", "-c", c, "-v", "0", "-y", "0", 
                     os.path.join(loadingFolder,"data_TRAIN_fold_SMALL"+str(n_big_fold)+'_'+str(n_fold)+".pkl"),
                     os.path.join(loadingFolder,"modelfile_"+c+"_"+str(n_big_fold)+'_'+str(n_fold)+".pkl")], stdout=outputLearn)
    return 1

#def sousProcessClassify(c,n_fold, loadingFolder):
#    #subprocess.call(["/media/lalil0u/New/software2/downloads/unpacked/svm-python-v204/svm_python_classify", "--m", "test", "-v", "3", "results/data_TEST_fold"+str(n_fold)+".pkl", "results/modelfile_"+c+"_"+str(n_fold)+".pkl"])
#    f = open(os.path.join(loadingFolder,"models", "modelfile_"+c+"_"+str(n_fold)+".pkl"), 'r')
#    mesPoids = pickle.load(f)
#    f.close()
#    #pdb.set_trace()
#    test_sol = classify.read_examples(os.path.join(loadingFolder,"data_TEST_fold"+str(n_fold)+".pkl"))
#    total_loss = []; test_contenu = []
#    for (x,y) in test_sol:
#        test_contenu.append(classify.contenu(y))
#        ybar = classify.classify_example(x, mesPoids)
#        l = classify.loss(y, ybar)#attention ce n'est plus la hamming loss ms une liste de cinq chiffres
#        total_loss = l if total_loss == [] else np.vstack((total_loss, l))
#    return total_loss, test_contenu

def wholeProcess(i, nb_big_folds, nb_fold, outputLearn, loadingFolder):
    outputTotal = open(os.path.join(loadingFolder,"outputTotal_C_10**"+str(i)+".txt"), 'w')
    
    #loss_all_folds_events = np.zeros(shape=(len(EVENTS),)); loss_all_folds = 0
    print>>outputTotal,"TIME TIME TIME TIME TIME", time.clock()
    print "number of error folds ", nb_big_folds
    print "number of model selection folds :", nb_fold
    
    if nb_fold==-1:
        print "LEARNING with C = ", 10**i, "and big fold number ", nb_big_folds
        if outputLearn == True:
            output = open(os.path.join(loadingFolder, "outputLearn"+"{:f}".format(10**i)+"_"+str(nb_big_folds)+".txt"), 'w')
        else:
            output = None
        print>>outputTotal,"LEARNING with C = ", 10**i, "and big fold number ", nb_big_folds
        
        try:
            sousProcessLearn("{:f}".format(10**i), nb_big_folds, nb_fold, output, loadingFolder)
        except IOError:
            print "IOError when calling sousProcessLearn, C=", 10**i, "n_big_fold = ", nb_big_folds
            nb_fold-=1

    else:
        for n_big_fold in range(nb_big_folds):
            for n_fold in range(nb_fold):
                print "LEARNING with C = ", 10**i, "and big fold number ", n_big_fold, "and small fold number", n_fold
                if outputLearn == True:
                    output = open(os.path.join(loadingFolder,"outputLearn"+"{:f}".format(10**i)+"_"+str(n_big_fold)+"_"+str(n_fold)+".txt"), 'w')
                else:
                    output = None
                print>>outputTotal,"LEARNING with C = ", 10**i, "and big fold number ", n_big_fold, "and small fold number", n_fold
                try:
                    sousProcessLearn("{:f}".format(10**i), n_big_fold, n_fold, output, loadingFolder)
                except IOError:
                    print "IOError when calling sousProcessLearn, C=", 10**i, "n_fold = ", n_fold
                    nb_fold-=1
                    continue

    outputTotal.close()    
    return 1

def foldingBig(mesSolu, trainT, testT, n_big_fold, small=True):
    if n_big_fold == -1:
        minMax = mesSolu.normalisation()
        output = open(os.path.join(loadingFolder,"minMax_data_all.pkl"), "w")  
        pickle.dump(minMax, output)
        output.close()
    
    else:
        train, test =trainT[n_big_fold], testT[n_big_fold]
        if small:
            kf = KFold(len(train), nb_folds, shuffle = True)
            trainST = []; testST = []
            for tr, te in kf:
                trainST.append(tr); testST.append(te)
            
        training_set = filter(lambda x : mesSolu.lstSolutions.index(x) in train, mesSolu.lstSolutions)
        test_set = filter(lambda x : mesSolu.lstSolutions.index(x) in test, mesSolu.lstSolutions)
        training_sol = joining.Solutions(solution = None, lstSolutions = training_set)
        test_sol = joining.Solutions(solution = None, lstSolutions = test_set)
        
        minMax = training_sol.normalisation()
        #test_sol.normalisation(minMax)
        output = open(os.path.join(loadingFolder,"minMax"+str(n_big_fold)+".pkl"), "w")  
        pickle.dump(minMax, output)
        output.close()
        output = open(os.path.join(loadingFolder,"data_TRAIN_fold_BIG"+str(n_big_fold)+".pkl"), "w")  
        pickle.dump(training_sol, output)
        output.close()
        
        if small:
            training_sol.denormalisation(minMax)
            
            for n_f in range(nb_folds):
                folding(training_sol, trainST, testST, n_f, n_big_fold) 
        
        td = {}
        for sol in test_sol.lstSolutions:
            if sol.plate not in td:
                td[sol.plate]={}
            if sol.well not in td[sol.plate]:
                td[sol.plate][sol.well]=[]
            td[sol.plate][sol.well].append(sol.index)
        output = open(os.path.join(loadingFolder,"data_TEST_fold_BIG"+str(n_big_fold)+".pkl"), "w")  
        pickle.dump(td, output)
        output.close()

    return 1

def folding(mesSolu, trainT, testT, n_fold, n_big_fold):
    if n_fold == -1:
        minMax = mesSolu.normalisation()
        output = open(os.path.join(loadingFolder,"minMax_data_all.pkl"), "w")  
        pickle.dump(minMax, output)
        output.close()
    
    else:
        train, test =trainT[n_fold], testT[n_fold]
        training_set = filter(lambda x : mesSolu.lstSolutions.index(x) in train, mesSolu.lstSolutions)
        test_set = filter(lambda x : mesSolu.lstSolutions.index(x) in test, mesSolu.lstSolutions)
        
        training_sol = joining.Solutions(solution = None, lstSolutions = training_set)
        test_sol = joining.Solutions(solution = None, lstSolutions = test_set)

        minMax = training_sol.normalisation()
        #test_sol.normalisation(minMax)
        output = open(os.path.join(loadingFolder,"minMax"+str(n_big_fold)+'_'+str(n_fold)+".pkl"), "w")  
        pickle.dump(minMax, output)
        output.close()
        output = open(os.path.join(loadingFolder,"data_TRAIN_fold_SMALL"+str(n_big_fold)+'_'+str(n_fold)+".pkl"), "w")  
        pickle.dump(training_sol, output)
        output.close()
        training_sol.denormalisation(minMax)
        
        td = {}
        for sol in test_sol.lstSolutions:
            if sol.plate not in td:
                td[sol.plate]={}
            if sol.well not in td[sol.plate]:
                td[sol.plate][sol.well]=[]
            td[sol.plate][sol.well].append(sol.index)
        output = open(os.path.join(loadingFolder,"data_TEST_fold_SMALL"+str(n_big_fold)+'_'+str(n_fold)+".pkl"), "w")  
        pickle.dump(td, output)
        output.close()

    return 1

if __name__ == '__main__':
    #global loadingFolder
    #pdb.set_trace()
    parser = argparse.ArgumentParser(description = "Structured learning for cell tracking")
    parser.add_argument('-r', '--reloading',default=0, type=int, dest ='reloading')#possible values : 0 to load data ab initio, 1 to load data from /results, s 
    #to load data from /results_small et b to load data from /results_big
    
#    parser.add_argument('-init', nargs = 1, default = -1, type = int, dest='c0')
#    parser.add_argument('-end', nargs = 1, default = 4, type = int, dest='c1')
    parser.add_argument('-o', '--output', default=0, type = int, dest='output')
    parser.add_argument('-w', '--who', default='folds', type = str, dest='training_set')
    parser.add_argument("-d", "--data_folder", default='/media/lalil0u/New/workspace2/Tracking/data/raw', 
                        type=str, dest="dataFolder", help="Give absolute path to data")
    parser.add_argument("-l", "--loading_folder", default='/media/lalil0u/New/workspace2/Tracking/results', 
                        type=str, dest="loadingFolder", help="Give absolute path to folder where to save outputs")
    args = parser.parse_args()

    reloading = args.reloading; #init= args.c0[0] ; end = args.c1[0]; 
    training_set = args.training_set
    loadingFolder=args.loadingFolder
    
    print 'loading folder', loadingFolder
    print 'allDataFolder ', args.dataFolder
    
    print "TIME TIME TIME TIME TIME", time.clock()
    
    if training_set == "all":
        if reloading==0:
        #DATA EXTRACTION
            mesSolu = gettingSolu(loadingFolder, args.dataFolder)

        #DATA NORMALIZATION
            folding(mesSolu, None, None, -1, -1)
        #DATA SAVING
            output = open(os.path.join(loadingFolder,"normalized_data_all.pkl"), "w")
            pickle.dump(mesSolu, output)
            output.close()
        elif reloading == 1:
            ins = open(os.path.join(loadingFolder,"normalized_data_all.pkl"), "r")
            mesSolu = pickle.load(ins)
            ins.close()
            #(c, n_big_fold, n_fold, outputLearn=None, loadingFolder='../results/')
        #Parallel(n_jobs=3, verbose=5)(delayed(sousProcessLearn)("{:f}".format(10**i), -1, -1, args.output, loadingFolder) for i in [1.5])
        sousProcessLearn("{:f}".format(10**1.5), -1, -1, args.output, loadingFolder)
            
    elif training_set == "folds":        
        if reloading==0:
        #DATA EXTRACTION
            mesSolu = gettingSolu(loadingFolder, args.dataFolder)
        #DATA FOLDING          
            
            kf = KFold(len(mesSolu.lstSolutions), nb_big_folds, shuffle = True)
            trainT = []; testT = []
            for train, test in kf:
                trainT.append(train); testT.append(test)
            Parallel(n_jobs=3, verbose=5)(delayed(foldingBig)(mesSolu, trainT, testT, n_big_fold, small=True) for n_big_fold in range(nb_big_folds))

        print "TIME TIME TIME TIME TIME", time.clock()
        print 'liste des puissances de C etudiees : ', c_list
        
        Parallel(n_jobs=3, verbose=5)(delayed(wholeProcess)(i, nb_big_folds, nb_folds, args.output, loadingFolder) for i in c_list)
        
        try:
            del mesSolu; gc.collect()
        except:
            pass
        
    elif training_set == "big_folds":        
        if reloading==0:
        #DATA EXTRACTION
            mesSolu = gettingSolu(loadingFolder, args.dataFolder)
        #DATA FOLDING          
            
            kf = KFold(len(mesSolu.lstSolutions), nb_big_folds, shuffle = True)
            trainT = []; testT = []
            for train, test in kf:
                trainT.append(train); testT.append(test)
            Parallel(n_jobs=3, verbose=5)(delayed(foldingBig)(mesSolu, trainT, testT, n_big_fold, small=False) for n_big_fold in range(nb_big_folds))

        print "TIME TIME TIME TIME TIME", time.clock()
        print 'liste des puissances de C etudiees : ', c_list
        
        
        Parallel(n_jobs=3, verbose=5)(delayed(wholeProcess)(2.5, k, -1, args.output, loadingFolder) for k in range(nb_big_folds))
        
        try:
            del mesSolu; gc.collect()
        except:
            pass


        print "TIME TIME TIME TIME TIME", time.clock()
