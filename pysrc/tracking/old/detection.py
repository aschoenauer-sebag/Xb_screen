import os, pdb, time, sys, gc
import numpy as np
from importPack import imp, EVENTS, FEATURE_NUMBER, nb_folds#, c_list
from dataPack import treatments, joining, cplexSolving, classify
import cPickle as pickle
import accuracy, tracking

def classify_example(x, weights):
    """Given a pattern x, return the predicted label."""

    number_events = len(x.events)
    length_feat= x.length_feat()
    L = []
    
    #weights :
    w = weights

    longueurs = []; M=None
    xTilde = x.featuresMatrix()
#SET UP LIKELIHOOD FUNCTION
    for k in range(number_events):
        Mf = np.transpose(xTilde[k])
        M=Mf if M==None else np.hstack((M, Mf))
      #  print Mf.shape #type :np.ndarray, taille: length_feat, length_hyp(k)
        
        if Mf.shape==(0,):
            longueurs.append(0)
            continue
        longueurs.append(Mf.shape[1])
        wC = w[length_feat*k:length_feat*(k+1)]
     #   print len(wC) #type : list longueur: length_feat
        L = np.hstack((L, np.dot(wC, Mf)))
    #print L.shape#normalement c'est length_hyp
    costFunc = L
    constraints = x.constraintsMatrix()
    #print constraints.shape#(t1+t2, length_hyp)
    
#FIND ARGMAX LIKELIHOOD FUNCTION
    output2, msg = cplexSolving.solve(costFunc, constraints, x.singletsSize)#c'est ybar qu'il faut retourner
    
    print "------------------------------------------####################MSG CPLEX"
    #print msg
    output = [int(j) for j in output2]
    ybar=[];last=0
    for k in range(number_events):
        ybar.append(list(output[last:last+longueurs[k]]))
        last +=longueurs[k]

#    print "ybar : "
#    for k in range(number_events):
#        print ybar[k]
#    print "loss :", loss(y, ybar, sparm)
    return ybar, costFunc, M

def sousProcessClassify(sol, loadingFolder, loadingFile= None, i =None, n_f =None):
    #subprocess.call(["/media/lalil0u/New/software2/downloads/unpacked/svm-python-v204/svm_python_classify", "--m", "test", "-v", "3", "results/data_TEST_fold"+str(n_fold)+".pkl", "results/modelfile_"+c+"_"+str(n_fold)+".pkl"])
    if i==None:
        f = open(loadingFolder+"modelfile_all.pkl", 'r')
    elif n_f==None:
        c = "{:f}".format(10**i)
        f = open(loadingFolder+"modelfile_all"+c+".pkl", 'r')
    else:
        c = "{:f}".format(10**i)
        f = open(loadingFolder+"modelfile_"+c+"_"+str(n_f)+".pkl", 'r')
    mesPoids = pickle.load(f)
    f.close()
    new_sol = []
    #pdb.set_trace()
    if loadingFile == None:
        for solu in sol.lstSolutions:
            new_sol.append((solu, solu.truthVec()))
    else:
        new_sol = classify.read_examples(loadingFile)
    r = {}
    for (x,y) in new_sol:
        r1=[];
        if x.plate not in r: r[x.plate]={}
        if x.well not in r[x.plate]: r[x.plate][x.well]={}
        try:
            ybar, cost, M = classify_example(x, mesPoids)
        except:
            print "pbl de NaN"
            x.truth = "PROBLEME DE NaN"
            continue
        else:
            for k in range(len(ybar)):
                r1.extend(ybar[k]) 
            x.truth = r1; dou=[]
            M=M.transpose()
            for k in range(len(cost)):
                dou.append((cost[k], M[k]))
            r[x.plate][x.well][x.index]=dou
            
    return new_sol, r

def formulation(hypotheses, costs):
    res = [[], [], []]
    app = [[],[]]
    ind = 0

    for hyp, t in hypotheses:
        c = costs[ind][0]
        if hyp[0]==-1:
            app[0].append(hyp[1])
            app[1].append(c)
        else:
            try:
                h0 = list(hyp[0])
            except TypeError:
                h0 = [hyp[0]]
            try:
                h1 = list(hyp[1])
            except TypeError:
                h1 = [hyp[1]]
                
            res[0].append(h0)
            res[1].append(h1)
            res[2].append(c)
        ind+=1

    res[0].append([-1])
    res[1].append(app[0])
    res[2].append(app[1])
    return res

def calcAccumulation(new_sol, newFrameLot, r, cost=True):
    stories = tracking.ordonnancement(new_sol) 
    acc = np.zeros(shape=(len(EVENTS), 3), dtype=list)

    for plate in newFrameLot.lstFrames:
        for well in newFrameLot.lstFrames[plate]:
            for index in newFrameLot.lstFrames[plate][well]:
                print plate, well, index
                try:
                    costs = r[plate][well][index]
                except KeyError:
                    continue
                sourceC, targetC = accuracy.filtre(newFrameLot.lstFrames[plate][well][index])
#                for b in range(len(sourceC)):
#                    print "training set :", sourceC[b], targetC[b]
                if index not in stories[plate][well]:
                    continue
                solC = stories[plate][well][index]
                if type(solC.truth)==str:
                    print "attention, il n'y a pas de prediction pour cet index"
                    continue
                costs2 = filter(lambda x: x[1]==1, zip(costs, solC.truth))

                predicted =formulation(filter(lambda x: x[1]==1, zip(solC.hypotheses, solC.truth)), costs2)
                
#                for k in range(len(predicted[0])):
#                    print "predictions :", predicted[0][k], predicted[1][k]
                i=0; upletsVus = []
                for uplet in sourceC:
                    t, ttarg = accuracy.movement(uplet, s = None, tC=targetC[i])
                    for u in uplet:
                        #moi je voudrais aussi l'index dans predicted
                        ind, (p, ptarg) = movement(u, predicted)
                        if p==1:
                            #c'est une apparition
                            try:    
                                ind = predicted[1][-1].index(ttarg)
                            except ValueError:
                                continue
                            else:
                                upletsVus.append(ttarg)

                 #       pdb.set_trace()
                        if t==p:
                            test = (ttarg==ptarg) if t!=1 else (ttarg in ptarg)
                            if test and t!=1:
                                try:
                                    acc[t][0].append(predicted[2][ind])
                                except AttributeError:
                                    acc[t][0]=[predicted[2][ind]]
                            elif test and t==1:
 #                               pdb.set_trace()
                                try:
                                    acc[t][0].append(predicted[2][-1][ind])
                                except AttributeError:
                                    acc[t][0]=[predicted[2][-1][ind]]
                        else:
                            print '-------------------------------------------------------pas pareil...'
                            if p!=-1:
                                try:
                                    acc[p][1].append(predicted[2][ind])
                                except AttributeError:
                                    acc[p][1]=[predicted[2][ind]]
                            else:
                                try:
#                                    pdb.set_trace()
                                    acc[p][1].append(predicted[2][-1][ind])
                                except AttributeError:
                                    acc[p][1]=[predicted[2][-1][ind]]   
                    i+=1
                #pdb.set_trace()
                for uplet in filter(lambda x: x!=-1 and x not in sourceC, predicted[0]):
                    for u in uplet:
                        ind, (p, ptarg) = movement(u, predicted)
                        try:
                            acc[p][2].append(predicted[2][ind])
                        except AttributeError:
                            acc[p][2]=[predicted[2][ind]]
                #if upletsVus!=[]: pdb.set_trace()
                for u in filter(lambda x: x not in upletsVus, predicted[1][-1]):
                    p=1; ind = predicted[1][-1].index(u)
                    try:
                        acc[p][2].append(predicted[2][ind])
                    except AttributeError:
                        acc[p][2]=[predicted[2][ind]]

    return acc

def movement(uplet, s):
    z = s
    tu = filter(lambda x: uplet in x, s[0])[0]
    index = s[0].index(tu)
    return index, accuracy.movement(tu, s= None, tC= z[1][z[0].index(tu)])

def accumulation():
    global nb_folds; nb_fold = -1; nom='accumulation'; loadingFolder = '../results/'
    #i: charger les data
    #res = np.zeros(shape=(len(EVENTS), 3), dtype=list)
    
    listD = os.listdir('/media/lalil0u/New/workspace2/Tracking/data/raw')
    first = True
    for plate in listD:
        outputTotal = open(nom+plate+".pkl", 'w')
        sol, newFrameLot = accuracy.gettingSolu(first, plate)
        fichier = open(loadingFolder+"minMax_data_all.pkl", "r")  
        minMax = pickle.load(fichier)
        fichier.close()
        sol.normalisation(minMax)
        first = False
    
        #ii. predire et changer la forme de new_sol pour qu'elles soient exploitables
        new_sol, r = sousProcessClassify(sol, loadingFolder, loadingFile= None, i=None)
        
        #iii. construire le doublet source, target pour chaque frame, avec uniquement les hypotheses concernant les individus du training set
        acc = calcAccumulation(new_sol, newFrameLot, r)
        del sol ;del newFrameLot ; del new_sol ; gc.collect()
#        for j in range(res.shape[0]):
#            for k in range(res.shape[1]):
#                if type(acc[j][k])!=int:
#                    try:
#                        res[j][k].extend(acc[j][k])
#                    except AttributeError:
#                        res[j][k]=acc[j][k]
        
        pickle.dump(acc,outputTotal)   
#        s = np.sum(r, 1); d=np.zeros(r.shape)
#        for k in range(len(EVENTS)):
#            d[k]=r[k]/float(s[k])*100
#        print d
        #print>>outputTotal, d
        outputTotal.close()

if __name__ == '__main__':
    accumulation()
