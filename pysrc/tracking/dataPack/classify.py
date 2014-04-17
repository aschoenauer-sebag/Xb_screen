import cPickle as pickle
import numpy as np
import cplexSolving
import pdb

def read_examples(filename):
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

def loss(y, ybar):
    #LOSS par rapport a chaque evenement
    number_eventsC = len(y)
    #sigma = 0
    #prod = 0
    result = [0 for x in range(number_eventsC)]
    for i in range(number_eventsC):
        #sigma=np.sum(y[i])
        result[i]=np.sum( np.array(y[i])*np.array(ybar[i]))
        #donc chaque loss est une liste de cinq chiffres
    return np.array(result) 

def contenu(y):
    #LOSS par rapport a chaque evenement
    number_eventsC = len(y)

    result = [np.sum(y[i]) for i in range(number_eventsC)]
    #donc chaque loss est une liste de cinq chiffres
    return np.array(result) 
#
#def create_conf_matrix(expected, predicted, n_classes):
#    m = [[0] * n_classes for i in range(n_classes)]
#    for pred, exp in zip(predicted, expected):
#        m[pred][exp] += 1
#    return m

def calc_accuracy(total_loss, test_contenu):
#    loss = loss(y, ybar); contenu = contenu(y)
    #pdb.set_trace()
    number_events = total_loss.shape[1]
    r1 = [0 for x in range(number_events)]
    a = np.sum(total_loss, 0) ; b = np.sum(test_contenu, 0)
    for i in range(number_events):
        if b[i]>0:
            r1[i]= a[i]/float(b[i])*100
        else:
            r1[i]=0
        
    r2 = np.sum(total_loss)/float(np.sum(test_contenu))*100
    return r1, r2

def classify_example(x, weights):
    """Given a pattern x, return the predicted label."""

    number_events = len(x.events)
    length_feat= x.length_feat()
    L = []
    
    #weights :
    w = weights

    longueurs = []
    xTilde = x.featuresMatrix()
#SET UP LIKELIHOOD FUNCTION
    for k in range(number_events):
        Mf = np.transpose(xTilde[k])
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
    
    print "--",
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
    return ybar

def classify_example_notInt(x, weights):
    """Given a pattern x, return the predicted label."""

    number_events = len(x.events)
    length_feat= x.length_feat()
    L = []
    
    #weights :
    w = weights
    #pdb.set_trace()
    print np.sum(w)
    longueurs = []
    xTilde = x.featuresMatrix()
#SET UP LIKELIHOOD FUNCTION
    for k in range(number_events):
        Mf = np.transpose(xTilde[k])
        print Mf.shape #type :np.ndarray, taille: length_feat, length_hyp(k)
        
        if Mf.shape==(0,):
            longueurs.append(0)
            continue
        longueurs.append(Mf.shape[1])
        wC = w[length_feat*k:length_feat*(k+1)]
        print len(wC) #type : list longueur: length_feat
        L = np.hstack((L, np.dot(wC, Mf)))
    print L.shape#normalement c'est length_hyp
    costFunc = L
    constraints = x.constraintsMatrix()
    print constraints.shape#(t1+t2, length_hyp)
    
#FIND ARGMAX LIKELIHOOD FUNCTION
    output, msg = cplexSolving.solve_notInt(costFunc, constraints, x.singletsSize)#c'est ybar qu'il faut retourner
    
    print "------------------------------------------####################MSG CPLEX"
    print msg
#    output = [int(j) for j in output2]
    ybar=[];last=0
    for k in range(number_events):
        ybar.append(list(output[last:last+longueurs[k]]))
        last +=longueurs[k]

    print "ybar : "
    for k in range(number_events):
        print ybar[k]
#    print "loss :", loss(y, ybar, sparm)
    return ybar
    
def read_model(filename, sparm):
    """Load the structure model from a file.
    
    Return the structmodel stored in the file at path filename, or
    None if the file could not be read for some reason.

    The default behavior is equivalent to
    'return cPickle.load(bz2.BZ2File(filename))'."""
#    sm= pickle.load(bz2.BZ2File(filename, 'r'))
    f = open(filename, 'r')
    w= pickle.load(f)
    f.close()
    #pdb.set_trace()
#    for i in sm.w: print i
    return w

def write_label(fileptr, y):#
    """Write a predicted label to an open file.

    Called during classification, this function is called for every
    example in the input test file.  In the default behavior, the
    label is written to the already open fileptr.  (Note that this
    object is a file, not a string.  Attempts to close the file are
    ignored.)  The default behavior is equivalent to
    'print>>fileptr,y'"""
    print>>fileptr,y
    
def eval_prediction(exnum, (x, y), ypred, sm, sparm, teststats):
    """Accumulate statistics about a single training example.
    
    Allows accumulated statistics regarding how well the predicted
    label ypred for pattern x matches the true label y.  The first
    time this function is called teststats is None.  This function's
    return value will be passed along to the next call to
    eval_prediction.  After all test predictions are made, the last
    value returned will be passed along to print_testing_stats.

    On the first call, that is, when exnum==0, teststats==None.  The
    default behavior is that the function does nothing."""
    if exnum==0: teststats = []
    print 'on example',exnum,'predicted',ypred,'where correct is',y
    teststats.append(loss(y, ypred, sparm))
    return teststats

def print_testing_stats(sample, sm, sparm, teststats):
    """Print statistics once classification has finished.
    
    This is called after all test predictions are made to allow the
    display of any summary statistics that have been accumulated in
    the teststats object through use of the eval_prediction function.

    The default behavior is that nothing is printed."""
    print teststats