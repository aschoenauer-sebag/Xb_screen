import numpy as np
import svmapi
import pdb
import cplexSolving



def read_examples(filename, sparm):
    """Reads and returns x,y example pairs from a file.
    
    This reads the examples contained at the file at path filename and
    returns them as a sequence.  Each element of the sequence should
    be an object 'e' where e[0] and e[1] is the pattern (x) and label
    (y) respectively.  Specifically, the intention is that the element
    be a two-element tuple containing an x-y pair."""
    global solutions
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
    sparm.c = 1 #?? TO SET
    sparm.epsilon = 0.01 #?? TO SET
    sm.size_psi = len(sample[0][0])
    print sm.size_psi
    
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
    number_events = len(y)
    sigma = 0
    prod = 0
    for i in range(number_events):
        sigma+=np.sum(y[i])
        prod +=np.sum( np.array(y[i])*np.array(ybar[i]))
    
    return (1 - prod/float(sigma)) 
    
#    VOIR LES CONTRAINTES SPECIALES ?
#def init_constraints(sample, sm, sparm):
#    """Initializes special constraints.
#
#    Returns a sequence of initial constraints.  Each constraint in the
#    returned sequence is itself a sequence with two items (the
#    intention is to be a tuple).  The first item of the tuple is a
#    document object.  The second item is a number, indicating that the
#    inner product of the feature vector of the document object with
#    the linear weights must be greater than or equal to the number
#    (or, in the nonlinear case, the evaluation of the kernel on the
#    feature vector with the current model must be greater).  This
#    initializes the optimization problem by allowing the introduction
#    of special constraints.  Typically no special constraints are
#    necessary.  A typical constraint may be to ensure that all feature
#    weights are positive.
#        POURQUOI ? pas les features weights ms les slack
#    Note that the slack id must be set.  The slack IDs 1 through
#    len(sample) (or just 1 in the combined constraint option) are used
#    by the training examples in the sample, so do not use these if you
#    do not intend to share slack with the constraints inferred from
#    the training data.
#        ???
#    The default behavior is equivalent to returning an empty list,
#    i.e., no constraints."""
#
#    if True:
#        # Just some example cosntraints.
#        c, d = svmapi.Sparse, svmapi.Document
#        # Return some really goofy constraints!  Normally, if the SVM
#        # is allowed to converge normally, the second and fourth
#        # features are 0 and -1 respectively for sufficiently high C.
#        # Let's make them be greater than 1 and 0.2 respectively!!
#        # Both forms of a feature vector (sparse and then full) are
#        # shown.
#        return [(d([c([(1,1)])],slackid=len(sample)+1),   1),
#                (d([c([0,0,0,1])],slackid=len(sample)+1),.2)]
#    # Encode positivity constraints.  Note that this constraint is
#    # satisfied subject to slack constraints.
#    constraints = []
#    for i in xrange(sm.size_psi):
#        # Create a sparse vector which selects out a single feature.
#        sparse = svmapi.Sparse([(i,1)])
#        # The left hand side of the inequality is a document.
#        lhs = svmapi.Document([sparse], costfactor=1, slackid=i+1+len(sample))
#        # Append the lhs and the rhs (in this case 0).
#        constraints.append((lhs, 0))
#    return constraints
    
def psi(x, y, sm, sparm):
    """Return a feature vector representing pattern x and label y.

    This is the combined feature function, which this returns either a
    svmapi.Sparse object, or sequence of svmapi.Sparse objects (useful
    during kernel evaluations, as all components undergo kernel
    evaluation separately).  There is no default behavior."""
    #x is a list of the features of all hypotheses of different events, in the order of EVENTS (importPack)
    #y is also a list, of the truth of all hypotheses of different events, same order
    #m in x is np.ndarray of length length_hyp(e), length_feat
    
#######si x ==sol.featuresMatrix()
#    number_events = len(x)
#    psi = np.array([])
#    pdb.set_trace()
#    for m in range(number_events):
#        fz = np.dot(np.transpose(x[m]), y[m])
#        psi = np.hstack((psi, fz))
#si x ==sol
    xTilde = x.featuresMatrix()
    number_events = len(xTilde)
    psi = np.array([])
    pdb.set_trace()
    for m in range(number_events):
        fz = np.dot(np.transpose(xTilde[m]), y[m])
        psi = np.hstack((psi, fz))

    return svmapi.Sparse(psi)    

########
#TO DO : 
#- implementer la fonction pr trouver la contrainte la plus violee
#=> pr cela il y a besoin d''implementer la fonction loss
#On fait appel au cplex pour trouver la solution
#
#apres il faudra bien voir tous les petits parametres caches et tester !

def find_most_violated_constraint_margin(x, y, sm, sparm):
    """Return ybar associated with x's most violated constraint.

    The find most violated constraint function for margin rescaling.
    The default behavior is that this returns the value from the
    general find_most_violated_constraint function."""
    #dans l'exemple donne il y a une iteration sur toutes les classifications possibles pour trouver la contrainte la plus violee
    #nous on ne peut pas faire ca donc on doit utiliser cplex
    number_events = len(y)
    length_feat= sm.size_psi
    sigma = 0
    Y = [] ; L = []
    
    #weights :
    w = sm.w
    
    longueurs = []
    for k in range(number_events):
        sigma+=np.sum(y[k])
        longueurs.append(len(y[k]))
        Y=np.hstack((Y, y[k]))
    
    Y=Y/sigma#type np.array, taille length_hyp,
    xTilde = x.featuresMatrix()
#SET UP COST FUNCTION
    for k in range(number_events):
        Mf = np.transpose(xTilde[k])
        print Mf.shape #type :np.ndarray, taille: length_feat, length_hyp(k)
        wC = w[length_feat*k:length_feat*(k+1)]
        print len(wC) #type : list longueur: length_feat
        L = np.hstack((L, np.dot(wC, Mf)))
    print L.shape#normalement c'est length_hyp
    costFunc = L-Y
    constraints = x.constraintsMatrix()
    print constraints.shape#(t1+t2, length_hyp)
    output = cplexSolving.solve(costFunc, constraints, x.singletsSize)#c'est ybar qu'il faut retourner
    ybar=[];last=0
    for k in range(number_events):
        ybar.append(list(output[last:last+longueurs[k]]))
        last +=longueurs[k]
    
    return ybar
    

####si x==sol.featuresMatrix()
#    for k in range(number_events):
#        Mf = np.transpose(x[k])
#        print Mf.shape #normalement : length_hyp(k), length_feat
#        wC = w[length_feat*k:length_feat*(k+1)]
#        print wC.shape #normalement length_feat
#        L = np.hstack((L, np.dot(wC, np.transpose(Mf))))
#    print L.shape#normalement c'est length_hyp
#    costFunc = L-Y
#    constraints = sol.constraintsMatrix()
####si x sol

