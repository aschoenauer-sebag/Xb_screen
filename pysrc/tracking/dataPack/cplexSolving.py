import cplex as c
import getpass

compteur=0

def solve(objective, constraints, t1):
    global compteur
    
    msg=""
    my_prob = c.Cplex()
    if getpass.getuser()!='lalil0u':
        #If i'm not here on my computer, I want that there is only one CPU that is used. 
        #This is to prevent Cplex from trying to use all CPUs available because it causes the cluster to crash when used on the cluster
        my_prob.parameters.threads.set(1)
    
    my_prob.set_log_stream(None); my_prob.set_results_stream(None)
    my_prob.objective.set_sense(my_prob.objective.sense.maximize)
#    print my_prob.get_problem_type() bien LP
#    my_prob.set_problem_type(my_prob.problem_type.LP)
    
    #set senses pour dire le sens des contraintes lin
    expected_size = objective.shape[0]
    #msg+="taille du vecteur output"+str(objective.shape)
    ub = [1 for x in range(expected_size)]

    #je passe la fonction objectif
#    print "upper bounds sur le vecteur output", ub
    my_prob.variables.add(obj = objective, ub = ub, types =[my_prob.variables.type.integer] * len(objective), names = ["z_"+str(i) for i in range(len(objective))])
    
    #je passe les contraintes
    #msg+="\n il y a "+str(constraints.shape[0])+" contraintes, et elles sont de taille "+str(constraints.shape[1])
    #msg+="\n taille de la contrainte"+str(constraints[0].shape)
    for i in range(constraints.shape[0]):
        bou = [[range(constraints.shape[1]),constraints[i]]]
        if i==0:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[0], names=["c"+str(i)])
        elif i==t1:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[0], names=["c"+str(i)])        
        else:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[1], names=["c"+str(i)])    
    #msg+= "\n nb de variables"+str( my_prob.linear_constraints.get_num())
    my_prob.solve()
    
    msg+= "\n Solution status = " +str(my_prob.solution.get_status())+":"
    # the following line prints the corresponding string
    msg+=str( my_prob.solution.status[my_prob.solution.get_status()])
    msg+= "\n Solution value  = "+str( my_prob.solution.get_objective_value())
    #slack = my_prob.solution.get_linear_slacks()
    #pi    = my_prob.solution.get_dual_values()
    x     = my_prob.solution.get_values()
    #dj    = my_prob.solution.get_reduced_costs()
    
    #msg+="\n Slack values :"+str(slack)+"\n"

#    for i in range(my_prob.linear_constraints.get_num()):
#        msg+="\n Row %d:  Slack = %10f " % (i, slack[i])
#    for j in range(my_prob.variables.get_num()):
#        msg+= "\n Column %d:  Value = %10f" % (j, x[j])
#
#    my_prob.write("7_26 index"+str(compteur)+".lp")
    compteur+=1
    return x, msg

def solveqp(objective, quadobj, constraints, t1):
    msg=""
    my_prob = c.Cplex()
    my_prob.set_log_stream(None)
    my_prob.objective.set_sense(my_prob.objective.sense.maximize)
#    print my_prob.get_problem_type() bien LP
#    my_prob.set_problem_type(my_prob.problem_type.LP)
    
    #set senses pour dire le sens des contraintes lin
    expected_size = objective.shape[0]
    #msg+="taille du vecteur output"+str(objective.shape)
    ub = [1 for x in range(expected_size)]

    #je passe la fonction objectif
#    print "upper bounds sur le vecteur output", ub
    my_prob.variables.add(obj = objective, ub = ub, types =[my_prob.variables.type.integer] * len(objective), names = ["z_"+str(i) for i in range(len(objective))])
    
    #je passe les contraintes
    #msg+="\n il y a "+str(constraints.shape[0])+" contraintes, et elles sont de taille "+str(constraints.shape[1])
    #msg+="\n taille de la contrainte"+str(constraints[0].shape)
    for i in range(constraints.shape[0]):
        bou = [[range(constraints.shape[1]),constraints[i]]]
        if i==0:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[0], names=["c"+str(i)])
        elif i==t1:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[0], names=["c"+str(i)])        
        else:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[1], names=["c"+str(i)])    
    #msg+= "\n nb de variables"+str( my_prob.linear_constraints.get_num())
    my_prob.solve()
    
    msg+= "\n Solution status = " +str(my_prob.solution.get_status())+":"
    # the following line prints the corresponding string
    msg+=str( my_prob.solution.status[my_prob.solution.get_status()])
    msg+= "\n Solution value  = "+str( my_prob.solution.get_objective_value())
    #slack = my_prob.solution.get_linear_slacks()
    #pi    = my_prob.solution.get_dual_values()
    x     = my_prob.solution.get_values()
    #dj    = my_prob.solution.get_reduced_costs()
    
    #msg+="\n Slack values :"+str(slack)+"\n"

#    for i in range(my_prob.linear_constraints.get_num()):
#        msg+="\n Row %d:  Slack = %10f " % (i, slack[i])
#    for j in range(my_prob.variables.get_num()):
#        msg+= "\n Column %d:  Value = %10f" % (j, x[j])

    my_prob.write("lpex1_Alice.lp")
    return x, msg

def solve_notInt(objective, constraints, t1):
    msg=""
    my_prob = c.Cplex()
    my_prob.objective.set_sense(my_prob.objective.sense.maximize)
#    print my_prob.get_problem_type() bien LP
#    my_prob.set_problem_type(my_prob.problem_type.LP)
    
    #set senses pour dire le sens des contraintes lin
    expected_size = objective.shape[0]
    #msg+="taille du vecteur output"+str(objective.shape)
    ub = [1 for x in range(expected_size)]

    #je passe la fonction objectif
#    print "upper bounds sur le vecteur output", ub
    my_prob.variables.add(obj = objective, ub = ub, names = ["z_"+str(i) for i in range(len(objective))])
    
    #je passe les contraintes
    #msg+="\n il y a "+str(constraints.shape[0])+" contraintes, et elles sont de taille "+str(constraints.shape[1])
    #msg+="\n taille de la contrainte"+str(constraints[0].shape)
    for i in range(constraints.shape[0]):
        bou = [[range(constraints.shape[1]),constraints[i]]]
        if i==0:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[0], names=["c"+str(i)])
        elif i==t1:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[0], names=["c"+str(i)])        
        else:
            my_prob.linear_constraints.add(lin_expr=bou, senses='E', rhs=[1], names=["c"+str(i)])    
    #msg+= "\n nb de variables"+str( my_prob.linear_constraints.get_num())
    my_prob.solve()
    
    msg+= "\n Solution status = " +str(my_prob.solution.get_status())+":"
    # the following line prints the corresponding string
    msg+=str( my_prob.solution.status[my_prob.solution.get_status()])
    msg+= "\n Solution value  = "+str( my_prob.solution.get_objective_value())
    #slack = my_prob.solution.get_linear_slacks()
    #pi    = my_prob.solution.get_dual_values()
    x     = my_prob.solution.get_values()
    #dj    = my_prob.solution.get_reduced_costs()
    
    #msg+="\n Slack values :"+str(slack)+"\n"

#    for i in range(my_prob.linear_constraints.get_num()):
#        msg+="\n Row %d:  Slack = %10f " % (i, slack[i])
#    for j in range(my_prob.variables.get_num()):
#        msg+= "\n Column %d:  Value = %10f" % (j, x[j])

#    my_prob.write("lpex1.lp")
    return x, msg

#a la fin un des conseils est d'enregistrer le pbl ds un fichier et d'essayer de le resoudre dans l'interactive optimizer de ibm