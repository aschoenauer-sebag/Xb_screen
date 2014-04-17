import os, pdb

def readPlateSetting(plateL, confDir, nb_row=None, nb_col=None, countEmpty=False, startAtZero = False):
    '''
    Function to go from csv file describing plate setup to a dictionary
    
    The csv file should not contain spaces in the cells
    '''
    result = {}
    for plate in plateL:
        result[plate]={}
        filename = os.path.join(confDir, "%s.csv" % plate)
        
        names=[]
        f=open(filename)
        lines=f.readlines(); f.close()
        
        lines=[line.strip("\n").split("\\") for line in lines]
        
        #loading experiment parameters
        if nb_col is None:
            try:
                a=int(lines[0][11])
            except:
                try:
                    a = int(lines[0][7])
                except:
                    raise AttributeError("Can't find the number of rows")
                else:
                    nb_col=8
            else:
                nb_col=12

            nb_row = len(lines)-1 #dire le nb de lignes dans le fichier - mais en fait on a toujours huit colonnes
            
        params=lines[0][nb_col:]
        k=0 if startAtZero else 1

        for line in lines[1:nb_row+1]:
            for el in line[:nb_col]:
                if el=='':
                    #if countEmpty on indique aussi que les puits vides sont vides. Sinon le numero des puits fait comme si les puits vides
                    #n'avaient jamais existe
                    
                    if countEmpty:
                        result[plate][k]={'Name': 'empty'}
                        k+=1
                    continue
                
                result[plate][k]={'Name': el}
                for i,param in enumerate(params):
                    result[plate][k][param]=line[nb_col+i]
                k+=1
                #name+='W{:<05}'.format(1+len(names))
    
    #wells=range(1,len(names)+1)
#    l=result[plate].keys(); l.sort()
#    for el in l:
#        print el, result[plate][el]
    return result, params