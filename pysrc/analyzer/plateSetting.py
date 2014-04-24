import os, pdb, datetime
from warnings import warn

from django.core.exceptions import ObjectDoesNotExist
from plates.models import Plate,Cond, Treatment, Well
from analyzer import *

def readPlateSetting(plateL, confDir, nb_row=None, nb_col=None, countEmpty=False, startAtZero = False,
                     plateName=None, dateFormat='%d%m%y', defaultMedium = "Complet",
                     addPlateWellsToDB=False, addCondToDB = False, addTreatmentToDB=False):
    '''
    Function to go from csv file describing plate setup to a dictionary
    
    The csv file should not contain spaces in the cells
    '''
    result = {}
    idL = {}
    for plate in plateL:
        result[plate]={}         
        idL[plate]={}
        filename = os.path.join(confDir, "%s.csv" % plate)
        
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
        
        if addPlateWellsToDB:
            p = Plate(name = plateName, date = datetime.datetime.strptime(plate, dateFormat), nb_col=nb_col, nb_row=nb_row)
            p.save()
        
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
                    
                    
                if addCondToDB:
                    try:
                        medium = result[plate][k]['Medium']
                    except KeyError:
                        medium=defaultMedium
                    try:
                        serum = int(result[plate][k]['Serum'])
                    except KeyError:
                        serum = None
                    
                    conds = Cond.objects.filter(medium = medium)
                    if len(conds)==0:
                        cond = addConds(medium, serum) 
                    elif len(conds)==1:
                        if conds[0].serum!=serum:
                            cond =addConds(medium, serum)
                        else:
                            cond = conds[0]
                    else:
                        conds = conds.filter(serum =serum)
                        if len(conds)==0:
                            cond =addConds(medium, serum)
                        elif len(conds)>1:
                            raise
                        else:
                            cond = conds[0]

                if addTreatmentToDB and 'Xenobiotic' in params:
                    xb = result[plate][k]['Xenobiotic']
                    dose = int(result[plate][k]['Doses'])
                    
                    treatments = Treatment.objects.filter(xb = xb)
                    if len(treatments)==0:
                        treatment = addTreatment(xb, dose)
                    elif len(treatments)==1:
                        if treatments[0].dose != dose:
                            treatment = addTreatment(xb, dose)
                        else:
                            treatment =treatments[0]
                    else:
                        treatments = treatments.filter(dose =dose)
                        if len(treatments)==0:
                            treatment = addTreatment(xb, dose)
                        elif len(treatments)>1:
                            raise
                        else:
                            treatment =treatments[0]
                        
                if addPlateWellsToDB:
                    try:
                        clone = result[plate][k]['Name'].split("_")[0]
                    except:
                        clone = None
                    w=Well(num=k, plate=p, treatment = treatment, clone = clone)
                    w.save()
            #this is how to add many to many relations
                    w.cond.add(cond); w.save()
                    idL[plate][k]=w.id
                k+=1
    if addPlateWellsToDB:
        return result, params, idL
    else:
        return result, params

def addConds(medium, serum):
    cond = Cond(medium = medium, serum = serum)
    cond.save()
    return cond

def addTreatment(xb, dose):
    is_ctrl = xb in CONTROLS.values()
    c=None
    if not is_ctrl:
        if xb not in CONTROLS:
            warn('This xenobiotic {} has no control defined yet.'.format(xb))
        else:
            controls = Treatment.objects.filter(xb = CONTROLS[xb])#.get(dose =dose)
            if len(controls)==0:
                c=addTreatment(CONTROLS[xb], dose)
            elif len(controls)==1:
                if controls[0].dose==dose:
                    c=controls[0]
                else:
                    c=addTreatment(CONTROLS[xb], dose)
            else:
                controls = controls.filter(dose=dose)
                if len(controls)==0:
                    c=addTreatment(CONTROLS[xb], dose)
                elif len(controls)>1:
                    raise
                else:
                    c=controls[0]
        
    treatment = Treatment(xb=xb, dose = dose, is_ctrl = is_ctrl, ctrl = c)
    treatment.save()
    return treatment


