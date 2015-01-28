
import os, pdb
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as p

from operator import itemgetter

from util import settings
from collections import defaultdict
from itertools import product
from scipy.stats import spearmanr
import rpy2.robjects as objects
locfit = objects.packages.importr('locfit')
globalenv = objects.globalenv

from util.plots import couleurs, markers
from analyzer import CONTROLS, quality_control, plates

def parameterStabilityEvaluation(well_num=50,n_list=list([0.4,0.5,0.6,0.7,0.8]), h_list=[10], deg_list=[2,3],
                                 pheno_list=['Anaphase_ch1', 'Apoptosis_ch1', 'Folded_ch1', 'Interphase_ch1', 'Metaphase_ch1',
                                           'Polylobbed_ch1', 'Prometaphase_ch1', 'WMicronuclei_ch1'],
                               loadingFolder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results', passed=None):
    '''
    Here we're looking at how well results are robust to parameters change (result r1) and how well a subset of
    parameters is able to separate experiments (variance of r2).
    '''
    param_num=len(n_list)*len(h_list)*len(deg_list)
    param_list=[]
    if passed is None:
        passed,_,_,_=quality_control.computingToRedo(threshold_flou=0.4, threshold_init_cell_count=20)
    selected=np.array(passed)[np.random.permutation(len(passed))[:well_num]]
    processedDict={};
    result=defaultdict(list)
    for pl,w in selected:
        if pl not in processedDict:
            print "Opening processed dictionary file for pl ", pl
            f=open(os.path.join(loadingFolder, "processedDictResult_P{}.pkl".format(pl)))
            processedDict[pl]=pickle.load(f); f.close()
        try:
            currCtrl=quality_control.usable_XBSC(CONTROLS[processedDict[pl][int(w)]['Xenobiotic']], 0, pl)
        except KeyError:
            if processedDict[pl][int(w)]['Xenobiotic'] in CONTROLS.values():
                currCtrl=quality_control.usable_XBSC(processedDict[pl][int(w)]['Xenobiotic'], 0, pl)
            else:
                raise KeyError
        
        for pheno in pheno_list:
            y_exp=np.array(processedDict[pl][int(w)][pheno])
            for n,h,deg in product(n_list, h_list, deg_list):
                if (n,h,deg) not in param_list:
                    param_list.append((n,h,deg))
                    
                r=[]
                for pl, ctrl in filter(lambda x: x[1]!=w,currCtrl):
                    y_ctrl = np.array(processedDict[pl][int(ctrl)][pheno])
                    _, r_exp=localReg(y_exp, n, h, deg, plot=False)
                    _, r_ctrl = localReg(y_ctrl, n, h, deg, plot=False)
                    r.append(np.max(np.abs(r_exp-r_ctrl)))
                    
                result[pheno].append(np.median(r))
    print "Done with local regressions"
    for pheno in pheno_list:
        np.reshape(result[pheno], newshape=(well_num, param_num))
        
    r1=[np.ones(shape=(param_num, param_num), dtype=float) for pheno in pheno_list]
    r2=[np.zeros(shape=(param_num,), dtype=list) for pheno in pheno_list]
    for i, pheno in enumerate(pheno_list):
        result[pheno]=np.reshape(result[pheno], newshape=(well_num, param_num))
        for k in range(param_num):
            r2[i][k]=[0]
            iter_ = np.nditer(result[pheno][1:,k], flags=['c_index']) #one-dim array: the order is the same in C or F
            while not iter_.finished :
                r2[i][k]= np.hstack((r2[i][k],[a-iter_[0]  for a in np.nditer(result[pheno][:iter_.index+1,k])]))
                iter_.iternext()  
            r2[i][k]=np.abs(r2[i][k][1:])/np.mean(result[pheno][:,k])
            
            for j in range(k):
                r1[i][k,j]=spearmanr(result[pheno][:,k],result[pheno][:,j])[0]
                r1[i][j,k]=r1[i][k,j]
    return passed, param_list,r1, r2

def globalParameterComparison(well_num=50, perc=20, n_list=list([0.4,0.5,0.6,0.7,0.8]), h_list=[10], deg_list=[2,3],
                               pheno_list=['Anaphase_ch1', 'Apoptosis_ch1', 'Folded_ch1', 'Interphase_ch1', 'Metaphase_ch1',
                                           'Polylobbed_ch1', 'Prometaphase_ch1', 'WMicronuclei_ch1'],
                               loadingFolder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results', passed=None):
    '''
    The idea here is to study how well a local regression is generalizable by doing modeling without a subset of
    points and then looking at the errors of prediction on those points
    
    It might not be such a good idea since we have very noisy curves
     
    '''
    
    if passed is None:
        passed,_,_,_=quality_control.computingToRedo(threshold_flou=0.4, threshold_init_cell_count=20)
    selected=np.array(passed)[np.random.permutation(len(passed))[:well_num]]
    
    y_list=[]
    
    processedDict={}
    for pl,w in selected:
        if pl not in processedDict:
            print "Opening processed dictionary file for pl ", pl
            f=open(os.path.join(loadingFolder, "processedDictResult_P{}.pkl".format(pl)))
            processedDict[pl]=pickle.load(f); f.close()
            
        for pheno in pheno_list:
            y_list.append(np.array(processedDict[pl][int(w)][pheno]))
            
    r=parameterComparison(y_list, n_list, h_list, deg_list, perc, plot=True)
    return r
    

def selectLocalReg(y,n_list=list(np.linspace(0.5,0.8,4)), h_list=[2,10,20], deg_list=[2,3], plot=True, ax=None):
    if ax is None:
        f=p.figure(); ax=f.add_subplot(111)
    for n,h,deg in product(n_list, h_list, deg_list):
        result, prediction = localReg(y, n, h, deg, plot=plot, ax=ax, color=couleurs[n_list.index(n)], marker=markers[n_list.index(n)])
    if ax is None:
        p.show()
        
    return

def parameterComparison(y_list, n_list=list([0.4,0.5,0.6,0.7,0.8]), h_list=[10], deg_list=[2,3], perc=10, plot=True):
    '''
    This function compares different parameter sets: how well does the local regression performs on points that
    were left out during modeling?
    
    - y_list: list of different data sets. Important because we need the same parameter set for different types of curves
    
    '''
    r={}
    for n,h,deg in product(n_list, h_list, deg_list):
        r[(n,h,deg)]=crossValidation(y_list, n, h, deg, perc=perc)
        
    if plot:
        f=p.figure(); ax=f.add_subplot(111)
        w=[r[el] for el in sorted(r)]
        ax.boxplot(w)
        ax.set_xticklabels([el for el in sorted(r)])
        ax.tick_params(axis='both', which='major', labelsize=8)
        p.show()
    return r

def crossValidation(y_list,n,h,deg, perc=10):
    r=[]
    number_cv=100/perc

    for y in y_list:
        m=np.mean(y)
        time_length=y.shape[0]
        
        aside = int(time_length/perc)
        
        for k in range(number_cv):
            perm=np.random.permutation(time_length)
            test_=objects.FloatVector(perm[:aside]); train_set=y[perm[aside:]]
            model = localReg(train_set, n, h, deg, plot=False)
            r.append(np.sum(((np.array(locfit.predict_locfit(model, test_, what='coef'))- y[test_])/m)**2))
        
    return r
        

def localReg(y, n, h, deg=3, plot=True, ax=None, color=None, marker=None):
    '''
    Using locfit to perform local regression. Formula http://cm.bell-labs.com/stat/doc/locfitscg.ps
    - n: nearest neighbour fraction for the local least squares criterion
    - h: size of the bandwith
    - deg: degree of polynome used in local regression
    '''
    time_length=y.shape[0]
    
    x=objects.FloatVector(range(time_length))
    y=objects.FloatVector(y)
    data_frame= objects.DataFrame({'x': x, 'y': y})
    
    globalenv['x'] = x
    globalenv['y'] = y
    globalenv['df'] = data_frame
    
    result = locfit.locfit(objects.Formula("y~lp(x, nn={}, h={}, deg={})".format(n,h,deg)), data=data_frame)
    prediction=locfit.predict_locfit(result, objects.FloatVector(range(time_length)), what="coef")
    if plot and ax is not None:
        ax.scatter(range(time_length), y, color='green',s=7)
        ax.scatter(range(time_length), prediction, color=color, marker=marker, label='Param n {} h {} deg {}'.format(n,h,deg), s=3)
        ax.legend(fontsize=7); ax.grid(True)
    
    return result, np.array(prediction)

class wellPhenoAnalysis(object):
    def __init__(self, settings, pl, w, compound, nn,h,deg,pheno_list,verbose,localRegMeasure=True, settings_file=None):
        if settings_file is not None:
            self.settings = settings.Settings(settings_file, globals())
        else:
            self.settings=settings
        self.pheno_list = self.settings.pheno_list if pheno_list is None else pheno_list
        
        self.plate = pl
        self.well=int(w)
        self.nn = nn
        self.h=h;self.deg=deg;self.localRegMeasure=localRegMeasure
        self.compound=compound
        self.verbose=verbose
        
    def _getData(self):
        f=open(os.path.join(self.settings.loadingFolder, 'processedDictResult_P{}.pkl'.format(self.plate)))
        d=pickle.load(f); f.close()
        
        return d[self.well], d
    
    def _getCtrlData(self, filter_, total_data):
        try:
            ctrl = CONTROLS[self.compound]
        except KeyError:
            if filter_ is None:
                raise AttributeError
            else:
                ctrl=self.compound
                
        ctrl_wells = quality_control.usable_XBSC(ctrl, 0, self.plate)
        if filter_ is not None:
            if self.verbose:
                print "Avant filtre", ctrl_wells, filter_

            ctrl_wells=np.array(ctrl_wells)[filter_]
            if self.verbose:
                print "Apres filtre", ctrl_wells
            
        return [total_data[int(w)] for pl,w in ctrl_wells]  

    def _analyze(self, exp_data, ctrl_data):
        result=np.zeros(shape=(len(self.pheno_list)))
        for i,pheno in enumerate(self.pheno_list):
            r=[]
            for data in ctrl_data:
                exp_data[pheno]=np.array(exp_data[pheno])
                data[pheno]=np.array(data[pheno])
                if self.localRegMeasure:
    #I look at the max between local regression curves over time
                    _,r1=localReg(exp_data[pheno], n=self.nn, h=self.h, deg=self.deg, plot=False)
                    _,r2=localReg(data[pheno], n=self.nn, h=self.h, deg=self.deg, plot=False)
                    r.append(np.max(np.abs(r1-r2)))
                else:
    #My intuition is that this is not going to be a good measure for all curves because they're super noisy when the number of cells is low
    #Hence I want to see what it gives if I look at the area between the two curves
    #I normalize because I know that for 201214 they last 166 frames and not 192
                    try:
                        r1=np.abs(exp_data[pheno]-data[pheno])[:192]
                    except ValueError:
                        r1=np.abs(exp_data[pheno][:166]-data[pheno][:166])
                    r.append(np.sum(r1)/r1.shape[0])
                    
            result[i]=np.median(r)
            
        return result
            
    def __call__(self, filter_=None):
        #i. Get the data
        exp_data, total_data = self._getData()
        
        #ii. Get control data
        control_data = self._getCtrlData(filter_, total_data)
        
        #iii. Measure the distance
        result = self._analyze(exp_data, control_data)
        
        return result

class xbPhenoAnalysis(object):
    #Un peu con j'aurais pu la faire heriter de la precedente
    def __init__(self, settings_file, compound, dose=0, nn=None, localRegMeasure=True, pheno_list=None, verbose=0):
        self.settings = settings.Settings(settings_file, globals())

        self.nn=self.settings.nn if nn is None else nn
        self.localRegMeasure=localRegMeasure
            
        self.pheno_list = self.settings.pheno_list if pheno_list is None else pheno_list
            
        self.compound = compound
        self.dose= dose
        self.verbose = verbose
        
    def parameters(self,exp_list):
        r={
           'localReg':self.localRegMeasure,
           'nn':self.nn,
           'h':self.settings.h,
           'deg':self.settings.deg,
           'pheno_list':tuple(self.pheno_list),
           'wells':tuple(exp_list)
           }

        return tuple(sorted(zip(r.keys(), r.values()), key=itemgetter(0)))
        
    def _getExperiments(self):
        if self.compound not in CONTROLS.values():
            return quality_control.usable_XBSC(self.compound, self.dose), None
        else:
            exp_list = []
            filter_=[]
            for plate in plates:
                l= quality_control.usable_XBSC(self.compound, self.dose, plate=plate)
                size_=2 if len(l)>3 else 1
                f_=np.random.permutation(len(l))
                filter_.append(f_[size_:]); exp_list.extend(['{}_{}'.format(el[0], el[1]) for el in np.array(l)[f_[:size_]]])
                if self.verbose:
                    print "Chosen for plate ", plate, np.array(l)[f_[:size_]]
                
            return exp_list, filter_

    def save(self, result,exp_list):
        if self.settings.outputFile.format(self.compound, self.dose) in os.listdir(self.settings.savingFolder):
            f=open(os.path.join(self.settings.savingFolder, self.settings.outputFile.format(self.compound, self.dose)), 'r')
            d=pickle.load(f); f.close()
        else:
            d={}
            
        f=open(os.path.join(self.settings.savingFolder, self.settings.outputFile.format(self.compound, self.dose)), 'w')    
        d[self.parameters(exp_list)]=result
        pickle.dump(d,f);f.close()
        
        return


    def __call__(self):
        #i.Need to get the corresponding experiments - separate two cases depending if it's a control or no
        exp_list, filter_ = self._getExperiments()

        #ii. Need to measure the distances
        r=np.zeros(shape=(len(exp_list), len(self.pheno_list)))
        for i,exp in enumerate(exp_list):
            pl,w=exp if len(exp)==2 else exp.split('_')
            print pl
            wellPA = wellPhenoAnalysis(self.settings, pl, w, compound=self.compound, 
                                       h=self.settings.h, deg=self.settings.deg,localRegMeasure=self.localRegMeasure,
                                       nn=self.nn, pheno_list=self.pheno_list, verbose=self.verbose)
            r[i]=wellPA(filter_=filter_[plates.index(pl)]) if filter_ is not None else wellPA() 
            
        #iii. Need to save it
        self.save(r, exp_list)
        
        return 1
        