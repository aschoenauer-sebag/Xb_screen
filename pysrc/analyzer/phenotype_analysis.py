
import os, pdb, getpass
import numpy as np
import cPickle as pickle
from optparse import OptionParser

from operator import itemgetter

from util import settings
from collections import defaultdict, Counter
from itertools import product
from scipy.stats import spearmanr
import rpy2.robjects as objects
if getpass.getuser()=='aschoenauer':
    locfit = objects.packages.importr('locfit', lib_loc="/cbio/donnees/aschoenauer/software/lib64/R")

    import matplotlib as mpl
    mpl.use('Agg')
    
elif getpass.getuser()=='lalil0u':
    locfit = objects.packages.importr('locfit')
    
import matplotlib.pyplot as p
globalenv = objects.globalenv

from tracking.histograms.summaries_script import progFolder, scriptFolder, pbsArrayEnvVar, pbsErrDir, pbsOutDir, path_command

from util.plots import couleurs, markers, plotBarPlot
from analyzer import CONTROLS, quality_control, plates, xbL, compoundL

def plotResults(result, who, dose_list, compound_list, outputFile, outputFolder, features=False):
    compounds = sorted(filter(lambda x: x in xbL or x=='Rien', Counter(compound_list).keys()))
    f,axes=p.subplots(len(compounds), len(result), sharey=True)
    for i,compound in enumerate(compounds):
        control=CONTROLS[compound]
        where_c = np.where(compound_list==control)
        print control, where_c
        ww_=[np.where((compound_list==compound)&(dose_list==k))[0] for k in range(11)]
        for j,pheno in enumerate(result):
            if not features:
                data=[result[pheno][where_c]]
                data.extend([result[pheno][el] for el in filter(lambda x: len(x)>0, ww_)]) 
            else:
                data=[np.hstack((el for el in result[pheno][where_c]))]
                for el in filter(lambda x: len(x)>0, ww_):
                    data.append(np.hstack((result[pheno][k] for k in el)) )
                
            if compound !='Rien':
                labels=['Ctrl']; labels.extend([dose for dose in sorted(Counter(dose_list[np.where(compound_list==compound)]).keys())])
            else:
                if not features:
                    data.append(result[pheno][np.where(compound_list=='DMSO')]);data.append(result[pheno][np.where(compound_list=='TGF')])
                else:
                    data.append(np.hstack((el for el in result[pheno][np.where(compound_list=='DMSO')])))
                    data.append(np.hstack((el for el in result[pheno][np.where(compound_list=='TGF')])))
                labels=['Nonane', 'Rien', 'DMSO', 'TGF']
            axes[i,j].boxplot(data); axes[i,j].set_title('{}'.format(pheno), fontsize=10)
            axes[i,j].set_ylabel(compound)
            axes[i,j].tick_params(axis='both', which='major', labelsize=5)
            axes[i,j].set_xticklabels(labels)
    p.show()
    
def plotResults2(result, who, dose_list, compound_list, outputFile, outputFolder, features=False):
    compounds = sorted(filter(lambda x: x in xbL or x=='Rien', Counter(compound_list).keys()))
    f,axes=p.subplots(len(compounds), len(result), sharey=True)
    for i,compound in enumerate(compounds):
        control=CONTROLS[compound]
        where_c = np.where(compound_list==control)
        ww_=[np.where((compound_list==compound)&(dose_list==k))[0] for k in range(11)]
        for j,pheno in enumerate(result):
            if not features:
                data=[result[pheno][where_c]]
                data.extend([result[pheno][el] for el in filter(lambda x: len(x)>0, ww_)])
            else:
                data=[np.hstack((el for el in result[pheno][where_c]))]
                for el in filter(lambda x: len(x)>0, ww_):
                    data.append(np.hstack((result[pheno][k] for k in el)) )
            
            if compound !='Rien':
                labels=['Ctrl']; labels.extend([dose for dose in sorted(Counter(dose_list[np.where(compound_list==compound)]).keys())])
            else:
                if  not features:
                    data.append(result[pheno][np.where(compound_list=='DMSO')]);data.append(result[pheno][np.where(compound_list=='TGF')])
                else:
                    data.append(np.hstack((el for el in result[pheno][np.where(compound_list=='DMSO')])))
                    data.append(np.hstack((el for el in result[pheno][np.where(compound_list=='TGF')])))
                labels=['Nonane', 'Rien', 'DMSO', 'TGF']
                
            plotBarPlot(np.array([[np.mean(el) for el in data]]), pheno, compound, labels, title=compound,name=None,target=None,\
                        stds=[np.std(el) for el in data],fig=f,ax=axes[i,j], sh=False, save=False )
#             axes[i,j].boxplot(data); axes[i,j].set_title('{} {}'.format(compound[:4], pheno.split('_')[0]), fontsize=10)
#             axes[i,j].tick_params(axis='both', which='major', labelsize=5)
#             axes[i,j].set_xticklabels(labels)
    p.show()
    return
    
def plotResults3(result, who, dose_list, compound_list, outputFile, outputFolder):
    f=p.figure()
    for i, pheno in enumerate(result):
        ax=f.add_subplot(1,8,i)
        ax.set_title(pheno)
        a=zip(result[pheno],compound_list, dose_list)
        a=sorted(a, key=itemgetter(0))
        for k,el in enumerate(a[250:]):
            ax.scatter(k, el[0], color=couleurs[compoundL.index(el[1])])
            ax.text(k, el[0], "{}{}".format(el[1][:2], el[2]), fontsize=10)
        ax.grid(True)
    p.show()
    return

def saveTextFile(localRegMeasure, result, who, dose_list, compound_list, outputFile, outputFolder):
    for pheno in result:
        f=open(os.path.join(outputFolder, 'ranking_pval_{}_localReg{}.txt'.format(pheno, localRegMeasure)), 'w')
        f.write('Distance\t Xb\t Dose\t Plate\t Well\n')
        a=zip(result[pheno],compound_list, dose_list, who)
        a=sorted(a, key=itemgetter(0), reverse=True)
        for el in a:
            f.write("{}\t {} \t {}\t {}\t {} \n".format(el[0], el[1], el[2], *el[3]))
        f.close()
    return

def loadResults(localRegMeasure=False, 
                loadingFolder='/media/lalil0u/New/projects/Xb_screen/dry_lab_results/MITOSIS/phenotype_analysis_up_down', 
                filename='phenoAnalysis_plateNorm_', 
                pheno_list = ['Anaphase_ch1', 'Apoptosis_ch1', 'Folded_ch1', 'Interphase_ch1', 'Metaphase_ch1',
                              'Polylobbed_ch1', 'Prometaphase_ch1', 'WMicronuclei_ch1', 'Frozen_ch1'],
                plot1=False, plot2=False, plot3=False, saveText=False):
    
    file_list=sorted(filter(lambda x: filename in x, os.listdir(loadingFolder)))
    result={pheno:[] for pheno in pheno_list}
    who=[]#for well list
    compound_list=[]; dose_list=[]
    for file_ in file_list:
        #try:
        f=open(os.path.join(loadingFolder, file_))
        d=pickle.load(f);f.close()
            
#         except IOError:
#             pdb.set_trace()
        #else:
        try:
            param_sets = [(i,dict(el)) for i,el in enumerate(d.keys())]
            index, param_d=filter(lambda x: x[1]['localReg']==localRegMeasure, param_sets)[0]
            wells=param_d["wells"] if len(param_d['wells'][0])==2 else [el.split('_') for el in param_d['wells']]
            available_phenos = param_d['pheno_list']
            arr=d[d.keys()[index]]
        except IndexError, KeyError:
            pdb.set_trace()
        else:
            for pheno in pheno_list:
                result[pheno] = arr[:,available_phenos.index(pheno)] if result[pheno]==[] else np.hstack((result[pheno], arr[:,available_phenos.index(pheno)]))
            who.extend(wells)
            compound_list.extend([file_.split('_')[-2] for k in range(len(wells))])
            dose_list.extend([file_.split('_')[-1].split('.')[0] for k in range(len(wells))])
    who=np.array(who); dose_list=np.array(dose_list, dtype=int); compound_list=np.array(compound_list) 
    if plot1:
        plotResults(result, who, dose_list, compound_list, outputFile=filename, outputFolder=loadingFolder)
    elif plot2:
        plotResults2(result, who, dose_list, compound_list, outputFile=filename, outputFolder=loadingFolder)
    elif plot3:
        plotResults3(result, who, dose_list, compound_list, outputFile=filename, outputFolder=loadingFolder)
    if saveText:
        saveTextFile(localRegMeasure, result, who, dose_list, compound_list, outputFile=filename, outputFolder=loadingFolder)
    return result, who, compound_list, dose_list
            
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
        ax.scatter(range(time_length), y, color=color,s=7)
        ax.scatter(range(time_length), prediction, color=color, marker=marker, label='Param n {} h {} deg {}'.format(n,h,deg), s=3)
        ax.legend(fontsize=7); ax.grid(True)
    
    return result, np.array(prediction)

class wellPhenoAnalysis(object):
    def __init__(self, settings_, pl, w, compound, nn,pheno_list,verbose,localRegMeasure=True, settings_file=None):
        if settings_ is None:
            self.settings = settings.Settings(settings_file, globals())
        else:
            self.settings=settings_
        self.pheno_list = self.settings.pheno_list if pheno_list is None else pheno_list
        
        self.plate = pl
        self.well=int(w)
        self.nn = nn
        
        self.h=self.settings.h
        self.deg=self.settings.deg
        self.localRegMeasure=localRegMeasure
        
        self.compound=compound
        self.verbose=verbose
        
    def _getData(self):
        f=open(os.path.join(self.settings.loadingFolder, 'processedDictResult_P{}.pkl'.format(self.plate)))
        d=pickle.load(f); f.close()
        
        return d[self.well], d
    
    def _getCtrlData_ctrlNorm(self, filter_, total_data):
        try:
            ctrl = CONTROLS[self.compound]
        except KeyError:
            if filter_ is None:
                raise AttributeError
            else:
                ctrl=self.compound
                
        ctrl_wells = self._getUsable(comp=ctrl, dose=0)
        if filter_ is not None:
            if self.verbose:
                print "Avant filtre", ctrl_wells, filter_

            ctrl_wells=np.array(ctrl_wells)[filter_]
            if self.verbose:
                print "Apres filtre", ctrl_wells
            
        return [total_data[int(w)] for pl,w in ctrl_wells]
    
    def _getCtrlData_plateNorm(self, filter_, total_data):
        '''
        We compute the percentage of phenotypes with respect to the number of cells that are not artefacts or out of focus cells.
        '''
        #Loading well list, that passed the QC - instead of recomputing it each time
        f=open(self.settings.ok_wells_asLIST, 'r')
        a=pickle.load(f); f.close()
        
        exp_list= [(pl, "{:>05}".format(w)) for pl, w in a[np.where(a[:,0]==self.plate)]]
        result=[]
    #i. Doing five permutations
        for k in range(5):
            currExp = np.array(exp_list)[np.random.permutation(range(len(exp_list)))[:self.settings.n_n]]
            r=defaultdict(list)
            for pheno in self.pheno_list:
                global_cell_count=None
                for pl,w in currExp:
     #ii. Going from counts+percentage/wells to percentages/well list
                    curr_object_count = np.array(total_data[int(w)]['object_count'], dtype=float)
                    curr_cell_count = np.array(total_data[int(w)]['cell_count'], dtype=float)
                    global_cell_count = curr_cell_count if global_cell_count is None \
                        else global_cell_count[:min(curr_cell_count.shape[0], global_cell_count.shape[0])]+curr_cell_count[:min(curr_cell_count.shape[0], global_cell_count.shape[0])]
                    
                    curr_pheno_count=np.array(total_data[int(w)][pheno])[:,0]*curr_object_count
                    if len(r[pheno])==0:
                        r[pheno]=curr_pheno_count
                    else:
                        r[pheno]= r[pheno][:min(r[pheno].shape[0], curr_pheno_count.shape[0])]
                        curr_pheno_count=curr_pheno_count[:min(r[pheno].shape[0], curr_pheno_count.shape[0])]
                        r[pheno]+=curr_pheno_count
                        
                r[pheno]/=global_cell_count[:,0]
            result.append(r)
        return result
        
    def _getCtrlData(self, filter_, total_data):
        '''
        Returning a list of control data arrays. The final result for the well will be the median of the results on each array of the list
        '''
        if self.settings.norm=='ctrl_neg':
            return self._getCtrlData_ctrlNorm(filter_, total_data)
        elif self.settings.norm=='plate':
            return self._getCtrlData_plateNorm(filter_, total_data)
        
    def _getUsable(self, comp=None, dose=None):
        try:
            f=open(self.settings.ok_wells_asDICT, 'r')
            ok_dict=pickle.load(f);f.close()
        except:
            if comp is None:
                return quality_control.usable_XBSC(self.compound, self.dose,confDir=self.settings.plate_setups_folder)
            else:
                return quality_control.usable_XBSC(comp, dose,confDir=self.settings.plate_setups_folder)
        else:
            if comp is None:
                return ok_dict[self.compound][self.dose]
            else:
                return ok_dict[comp][dose]
        
    def _analyze(self, exp_data, ctrl_data):
        result=np.zeros(shape=(len(self.pheno_list)))
        for i,pheno in enumerate(self.pheno_list):
            r=[]
            
            for j,data in enumerate(ctrl_data):
                local_exp_data=np.array(exp_data[pheno], dtype=float)[:,0]*np.array(exp_data['object_count'], dtype=float)/(np.array(exp_data['cell_count'], dtype=float)[:,0])
                data[pheno]=np.array(data[pheno])
                if self.localRegMeasure:
                    if j==0:
                        f=p.figure(figsize=(12,12)); ax=f.add_subplot(111)
                        plot=True
    #I look at the max between local regression curves over time, but stop at 48h not too add the bias of different durations
                    _,r1=localReg(local_exp_data[:192], n=self.nn, h=self.h, deg=self.deg, plot=plot, ax=ax, color=self.settings.COLORD[pheno], marker='o')
                    r1[np.where(r1<0)]=0
                    plot=False
                    
                    _,r2=localReg(data[pheno][:192], n=self.nn, h=self.h, deg=self.deg, plot=True, ax=ax, color=self.settings.plate_colors[self.plate], marker='*')
                    r2[np.where(r2<0)]=0
                    try:
                        arr = r1-r2
                    except ValueError:
                        arr = r1[:min(r1.shape[0], r2.shape[0])]-r2[:min(r1.shape[0], r2.shape[0])]
                    ind= np.argmax(np.absolute(arr))
                    r.append(arr[ind])

                else:
    #My intuition is that this is not going to be a good measure for all curves because they're super noisy when the number of cells is low
    #Hence I want to see what it gives if I look at the area between the two curves
    #I normalize because I know that for 201214 they last 166 frames and not 192
                    try:
                        r1=(exp_data[pheno]-data[pheno])[:192]
                    except ValueError:
                        r1=(exp_data[pheno][:min(exp_data[pheno].shape[0], data[pheno].shape[0])]-data[pheno][:min(exp_data[pheno].shape[0], data[pheno].shape[0])])
                    r.append(np.sum(r1)/r1.shape[0])
                    if np.isnan(r[-1]):
                        pdb.set_trace()
            if self.localRegMeasure:
                ax.set_ylim(self.settings.plot_ylim[i])
                ax.set_xlim(0,192)
                p.savefig(os.path.join(self.settings.savingFolder, self.settings.imageName.format(pheno, self.plate, self.well)))
                    
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
    
    #Rq: concretement les deux mesures se valent (correlation de Spearman 0.913)
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
    
    def _getUsable(self):
        try:
            f=open(self.settings.ok_wells_asDICT, 'r')
            ok_dict=pickle.load(f);f.close()
        except:
            return quality_control.usable_XBSC(self.compound, self.dose,confDir=self.settings.plate_setups_folder)
        else:
            return ok_dict[self.compound][self.dose]
        
    def _getExperiments(self):
        if self.settings.norm=='neg_ctrl':
            if self.compound not in CONTROLS.values():
                return self._getUsable(), None
            else:
                exp_list = []
                filter_=[]
                for plate in plates:
                    l= self._getUsable()
                    size_=2 if len(l)>3 else 1
                    f_=np.random.permutation(len(l))
                    filter_.append(f_[size_:]); exp_list.extend(['{}_{}'.format(el[0], el[1]) for el in np.array(l)[f_[:size_]]])
                    if self.verbose:
                        print "Chosen for plate ", plate, np.array(l)[f_[:size_]]
                    
                return exp_list, filter_
        elif self.settings.norm=='plate':
            return self._getUsable(), None
        else:
            raise ValueError

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
        if not os.path.isdir(self.settings.savingFolder):
            os.mkdir(self.settings.savingFolder)
        #i.Need to get the corresponding experiments - separate two cases depending if it's a control or no
        exp_list, filter_ = self._getExperiments()

        #ii. Need to measure the distances
        r=np.zeros(shape=(len(exp_list), len(self.pheno_list)))
        for i,exp in enumerate(exp_list):
            pl,w=exp if len(exp)==2 else exp.split('_')
            print pl
            wellPA = wellPhenoAnalysis(self.settings, pl, w, compound=self.compound, 
                                       localRegMeasure=self.localRegMeasure,
                                       nn=self.nn, pheno_list=self.pheno_list, verbose=self.verbose)
            r[i]=wellPA(filter_=filter_[plates.index(pl)]) if filter_ is not None else wellPA() 
            
        #iii. Need to save it
        self.save(r, exp_list)
        
        return 1
    
    
def script(condition_list_file='../data/siRNA_xb.pkl', baseName="phenotype_analysis"):
    f=open(condition_list_file)
    l=pickle.load(f)
    f.close()
    
    l.extend(['Nonane_0', 'DMSO_0'])
    
    head = """#!/bin/sh
cd %s""" %progFolder
    
    for i,condition in enumerate(l):
        script_name = os.path.join(scriptFolder, baseName+'{}.sh'.format(i+1))
        script_file = file(script_name, "w")
        cmd="\npython analyzer/phenotype_analysis.py -x {} -d {}\n".format(*condition.split('_'))
                
        script_file.write(head + cmd)
        script_file.close()
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
    sub_cmd = 'qsub -hold_jid  -t 1-%i %s' % (i+1, array_script_name)

    print sub_cmd
    
        
if __name__ == '__main__':
    verbose=0
    description =\
'''
%prog - 

'''
    parser = OptionParser(usage="usage: %prog [options]",
                         description=description)
    
    parser.add_option("-f", "--settings_file", dest="settings_file", default='analyzer/settings/settings_phenoAnalysis_thalassa.py',
                      help="Settings_file")

    parser.add_option("-x", dest="xenobiotic",type=str,
                      help="The xenobiotic which you are interested in")
    
    parser.add_option("-d", dest="dose", type=int,
                      help="The dose which you are interested in")

    parser.add_option("-n", dest="nn", default=0.6, type=float,
                      help="Number of neighbours for local regression")
    parser.add_option("-v", dest="verbose", default=0, type=int,
                      help="Verbosity level")
        
    (options, args) = parser.parse_args()
    
    x=xbPhenoAnalysis(options.settings_file, compound=options.xenobiotic, dose=options.dose, nn=options.nn, verbose=options.verbose, localRegMeasure=True)    
    x()
    x=xbPhenoAnalysis(options.settings_file, compound=options.xenobiotic, dose=options.dose, nn=options.nn, verbose=options.verbose, localRegMeasure=False)
    x()
