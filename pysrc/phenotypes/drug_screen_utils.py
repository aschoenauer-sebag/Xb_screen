import matplotlib.pyplot as p
import cPickle as pickle
import numpy as np
import pdb,os
from _collections import defaultdict

from util.make_movies_mito_cbio import ColorMap

couleurs=['green', 'red', 'yellow', 'blue', 'purple', 'orange', 'black', 'cyan']
couleurs.append("#9e0142")
couleurs.append("#5e4fa2")
DRUGS=['Acyclovir',
     'Adenosine3',
     'Aminopurine6benzyl',
     'Anisomycin',
     'Azacytidine5',
     'Camptothecine(S,+)',
     'Daunorubicinhydrochloride',
     'Dexamethasoneacetate',
     'Doxorubicinhydrochloride',
     'Epiandrosterone',
     'Etoposide',
     'Hesperidin',
     'Idoxuridine',
     'JNJ7706621',
     'MLN8054',
     'Methotrexate',
     'Nocodazole',
     'Paclitaxel',
     'R763',
     'Ribavirin',
     'Sulfaguanidine',
     'Sulfathiazole',
     'Thalidomide',
     'VX680',
     'Zidovudine,AZT']

CLASSES=['Interphase',
     'Large',
     'Elongated',
     'Binucleated',
     'Polylobed',
     'Grape',
     'Metaphase',
     'Anaphase',
     'MetaphaseAlignment',
     'Prometaphase',
     'ADCCM',
     'Apoptosis',
     'Hole',
     'Folded',
     'SmallIrregular']

def plotPrep(file_='/media/lalil0u/New/projects/drug_screen/results/MDS_Mitocheck_DS_distances_cost2_10.pkl'):
    '''
    File without plate 4
'''
    f=open(file_, 'r')
    r,who=pickle.load(f); f.close()
    
    f=open('../data/mitocheck_exp_hitlist_perPheno.pkl')
    hitlistperpheno=pickle.load(f)
    f.close()
    
    colors=[]
    type_=defaultdict(list)
    for phenotype in hitlistperpheno:
        for el in hitlistperpheno[phenotype]:
            type_[el].append(phenotype)
            
    phenotypes=hitlistperpheno.keys()
    for el in who[:6214]:
        colors.append(couleurs[phenotypes.index(type_[el][0])])
        
    f=open('/media/lalil0u/New/projects/drug_screen/results/well_drug_dose.pkl')
    genes, drugs, doses, doses_cont, _=pickle.load(f)
    f.close()
    
    exposure=[]
    ord_doses=[]
    for exp in who[:6214]:
        if len(type_[exp])<2:
            exposure.append(type_[exp][0])
        else:
            exposure.append('NO')
    for exp in who[6214:]:
        exposure.append(drugs[exp])
        ord_doses.append(doses_cont[exp])
        
        if drugs[exp]=='empty':
            colors.append('green')
        else:
            colors.append('red')
        
    return colors, drugs, r, who, phenotypes, np.array(exposure), np.array(ord_doses)


def globalPlot(colors, drugs,r, who, phenotypes):
    f,axes=p.subplots(1,2,sharex=True, sharey=True)
    for i in range(6214):
        axes[0].scatter(r[i,0], r[i,1],color=colors[i], s=5)
    for i in range(6214, r.shape[0]):
        if drugs[who[i]]=='empty':
            axes[0].scatter(r[i,0], r[i,1],color='grey', s=5, marker='+')
    axes[0].scatter(0,0,color='grey',label='DS control',marker='+')
    for pheno in phenotypes:
        axes[0].scatter(0,0,color=couleurs[phenotypes.index(pheno)], label=pheno, s=5)
    for i in range(6214, r.shape[0]):
        axes[1].scatter(r[i,0], r[i,1],color=colors[i], s=5)
    axes[1].scatter(0,0,color='red',label='DS')
    axes[1].scatter(0,0,color='green',label='DS ctrl')
    axes[1].legend()
    axes[0].legend()
    p.show()

def distinctDrugPlots(colors, drugs,r, who, phenotypes, exposure):
    f,axes=p.subplots(3,4,sharex=True, sharey=True)
    for i in range(6214):
        if exposure[i] in phenotypes[:4]:
            axes[0,phenotypes.index(exposure[i])].scatter(r[i,0], r[i,1],color=colors[i], s=5)
        elif exposure[i] in phenotypes[4:]:
            axes[1,phenotypes.index(exposure[i])-4].scatter(r[i,0], r[i,1],color=colors[i], s=5)
    
    for i,pheno in enumerate(phenotypes):
        if i<4:
            axes[0,i].scatter(0,0,color=couleurs[phenotypes.index(pheno)], label=pheno, s=5)
        else:
            axes[1,i-4].scatter(0,0,color=couleurs[phenotypes.index(pheno)], label=pheno, s=5)
    
    for i,el in enumerate(who[6214:]):
        if drugs[el] in DRUGS[:9]:
            axes[2,1].scatter(r[6214+i,0], r[6214+i,1],color=couleurs[DRUGS.index(drugs[el])], s=5)
        elif drugs[el] in DRUGS[9:18]:
            axes[2,2].scatter(r[6214+i,0], r[6214+i,1],color=couleurs[DRUGS.index(drugs[el])-9], s=5)
        elif drugs[el] in DRUGS[18:]:
            axes[2,3].scatter(r[6214+i,0], r[6214+i,1],color=couleurs[DRUGS.index(drugs[el])-18], s=5)
        else:
            axes[2,0].scatter(r[i,0], r[i,1],color='grey', s=5, marker='+')
            
    for i,drug in enumerate(DRUGS):
        if i<9:
            axes[2,1].scatter(0,0, color=couleurs[DRUGS.index(drug)], s=5, label=drug)
        elif i>=18:
            axes[2,3].scatter(0,0, color=couleurs[DRUGS.index(drug)-18], s=5, label=drug)
        else:
            axes[2,2].scatter(0,0, color=couleurs[DRUGS.index(drug)-9], s=5, label=drug)
            
    for ax in axes:
        for el in ax:
            el.legend(prop={'size':8})
    p.show()
    
def distinctDrugBoxplots_PERC(who, exposure,doses, perc, phenotypes):
    '''
    We suppose that plate 4 is already removed from r and who
    Here we're directly looking at phenotypes percentages over the whole movie, compared with controls and Mitocheck hits
    Percentages are in the file /media/lalil0u/New/projects/drug_screen/results/all_Mitocheck_DS_phenohit.pkl
'''
    cm = ColorMap()
    cr = cm.makeColorRamp(256, ["#FFFF00", "#FF0000"])
    degrade = [cm.getColorFromMap(x, cr, 0, 10) for x in range(11)]
    
    f,axes=p.subplots(4,9, sharex=True, sharey=True)
    
    for i, pheno in enumerate(phenotypes):
        axes.flatten()[i].boxplot([perc[np.where(exposure==pheno),k] for k in range(perc.shape[1])])
        axes.flatten()[i].set_title(pheno)
        axes.flatten()[i].set_ylim(-0.05,0.9)
        
    dd=np.zeros(shape=6214); dd.fill(-1)
    doses=np.hstack((dd, doses))
    
    for j, drug in enumerate(DRUGS):
        for dose in range(10):
            where_=np.where((exposure==drug)&(doses==dose))[0]
            if where_.shape[0]>0: 
                for k in range(perc.shape[1]):
                    axes.flatten()[8+j].scatter([k+1 for x in range(where_.shape[0])], perc[where_,k], color=degrade[dose], alpha=0.5, s=5)
        axes.flatten()[8+j].set_title(drug)
    print pheno
    where_=np.where(exposure=='empty')[0]
    axes.flatten()[8+j+1].boxplot([perc[where_,k] for k in range(perc.shape[1])])
    axes.flatten()[8+j+1].set_title('Control')
    axes.flatten()[8+j+1].set_xticklabels(CLASSES, rotation='vertical')
    p.show()
    
def distinctDrugBoxplots_PHENOSCORE(who_ps, res, ctrl_points, folder='/media/lalil0u/New/projects/drug_screen/results/'):
    '''
    Here we're looking at phenotypic scores, compared with controls and Mitocheck hits
'''
    cm = ColorMap()
    cr = cm.makeColorRamp(256, ["#FFFF00", "#FF0000"])
    degrade = [cm.getColorFromMap(x, cr, 0, 10) for x in range(11)]

    f=open('/media/lalil0u/New/projects/drug_screen/results/well_drug_dose.pkl')
    _, drugs, _, doses_cont, _=pickle.load(f)
    f.close()
    
    exposure=[]
    doses=[]

    for exp in who_ps:
        exposure.append(drugs["{}--{:>05}".format(exp.split('--')[0], int(exp.split('--')[1]))])
        doses.append(doses_cont["{}--{:>05}".format(exp.split('--')[0], int(exp.split('--')[1]))])
    exposure=np.array(exposure); doses=np.array(doses)
    for drug in DRUGS:
        f,axes=p.subplots(4,4, sharex=True, figsize=(24,12))
        for k,class_ in enumerate(CLASSES):
            axes.flatten()[k].boxplot(ctrl_points[:,k])
            for dose in range(10):
                where_=np.where((exposure==drug)&(doses==dose))[0]
                x=np.random.normal(1, 0.05, size=where_.shape[0])
                if where_.shape[0]>0: 
                    axes.flatten()[k].scatter(x, res[where_,k], color=degrade[dose], alpha=0.8, s=6)
            
            axes.flatten()[k].set_title(class_)
        
        p.title(drug)
        p.savefig(os.path.join(folder, 'phenoscore_{}.png'.format(drug)))
        
    p.close('all')
        
def distinctPhenoPlot(res, ctrl_points):
    '''
   Quick plots to see distributions of phenotypic scores as a function of phenotype, with different colors for control and experiment 
'''
    
    colors=['red' for k in range(res.shape[0])]
    colors.extend(['green' for k in range(ctrl_points.shape[0])]); colors=np.array(colors)
    res=np.vstack((res, ctrl_points))
    
    f,axes=p.subplots(4,4)
    for k in range(res.shape[1]):
        ord=np.argsort(res[:,k])
        loc_col=colors[ord]
        axes.flatten()[k].boxplot(res[ord, k][np.where(loc_col=='green')])
        x=np.random.normal(1, 0.08, size=np.where(loc_col=='red')[0].shape[0])
        axes.flatten()[k].plot(x, res[ord, k][np.where(loc_col=='red')],'r.', color='red', alpha=0.2)
        
        axes.flatten()[k].set_title(CLASSES[k])
        
    p.show()
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        