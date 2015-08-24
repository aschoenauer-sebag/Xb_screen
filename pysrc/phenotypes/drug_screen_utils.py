import matplotlib.pyplot as p
import matplotlib as mpl
import cPickle as pickle
import numpy as np
import pdb,os
from _collections import defaultdict

from util.make_movies_mito_cbio import ColorMap
from scipy.stats.stats import scoreatpercentile

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

lim_Mito=6452


'''
Just some things not to forget when doing the distances on phenotypic scores:
- we wanted to take out Anapahses and Interphases
- we need to normalize the different phenotypic scores so that they're comparable for different phenotypes
'''

def plotExternalConsistency(corr_dict, labels, cmap=mpl.cm.bwr):
    norm = mpl.colors.Normalize(0.5,1)
    corr=np.vstack((corr_dict[el] for el in corr_dict))
    
    f=p.figure()
    ax=f.add_subplot(111)
    ax.matshow(corr, cmap=cmap, norm=norm)
    p.xticks(range(len(labels)), labels, rotation='vertical')
    p.yticks(range(len(corr_dict)), [el for el in corr_dict])
    
    p.show()


def plotInternalConsistency(M, tick_labels, cmap=mpl.cm.bwr):
    
    norm = mpl.colors.Normalize(0,scoreatpercentile(M.flatten(), 95))
    f=p.figure()
    ax=f.add_subplot(111)
    ax.matshow(M, cmap=cmap, norm=norm)
    p.yticks(range(0,M.shape[0],2),tick_labels[::2])
    ax.tick_params(labelsize=6)
    p.show()
    
    return


def computing_Mito_pheno_scores(pheno_scores_dict, exp_list):
    r=np.zeros(shape=(len(exp_list), len(CLASSES)))
    
    for i,exp in enumerate(exp_list):
        for j, class_ in enumerate(CLASSES):
            if class_=='Binucleated':
                r[i,j]=pheno_scores_dict[exp]['max_{}'.format('Shape1')]
            elif class_=='Polylobed':
                r[i,j]=pheno_scores_dict[exp]['max_{}'.format('Shape3')]
            else:
                r[i,j]=pheno_scores_dict[exp]['max_{}'.format(class_)]
                
    return r


def plotPrep(file_='/media/lalil0u/New/projects/drug_screen/results/MDS_Mitocheck_DS_distances_0.1.pkl'):
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
    for el in who[:lim_Mito]:
        colors.append(couleurs[phenotypes.index(type_[el][0])])
        
    f=open('/media/lalil0u/New/projects/drug_screen/results/well_drug_dose.pkl')
    genes, drugs, doses, doses_cont, _=pickle.load(f)
    f.close()
    
    exposure=[]
    ord_doses=[]
    for exp in who[:lim_Mito]:
        if len(type_[exp])<2:
            exposure.append(type_[exp][0])
        else:
            exposure.append('NO')
    for exp in who[lim_Mito:]:
        e='{}--{:>05}'.format(exp.split('--')[0], int(exp.split('--')[1]))
        exposure.append(drugs[e])
        ord_doses.append(doses_cont[e])
        
        if drugs[e]=='empty':
            colors.append('green')
        else:
            colors.append('red')
        
    return colors, drugs, r, who, phenotypes, np.array(exposure), np.array(ord_doses)


def globalPlot(colors, exposure,r, who, phenotypes):
    '''
   MDS visu 
'''
    f,axes=p.subplots(1,2,sharex=True, sharey=True)
    for i in range(lim_Mito):
        axes[0].scatter(r[i,0], r[i,1],color=colors[i], s=5)
    for i in range(lim_Mito, r.shape[0]):
        if exposure[i]=='empty':
            axes[0].scatter(r[i,0], r[i,1],color='grey', s=5, marker='+')
    axes[0].scatter(0,0,color='grey',label='DS control',marker='+')
    for pheno in phenotypes:
        axes[0].scatter(0,0,color=couleurs[phenotypes.index(pheno)], label=pheno, s=5)
    for i in range(lim_Mito, r.shape[0]):
        axes[1].scatter(r[i,0], r[i,1],color=colors[i], s=5)
    axes[1].scatter(0,0,color='red',label='DS')
    axes[1].scatter(0,0,color='green',label='DS ctrl')
    axes[1].legend()
    axes[0].legend()
    p.show()

def distinctDrugPlots(colors, r, who, phenotypes, exposure):
    '''
       MDS visu 
'''
    f,axes=p.subplots(3,5)
    for i in range(lim_Mito):
        if exposure[i] in phenotypes[:5]:
            axes[0,phenotypes.index(exposure[i])].scatter(r[i,0], r[i,1],color=colors[i], s=5)
        elif exposure[i] in phenotypes[5:]:
            axes[1,phenotypes.index(exposure[i])-5].scatter(r[i,0], r[i,1],color=colors[i], s=5)
    
    for i,pheno in enumerate(phenotypes):
        if i<5:
            axes[0,i].scatter(0,0,color=couleurs[phenotypes.index(pheno)], label=pheno, s=5)
        else:
            axes[1,i-5].scatter(0,0,color=couleurs[phenotypes.index(pheno)], label=pheno, s=5)
    
    for i,el in enumerate(who[lim_Mito:]):
        j=i+lim_Mito
        if exposure[j] in DRUGS[:7]:
            axes[2,1].scatter(r[j,0], r[j,1],color=couleurs[DRUGS.index(exposure[j])], s=5)
        elif exposure[j] in DRUGS[7:14]:
            axes[2,2].scatter(r[j,0], r[j,1],color=couleurs[DRUGS.index(exposure[j])-7], s=5)
        elif exposure[j] in DRUGS[14:21]:
            axes[2,3].scatter(r[j,0], r[j,1],color=couleurs[DRUGS.index(exposure[j])-14], s=5)
        elif exposure[j] in DRUGS[21:]:
            axes[2,4].scatter(r[j,0], r[j,1],color=couleurs[DRUGS.index(exposure[j])-21], s=5)
        else:
            axes[2,0].scatter(r[i,0], r[i,1],color='grey', s=5, marker='+')
            
    for i,drug in enumerate(DRUGS):
        if i<7:
            axes[2,1].scatter(0,0, color=couleurs[i], s=5, label=drug)
        elif i>=21:
            axes[2,4].scatter(0,0, color=couleurs[i-21], s=5, label=drug)
        elif i>=7 and i<14:
            axes[2,2].scatter(0,0, color=couleurs[i-7], s=5, label=drug)
        else:
            axes[2,3].scatter(0,0, color=couleurs[i-14], s=5, label=drug)
            
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
        
    dd=np.zeros(shape=lim_Mito); dd.fill(-1)
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
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        