import string, getpass, datetime, os, pdb
from warnings import warn
from collections import defaultdict

from scipy.stats import nanmean, nanstd
import numpy as np
import cPickle as pickle
import vigra.impex as vi

if getpass.getuser()=='aschoenauer':
    import matplotlib as mpl
    mpl.use('Agg')
    show = False
elif getpass.getuser()=='lalil0u':
    import matplotlib as mpl
    show=True

from util import settings
from tracking.importPack.imp import importTargetedFromHDF5

from analyzer import COLORD_DS, COLORD_XB
from plateSetting import readPlateSetting

from util.plots import makeColorRamp, couleurs, basic_colors
from util.movies import makeMovieMultiChannels, makeMovie
import brewer2mpl
import matplotlib.pyplot as p

from plates.models import Well


class HTMLGenerator():
    def __init__(self, settings_file):
        self.settings = settings.Settings(settings_file, globals())
        
        self.ap = ArrayPlotter(plotDir=self.settings.plot_dir, 
                               legendDir=self.settings.plot_dir) 
        self.wp = WellPlotter(plotDir=self.settings.plot_dir)
        
    def intensity_qc(self, plateL):
        result=defaultdict(dict)
        for plate in plateL:
            well_setup=self.well_lines_dict[plate].flatten()
            result[plate]=defaultdict(dict)
            for well in well_setup[np.where(well_setup>0)]:
                #even if there are missing columns don't forget that the image directory uses Zeiss counting
                imgDir= os.path.join(self.settings.raw_data_dir, plate, 'W{:>05}'.format(well))
                try:
                    l=os.listdir(imgDir)
                except OSError:
                    print "No image for well ", np.where(well_setup==well)[0][0]+1, 'plate ', plate
                    continue
                else:
                    if l==[]:
                        print "No image for well ", np.where(well_setup==well)[0][0]+1, 'plate ', plate
                        continue
                    else:
                        if not os.path.isdir(os.path.join(self.settings.plot_dir, plate)):
                            os.mkdir(os.path.join(self.settings.plot_dir, plate))
                        if 'mean_intensity_{}--W{:>05}.png'.format(plate, np.where(well_setup==well)[0][0]+1) in os.listdir(os.path.join(self.settings.plot_dir,plate)):
                            continue
                        
                        mean_intensity=np.array([np.mean(vi.readImage(os.path.join(imgDir, im))) for im in sorted(filter(lambda x: 'c00002' in x, l))])
                        dev_intensity=np.std(mean_intensity)
                        local_mean=[np.mean(mean_intensity[k:k+10]) for k in range(mean_intensity.shape[0]-10)]
                        local_mean.extend(local_mean[-10:]); local_mean=np.array(local_mean)
         #so here the images to delete are the ones where the intensity is high above the intensity before, or high below 
         #the frame numbers go from 0 to the end of the movie
                        toDel=np.where(np.abs(mean_intensity-local_mean)>3*dev_intensity)[0]
                        toDelFin=toDel
                        
                        f=p.figure(); ax=f.add_subplot(111)
                        ax.scatter(range(len(mean_intensity)), mean_intensity, color='green', label='Mean frame intensity')
                        ax.plot(range(len(mean_intensity)), local_mean+3*dev_intensity, color='red', label='Sliding average+3sigmas')
                        ax.plot(range(len(mean_intensity)), local_mean-3*dev_intensity, color='red', label='Sliding average-3sigmas')
                        ax.set_ylim((0,1000)); ax.legend()
                        p.savefig(os.path.join(self.settings.plot_dir,plate, 'mean_intensity_{}--W{:>05}.png'.format(plate, np.where(well_setup==well)[0][0]+1)))
                        p.close(f)
                        result[plate][well]=toDelFin
        if result[plateL[0]]!={}:
            try:
                f=open(self.settings.intensity_qc_filename, 'r')
                d=pickle.load(f);f.close()
            except IOError:
                d=result
            else:
                for plate in result:
        #doing this on a per plate basis because otherwise if we have already done some of the wells for a given plate it is going to be erased in a global update command
                    try:
                        d[plate].update(result[plate])
                    except KeyError:
                        d.update(result)

            f=open(self.settings.intensity_qc_filename, 'w')
            pickle.dump(d,f);f.close()
        
        return
        
    def targetedDataExtraction(self, plateL, featureL):
        '''
       If there are more than one position per well they are all taken into account here
'''
        newFrameLot = None
        dataFolder = self.settings.raw_result_dir  
        print "Looking for features ", featureL
        for plate in plateL:
            listW = sorted(filter(lambda x: '.hdf5' in x or '.ch5' in x, os.listdir(os.path.join(dataFolder, plate, 'hdf5'))))
            for filename in listW[:5]:
                well=filename.split('.')[0]
                    
                filename = os.path.join(dataFolder, plate,"hdf5", filename)
                try:
                    featureL,classes, frameLotC= importTargetedFromHDF5(filename, plate, well,featureL,primary_channel_name=self.settings.primary_channel_name,
                                                                        secondary=self.settings.secondaryChannel)
                except ValueError:
                    print "Error at loading from hdf5 ", plate, well
                    continue
                if newFrameLot == None:
                    newFrameLot = frameLotC 
                else: newFrameLot.addFrameLot(frameLotC)
        return featureL,classes, newFrameLot 
    
    def featureComparison(self, table, index):
        primary = table[:,index]; secondary =table[:,self.FEATURE_NUMBER+index]+0.00001
        return np.divide(primary,secondary)
    
    def featureHist(self, array, bins, binOfInterest):
        h1,_ = np.histogram(array, bins = bins)
        return h1[binOfInterest]/float(np.sum(h1))
        
    def formatData(self, frameLot, resD, featureL,featureChannels):
        '''
       Need to modify it to take into account all positions for the same wells
'''
        featureL = list(featureL)
        if self.classes is not None:
            featureChannels.extend([0 for k in range(len(filter(lambda x: x not in featureL, self.classes)))])
            featureL.extend(filter(lambda x: x not in featureL, self.classes))
            self.settings.well_features.append('Nuclear_morphologies')
            
        print'--------------', featureL
        
        for plate in frameLot.lstFrames:
            if plate not in resD:
                print plate, resD.keys()#just checkin
                warn("Plate {} not in result dictionary from CSV file extraction".format(plate))
                continue
            for well in frameLot.lstFrames[plate]:
                secondary=frameLot.lstFrames[plate][well][frameLot.lstFrames[plate][well].keys()[0]].secondary
                
                if int(well[:-3]) not in resD[plate]:
                    warn("Well {} not in result dictionary from CSV file extraction".format(well))
                    continue
                
                result={'object_count':[], 'cell_count':[]}
                if secondary:
                    result.update({'red_only':[], 'circularity':[]})
                
                argL = zip(featureL, featureChannels)
                result.update({'{}_ch{}'.format(arg[0], arg[1]+1):[] for arg in argL})
                
                for frame_nb in frameLot.lstFrames[plate][well]:
                    frame = frameLot.lstFrames[plate][well][frame_nb]
                    
                #for each frame, the cell count = object count - [artefacts + out of focus objects]. Only possible if the classification was computed
                    if self.classes is not None:
                        bincount = np.bincount(frame.classification, minlength=len(self.classes))
                        if np.sum(bincount)==0:
                            print "No info frame {} well {} plate {}".format(frame_nb, well, plate)
                            continue
                        out_of_focus = bincount[np.where(self.classes==self.settings.focus_classname)] if self.settings.focus_classname in self.classes else 0
#NB: no more smallUnidentified after classifier update on 2/20/15
                        artefacts = bincount[np.where(self.classes==self.settings.artefact_classname)] if self.settings.artefact_classname in self.classes else 0
                        result["cell_count"].append( frame.centers.shape[0] - (out_of_focus + artefacts))
                        
                #for each frame, the object count = all objects = cells+out of focus objects + artefacts
                    result["object_count"].append(frame.centers.shape[0])
                    
                #if second channel was processed, computing percentage of red only cells
                    if secondary:
                        try:
                            red_only_dist =self.featureComparison(frame.features, featureL.index("roisize"))
                            result["red_only"].append(self.featureHist(red_only_dist, bins = [0, 0.95, 1], binOfInterest = 1))
                        except:
                            raise ValueError
                #other parameters are specified in the setting file    
                    for arg in argL:
                        if arg[0]=='irregularity':
                            result['{}_ch{}'.format(arg[0], arg[1]+1)].append(frame.features[:,arg[1]*self.FEATURE_NUMBER+featureL.index(arg[0])]+1)
                            continue
                        
                        if self.classes is not None and arg[0] in self.classes:
                            result['{}_ch1'.format(arg[0])].append(bincount[np.where(self.classes==arg[0])]/float(np.sum(bincount)))
                        else:
                            try:
                                result['{}_ch{}'.format(arg[0], arg[1]+1)].append(frame.features[:,arg[1]*self.FEATURE_NUMBER+featureL.index(arg[0])])
                            except:
                                if not secondary and arg[1]+1==2:
                                    continue
                                else:
                                    raise
                
                #if the second channel was processed, looking at round cells from a cytoplasmic point of view because it's the cells that are dying (or dividing but here we neglect that)
                    if "circularity_ch2" in result:
                    #circularity at well level
                        result['circularity'].append(self.featureHist(result["circularity_ch2"][-1], bins = [1, 1.5, 5], binOfInterest = 0))
                        
                result["initObjectCount"]=result["object_count"][0]
                result["initCellCount"]=result["cell_count"][0]
                result["proliferation"]=result["cell_count"][-1]/float(result["cell_count"][0])

                if secondary:
                    result["initNucleusOnly"]=result["red_only"][0]
                    result["endNucleusOnly"]=result["red_only"][-1]
                    result['initCircularity']= self.featureHist(result["circularity_ch2"][0], bins = [1, 1.4, 5], binOfInterest = 0)
                    result['endCircularity'] = self.featureHist(result["circularity_ch2"][-1], bins = [1, 1.4, 5], binOfInterest = 0)
             
                    result['initCircMNucleus']= result['initCircularity']-result["initNucleusOnly"]
                    result['endCircMNucleus']= result['endCircularity']-result["endNucleusOnly"]            
                    result["death"]=result['endCircMNucleus']/float(result['initCircMNucleus']+0.0001)
                    
                if self.classes is not None:
            #computing the percentage of out of focus nuclei on the last image
                    result['endFlou']=result[self.settings.focusFeature][-1]
                
                for el in self.classes:
                    result['{}_ch1'.format(el)]=np.array(result['{}_ch1'.format(el)])[:,0]
                    
                if int(well[-2:])==1:
                    resD[plate][int(well[:-3])].update(result)
                else:
                    resD[plate][int(well[:-3])]['pos_02']=result
        return 1
    
    def count_qc(self, failed_qc, plate, resD):
        qc_init_cell = self.settings.qc_init_cell
        qc_end_OOF = self.settings.qc_end_OOF
        
        try:
            resCour = resD[plate]
        except KeyError:
            print plate, ' not in result dictionary.'
            return
        
        for well in resCour:
            if 'object_count' not in resCour[well] and 'cell_count' not in resCour[well]:
                #No info to do the QC
                continue
            
            if 'cell_count' not in resCour[well]:
                #do count QC on object_count because we don't have the classification
                count=np.mean(resCour[well]['object_count'][:10])
                if 'pos_02' in resCour[well]:
                    #we take the two positions into account
                    count+=np.mean(resCour[well]['pos_02']['object_count'][:10])
                    
                if count<qc_init_cell:
                    print "Failing initial object count QC ", plate, well
                    failed_qc[plate].append(well)
                    continue
                
            else:
                #do count QC on cell_count
                count=np.mean(resCour[well]['cell_count'][:10])
                focus =np.mean(resCour[well][self.settings.focusFeature][-10:])
                
                if 'pos_02' in resCour[well]:
                    count+=np.mean(resCour[well]['pos_02']['cell_count'][:10])
                    focus+=np.mean(resCour[well]['pos_02'][self.settings.focusFeature][-10:])
                    
                if count<qc_init_cell:
                    print "Failing initial cell count QC ", plate, well
                    failed_qc[plate].append(well)
                    continue
                
                if focus>qc_end_OOF:
                    print "Failing out of focus end count QC ", plate, well
                    failed_qc[plate].append(well)
                    continue
                
        return 1
        
        
    def saveResults(self,failed_qc, plate, resD):
        '''
        Saving results for data extraction as well as quality control results
        '''
        try:
            resCour = resD[plate]
        except KeyError:
            print plate, ' not in result dictionary.'
            return
        resCour['FAILED QC'] = failed_qc[plate]
        
        if not os.path.isdir(self.settings.result_dir):
            os.mkdir(self.settings.result_dir)
        
        f=open(os.path.join(self.settings.result_dir, 'processedDictResult_P{}.pkl'.format(plate)), 'w')
        pickle.dump(resCour, f)
        f.close()
        return
    
    def plotALaDemande(self, outputInterest, parametersInterest, plateL=None, divideByInitVal=False, plotMoy=True):
        '''
        Here the idea is to plot cellCount or any well feature called outputInterest, as a function of some parameters of interest,
        parametersInterest.
        
        Input:
        - outputInterest : string
        - parametersInterest: list
        - divideByInitVal: bool, say if the values should be divided by the initial value (makes sense if it's cell count for example)
        - plotMoy: bool, say if you want one curve per condition, or as many curves as there are wells
        
        '''
        if plateL is None:
            results = filter(lambda x: 'processedDictResult' in x, os.listdir(self.settings.result_dir))
        else:
            results = filter(lambda x: x in ['processedDictResult_P{}.pkl'.format(plate) for plate in plateL],os.listdir(self.settings.result_dir)) 
        if results==[]:
            self.__call__(plateL)
            return self.plotALaDemande(outputInterest, parametersInterest, plateL)
        else:
            f=p.figure(figsize=(24,13))
            ax=f.add_subplot(111)
            paramIntVal=[]
            paramIntDict={}
            for result in results:
                f=open(os.path.join(self.settings.result_dir, result), 'r')
                resCour = pickle.load(f); f.close()
                plate = result.split("_")[1][1:6]
                for well in resCour:
                    try:
                        legendName = ''
                        if 'CloneA' in resCour[well]['Name']:
                            print 'Not taking Clone A into account'
                            continue
                        if 'Complet' in resCour[well]['Medium']:
                            print 'Not taking whole medium into account'
                            continue
                        for parameterInterest in parametersInterest:
                            legendName +=' '+resCour[well][parameterInterest].strip(' ').split('_')[0]
                        legendName+=' plate {} wells'.format(max(resCour.keys()))
                        if legendName not in paramIntVal:
                            paramIntVal.append(legendName)
                            legend = True
                            paramIntDict[legendName]=None
                        else:
                            legend=False
                    except KeyError:
                        print 'Parameter {} not in result for well {}'.format(parameterInterest, well)
                        #pdb.set_trace()
                    else:
                        try:
                            if divideByInitVal:
                                val = np.array(resCour[well][outputInterest])/float(resCour[well][outputInterest][0])
                            else:
                                val = np.array(resCour[well][outputInterest])
                            if val.shape[0]<180:
                                print well, ' taille inf a 150'
                                continue

                        except KeyError:
                            print 'Output {} not in results for well {}, plate {}'.format(outputInterest, well, plate)
                            if legend:
                                paramIntVal.pop()
                        else:
#                            if np.any(np.isnan(val)):
#                                pdb.set_trace()
                            
                            if plotMoy:
                                if paramIntDict[legendName] is None:
                                    paramIntDict[legendName]=val
                                else:
                                    try:
                                        paramIntDict[legendName]=np.vstack((paramIntDict[legendName], val))
                                    except:
                                        if val.shape[0]>paramIntDict[legendName].shape[1]:
                                            paramIntDict[legendName]=np.vstack((paramIntDict[legendName], val[:paramIntDict[legendName].shape[1]]))
                                        elif val.shape[0]<paramIntDict[legendName].shape[1]:
                                            paramIntDict[legendName]=np.vstack((paramIntDict[legendName][:,:val.shape[0]], val))
                            else:                 
                                if legend:
                                    ax.scatter(range(len(resCour[well]['object_count'])), val, 
                                        color=couleurs[paramIntVal.index(legendName)],
                                        label = legendName)
                                else:
                                    ax.scatter(range(len(resCour[well]['object_count'])), val,
                                        color=couleurs[paramIntVal.index(legendName)])
                                ax.text(len(resCour[well]["object_count"]), val[-1], 
                                        'p{} w{}'.format(plate, int(well)), fontsize=10)

            if plotMoy:
                length = np.min(np.array([paramIntDict[legendName].shape[1] for legendName in paramIntDict]))
                for paramVal in paramIntVal:
                    if paramIntVal.index(paramVal)>=len(couleurs):
                        couleurs.extend(basic_colors)
                    ax.errorbar(range(length), nanmean(paramIntDict[paramVal], 0)[:length], nanstd(paramIntDict[paramVal], 0)[:length],
                                color=couleurs[paramIntVal.index(paramVal)],
                                label = paramVal)
            
                ax.set_xlim(0,200)
                if divideByInitVal:
                    ax.set_ylim(0.5,3)
                else:
                    ax.set_ylim(0,1)
                ax.grid(True)
                ax.legend()
                ax.set_title("Plot param {} output {}/value at t=0".format(parametersInterest, outputInterest))
                if show:
                    p.show()
                else:
                    legendName = ''
                    for parameterInterest in parametersInterest:
                            legendName +=' '+parameterInterest
                    p.savefig(os.path.join(self.settings.result_dir, 'im_{}_{}.png'.format(outputInterest, legendName)))
            

    
    def generatePlatePlots(self, plate, resD, failed_qc):
        
        try:
            resCour = resD[plate]
        except KeyError:
            print plate, ' not in result dictionary.'
            return
        

        settingL = [
            {'min': self.settings.density_plot_settings['min_count'], 'max': self.settings.density_plot_settings['max_count'], 'int_labels': True},
            {'min': self.settings.density_plot_settings['min_count'], 'max': self.settings.density_plot_settings['max_count'], 'int_labels': True},
            {'min': self.settings.density_plot_settings['min_proliferation'], 'max': self.settings.density_plot_settings['max_proliferation'], 'int_labels': False},
            
            {'min': self.settings.density_plot_settings['min_circularity'], 'max': self.settings.density_plot_settings['max_circularity'], 'int_labels': False},
            {'min': self.settings.density_plot_settings['min_circularity'], 'max': self.settings.density_plot_settings['max_circularity'], 'int_labels': False},
            
            {'min': self.settings.density_plot_settings['min_circularity'], 'max': self.settings.density_plot_settings['max_circularity'], 'int_labels': False},
            {'min': self.settings.density_plot_settings['min_circularity'], 'max': self.settings.density_plot_settings['max_circularity'], 'int_labels': False},
            
            {'min': self.settings.density_plot_settings['min_circularity'], 'max': self.settings.density_plot_settings['max_circularity'], 'int_labels': False},
            {'min': self.settings.density_plot_settings['min_circularity'], 'max': self.settings.density_plot_settings['max_circularity'], 'int_labels': False},
            {'min': self.settings.density_plot_settings['min_class'], 'max': self.settings.density_plot_settings['max_class'], 'int_labels': False}
            ]
        
        plotDir = os.path.join(self.settings.plot_dir, self.settings.newPlateName(plate))
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)
#        try:
#            missing_col = self.settings.missing_cols[plate]
#        except KeyError:
#            missing_col = []
        
        ###printing plate setup
        print "working on plate setup"
        filename = '%s--%s.png' % ('plate', self.settings.newPlateName(plate))
        texts = self.ap.plotPlate(resCour,self.params, filename, where_wells = self.well_lines_dict[plate],plotDir=plotDir, title='{} {}'.format(plate, 'plate'),
                                  qc_info=failed_qc[plate])
        
        for pheno, setting in zip(['initObjectCount', 'initCellCount', 'proliferation', 
                                   'initCircularity', 'endCircularity',  
                                   "initNucleusOnly", "endNucleusOnly",
                                   'initCircMNucleus', 'endCircMNucleus', 'endFlou'], settingL):
            if pheno not in set([el for k in resCour for el in resCour[k]]):
                print "No {} for plate {}".format(pheno, plate)
                continue
            
            print "working on ", pheno
            filename = '%s--%s.png' % (pheno, self.settings.newPlateName(plate))
            
            data = self.ap.prepareData(resCour, self.ap.getPhenoVal(pheno), self.well_lines_dict[plate])
            
            try:
                if np.isnan(data):
                    continue
            except:
                pass
            #absent ce sont les puits ou on n'a pas de data mais on a des images, ie on devrait avoir des infos
            absent = np.where(data==-1)
            data[absent]=0
            self.ap.plotArray(data, filename, plotDir=plotDir,
                title='{} {}'.format(plate, pheno),
                absent= absent,
                int_labels = setting['int_labels'],
                norm = True,
                min_val=setting['min'], max_val=setting['max'], 
                label_colors = False, #labeledPosD=labeledPosD,
                legend_plot = True, legend_filename='legend_%s' % filename,
                texts=texts)
            p.close('all')
#TODO investigate this normalization story
        return 
    
    def generateWellPlots(self, plate, resD, colorDict):
        try:
            resCour = resD[plate]
        except KeyError:
            print plate, ' not in result dictionary.'
            return
        well_setup=self.well_lines_dict[plate].ravel()
        plotDir = os.path.join(self.settings.plot_dir, self.settings.newPlateName(plate))
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        for well in resCour:
            #we want to plot cell count and circularity evolutions
            print "working on well ", well
            final_well_num = np.where(well_setup==well)[0][0]+1
            for pheno in self.settings.well_features:
        #So here I enable to take a list of data tables and labels into account for them to be plotted on the same graph
                data=[];labels=[pheno];colors=['blue']
                print "working on ", pheno
                filename = '{}_{}--W{:>05}.png'.format(pheno,self.settings.newPlateName(plate), final_well_num)
                if pheno in list(resCour[well].keys()):
                    data.append(self.wp.prepareData(resCour[well][pheno]))
                    figsize=(8,6)
                elif pheno=='Nuclear_morphologies':
                    labels=self.classes
                    colors=[]
                    figsize=(20,10)
                    for class_ in self.classes:
                        data.append(self.wp.prepareData(resCour[well]['{}_ch1'.format(class_)]))
                        try:
                            colors.append(colorDict['{}'.format(class_)])
                        except KeyError:
                            colors.append('grey')
                else:
                    continue
            
                if pheno =='circularity':
                    data.append(self.wp.prepareData(resCour[well]['red_only']))
                    labels.append("Fraction of cells simply marked")
                    colors.append('red')
                elif pheno=='cell_count':
                    data.append(self.wp.prepareData(resCour[well]['object_count']))
                    labels.append('Object count')
                    colors.append('red')

                try:
                    min_=self.settings.well_plot_settings[pheno][0];max_=self.settings.well_plot_settings[pheno][1]
                except KeyError:
                    max_=None; min_=None
                    
                if pheno=='cell_count':
                    add_line=30
                elif pheno=='Flou_ch1' or pheno=='Focus_ch1':
                    add_line=0.5
                else:
                    add_line=None
                try:
                    self.wp.plotEvolution(final_well_num, filename, plotDir=plotDir,
                        max_=max_, min_=min_, add_line=add_line,
                        title='Plate {}, well {}, evolution of {}'.format(plate, final_well_num, pheno),
                        labels=labels, data=data, colors=colors,figsize=figsize)
                except IndexError:
                    pass
                p.close('all')
        return
    
    def generateMovies(self, plate, well_setup):
        well_setup=well_setup.flatten()
        
        corresp_existing_names = np.array([(x, int(self.settings.whereWell(x))) for x in sorted(os.listdir(os.path.join(self.settings.raw_data_dir, plate)))
                                     if os.path.isdir(os.path.join(self.settings.raw_data_dir, plate,x))])
        for well in well_setup[np.where(well_setup>0)]:
            #even if there are missing columns don't forget that the image directory uses Zeiss counting
            try:
                ind = np.where(corresp_existing_names[:,1]==str(well))[0][0]
                #wellFolder = filter(lambda x: os.path.isdir(os.path.join(self.settings.raw_data_dir, plate,x)) and well==int(self.settings.whereWell(x)), 
                #                   os.listdir(os.path.join(self.settings.raw_data_dir, plate)))[0]
            except IndexError:
                print "No image for well ", np.where(well_setup==well)[0][0]+1, 'plate ', plate
                continue
            else:
                wellFolder = corresp_existing_names[ind][0]
            
            imgDir= os.path.join(self.settings.raw_data_dir, plate, wellFolder)
            try:
                l=os.listdir(imgDir)
            except OSError:
                print "No image for well ", np.where(well_setup==well)[0][0]+1, 'plate ', plate
                continue
            else:
                if l==[]:
                    print "No image for well ", np.where(well_setup==well)[0][0]+1, 'plate ', plate
                    continue
                else:
                    secondary=np.any(np.array(['c00001' in el for el in l]))
                    if secondary:
                        makeMovieMultiChannels(imgDir=imgDir, outDir=self.settings.movie_dir, plate=self.settings.newPlateName(plate), well=np.where(well_setup==well)[0][0]+1, 
                                           ranges=[(50,1000),(200,3000)], secondary=secondary)
                    else:
                        makeMovie(imgDir, outDir=self.settings.movie_dir, plate=self.settings.newPlateName(plate),
                                    well=np.where(well_setup==well)[0][0]+1,
                                    filtre = (lambda x: 'P00001' in x),
                                    clef=(lambda x:int(x.split('--')[3][1:])),
                                #because the images for the drug screen start at 2**15
                                    offset=2**15)
        return    
    
    def changeDBWellNumbers(self,plate, well_setup, idL):
        well_setup =well_setup.flatten()
        try:
            idCour = idL[plate]
        except KeyError:
            print plate, ' not in result dictionary.'
            return
        
        for Zeiss_well_num in idCour:
            db_id = idCour[Zeiss_well_num]
            final_well_num = np.where(well_setup == Zeiss_well_num)[0][0]+1
            w=Well.objects.get(pk=db_id)
            w.num=final_well_num; w.save()
            
        return
        
        
    def __call__(self, plateL, featureL = None, featureChannels =None, saveDB=True, fillDBOnly=False, doQC=False, moviesOnly=False, doMovies=True):
        '''
        Default behaviour regarding movie generation: if already done they are not recomputed
        '''
    #Getting features of interest
        if featureL is None:
            featureL = self.settings.featuresOfInterest
            featureChannels = self.settings.featureChannels
        assert(len(featureL)==len(featureChannels)), "Not same number of features of interest and channels to consider them"
        self.FEATURE_NUMBER = len(featureL)
        
        if self.settings.xbscreen:
            colorDict=COLORD_XB 
            #Adding roisize to be able to assess number of cells that only have fluorescent nucleus
            if 'roisize' not in featureL:
                featureL.append('roisize')
                self.FEATURE_NUMBER +=1
        else:
            colorDict= COLORD_DS

        failed_qc=defaultdict(list)
            
        #first we get the plate setting, not counting empty wells. In particular it means that 
        #wells are counted according to their Zeiss number and not number on the plate.
        #After that step, resD contains only information regarding well treatments
        print ' *** reading plate settings and saving info for plate, well, condition and treatment in db for %s ***' % plateL
        
        #idL is a dictionary {plate:{well_Zeiss_num:well_id in db}} to renumber the wells in case
        #some wells were on the plate were not taken into account while doing plate setup on Zeiss software.
        #in resD, wells are numbered according to Zeiss numbering
        resD, self.params, self.well_lines_dict, idL= readPlateSetting(plateL, self.settings.confDir, 
                                             addPlateWellsToDB=saveDB,
                                             startAtZero = self.settings.startAtZero,
                                             plateName = self.settings.name,
                                             dateFormat= self.settings.date_format,
                                             xbscreen=self.settings.xbscreen
                                             )

        if doQC:
            print ' *** performing well intensity quality control ***'
            self.intensity_qc(plateL)

        if not fillDBOnly:
            if not moviesOnly:
                print ' *** get result dictionary ***'
                featureL,self.classes, frameLot = self.targetedDataExtraction(plateL, featureL)
                
                self.formatData(frameLot, resD, featureL, featureChannels)
                
                #for xb screen only self.classes=np.delete(self.classes, np.where(self.classes=='Polylobbed')[0])
            for plate in plateL:
                if not moviesOnly:
                    print ' *** Quality control for %s: ***' % plate
                    self.count_qc(failed_qc, plate, resD)

                    print ' *** generate density plots for %s: ***' % plate
                    #well_setup is an array representing the plate with the Zeiss well numbers
                    #on the plate. So well_setup.flatten() gives the following: on absolute position
                    #final well_number there is either a zero either the Zeiss well number.
                    
                    #This function is used to generate all plate plots
                    self.generatePlatePlots(plate, resD, failed_qc)
                    
                    print ' *** generate well plots ***'
                    self.generateWellPlots(plate, resD, colorDict)
    
                    print ' *** saving result dictionary ***'
                    self.saveResults(failed_qc, plate, resD)
    
                    if len(np.where(self.well_lines_dict[plate]==-1)[0])>1 and saveDB:
                        print ' *** changing well numbers in db, plate ', plate
                        self.changeDBWellNumbers(plate, self.well_lines_dict[plate], idL)
                if doMovies:
                    #try:
                    print ' *** generate movies ***'
                    self.generateMovies(plate, self.well_lines_dict[plate])
#              
#                     except: 
#                         print ' ERROR while working for %s' % plate
#                         continue
        return
    
    
    
class ArrayPlotter():
    def __init__(self, plotDir, legendDir):
        self.plotDir = plotDir
        if not os.path.isdir(self.plotDir):
            os.makedirs(self.plotDir)
            
        self.legendDir = legendDir
        return
    
    class getPhenoVal(object):

        def __init__(self, phenoClass):
            self.phenoClass = phenoClass

        def __call__(self, valuesD):
            try:
                return valuesD[self.phenoClass]
            except KeyError:
                return -1


    def prepareData(self, resD, getLabel,where_wells):
        '''
        Here we are preparing the data for plotting. 
        
        In the array there's the data from the hdf5 file. Where_wells says where physically on the plate we recorded images.
        I put -2 for those wells, -1 for where I miss data.
        '''
        
        expL = where_wells[np.where(where_wells>-1)]
        hitData=np.empty(shape = where_wells.shape); hitData.fill(-2)
    
        for exp in expL:
            id = exp
            if resD.has_key(id):
                label = getLabel(resD[id])
            else:
                label = -1
            try:
                hitData[np.where(where_wells==exp)]=label
            except TypeError:
                hitData[np.where(where_wells==exp)]=-1
        try:
            print np.min(hitData[np.where(hitData>-1)]), np.max(hitData)
        except:
            return np.nan
        return hitData    
  
    def plotPlate(self, res,params, filename, where_wells, plotDir=None, title='', show=False, qc_info=[]):
        
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
        texts = []
        fig = p.figure(figsize=(12,7))
        ax=fig.add_subplot(111)
        cmap = brewer2mpl.get_map('RdPu', 'Sequential', 9).mpl_colormap
        
        nrow = where_wells.shape[0]
        ncol = where_wells.shape[1]
        
        ax.pcolormesh(np.ones(shape = (nrow, ncol)),
                      cmap=cmap,
                      edgecolors = 'black')

        xticks = np.array(range(ncol))
        yticks = np.array(range(nrow))
        
        yticklabels=list(string.uppercase[:nrow]); yticklabels.reverse()
        xticklabels=range(1,ncol+1)
        
        ax.set_xticks(xticks+0.5)
        ax.set_xticklabels(xticklabels)
        
        ax.set_yticks(yticks+0.5)
        ax.set_yticklabels(yticklabels)
        
        ax.set_title(title)
#
#        if len(res.keys())<nrow*ncol:
#            if missing_col !=[]:
#                exp=np.reshape(range(1, (ncol-len(missing_col))*nrow+1), (nrow, ncol-len(missing_col)))
#                exp=(np.hstack((np.zeros(shape=(self.nb_row, len(missing_col))), exp))).astype(int)
#            else:
#                raise ValueError
#        else:
#            exp = range(1, nrow*ncol+1)
#            exp=np.reshape(exp, (nrow, ncol))
        #exp=np.flipud(exp)
        for well in filter(lambda x: x!=-1, res.keys()):
            a,b=np.where(where_wells==well)
            for x,y in zip(a,b):
                #restricting the name to max 6 char otherwise it's not lisible
                #txt = res[well]['Name'][:6]
                txt='{}'.format(res[well]['Xenobiotic'][:4])
                txt+='\n{}'.format(res[well]['Dose'])
                ax.text(y+0.1, nrow-x-1+0.1, txt, fontsize=6)
                if well in qc_info:
                    ax.text(y+0.1, nrow-x-1+0.1, 'QC', fontweight = 'bold', fontsize=10)
                    txt='QC'
                #txt+='\n{} {}'.format(res[well]['Medium'][:6],res[well]['Serum'])
                
                texts.append((y, nrow-x-1, txt))
        if show:
            p.show()
        else:
            p.savefig(full_filename)
        p.close(fig)
        print where_wells
        return texts
        
        
    def plotArray(self, data, filename, plotDir=None, labeledPosD=None,
                  title='', label_colors = True, absent =None,
                  legend_plot = False, legend_filename=None, texts=None,
                  control_color='black', numeric_as_int = False,
                  norm = False, min_val=0, max_val=1,
                  grid=False,int_labels=False, show=False):
        
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
        #print title, np.min(data), np.max(data)
        nrow = data.shape[0]
        ncol = data.shape[1]
        fig = p.figure(figsize=(12,7))
        ax=fig.add_subplot(111)
        
        if numeric_as_int:
#TODO
            pass
        else:
            nb_colors = 256
        
        if label_colors:
            #meaning we're not dealing with numerical data
#TODO
            pass
        else:
            #the figure actually gives the number of colors that is used. Max: 9. It has nothing to do
            #with the number of rgb quantization levels, which is 256, of the corresponding matplotlib colormap
            cmap = brewer2mpl.get_map('RdPu', 'Sequential', 9).mpl_colormap
            min_ = min_val if norm else np.min(data)
            max_ = max_val if norm else np.max(data) 

            ax.pcolormesh(np.flipud(data),
                      cmap=cmap,
                      vmin = min_,
                      vmax=max_,
                      edgecolors = 'black')
            

        xticks = np.array(range(ncol))
        yticks = np.array(range(nrow))
        
        yticklabels=list(string.uppercase[:data.shape[0]]); yticklabels.reverse()
        xticklabels=range(1,data.shape[1]+1)
        
        ax.set_xticks(xticks+0.5)
        ax.set_xticklabels(xticklabels)
        
        ax.set_yticks(yticks+0.5)
        ax.set_yticklabels(yticklabels)
        
        ax.set_title(title)
    #as we need to turn around the data to plot it we should also turn around the info of wells where data is absent
        abs = np.zeros(shape = data.shape)
        abs[absent] = 1
        absent=np.where(np.flipud(abs)==1)
        for el in zip(absent[0], absent[1]):
            ax.text(el[1]+0.5, el[0]+0.5, 'm', fontweight = 'bold', fontsize=10)
            
        no_images = np.where(np.flipud(data)==-2)
        for el in zip(no_images[0], no_images[1]):
            ax.text(el[1]+0.5, el[0]+0.5, 'no im.',fontsize=8)

        if grid: ax.grid(True)
        if texts is not None:
            for el in filter(lambda x: (x[0], x[1]) not in zip(absent[1], absent[0]), texts):
                if 'QC' not in el[2]:
                    ax.text(el[0]+0.1, el[1]+0.1, el[2],fontsize=6)
                else:
                    ax.text(el[0]+0.5, el[1]+0.5, 'QC', fontweight = 'bold', fontsize=10)
        if show:
            p.show()
        else:
            p.savefig(full_filename)
        p.close(fig)
        
        if legend_plot:
            self.plotNumLegend(cmap,max_, min_, title, filename=filename, legendDir=plotDir, int_labels=int_labels)
            
            
            
            
    def plotNumLegend(self, cmap,max_, min_, title, filename, legendDir=None, int_labels=False):
        filename = 'legend_{}'.format(filename)
        if legendDir is None:
            legendDir = self.legendDir
        nb_breaks = 5
        tickL = np.array([float(x)/(nb_breaks - 1) for x in range(nb_breaks)])
        
        if int_labels:
            colBreaksL = np.array([min_+ float(x) / (nb_breaks - 1) * (max_ - min_) for x in range(nb_breaks)], dtype=int)
        else:
            colBreaksL = np.array([min_ + float(x) / (nb_breaks - 1) * (max_ - min_) for x in range(nb_breaks)])

        full_filename = os.path.join(legendDir, filename)
        
        fig = p.figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, orientation='horizontal')
        cb.set_ticks(tickL)
        cb.set_ticklabels(["{0:.2f}".format(a) for a in colBreaksL])
        ax.set_title("Legend {}".format(title))
        ax.grid(True)
        p.savefig(full_filename)
        p.close(fig)
        return
    
    
class WellPlotter():
    def __init__(self, plotDir):
        self.plotDir = plotDir
        if not os.path.isdir(self.plotDir):
            os.makedirs(self.plotDir)
        return


    def prepareData(self, dataL):
        if dataL!=[]:
            return np.array(dataL)
        else:        
            return None    
        
        
    def plotEvolution(self, well_num, filename, plotDir=None, max_=None, min_=None,
                      add_line=None, title='', ctrl_data=None, show=False, **kwargs):
        data=kwargs['data']
        labels=kwargs['labels']
        colors=kwargs['colors']
        if data is None:
            return 
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
        print title
        try:
            fig, ax = p.subplots(1, figsize=kwargs['figsize'])
        except KeyError:
            fig, ax = p.subplots(1)
        if data is not None:
            for i,data_set in enumerate(data):
                ax.scatter(range(data_set.shape[0]), data_set, label="{}".format(labels[i]), color = colors[i], s=10)
        if ctrl_data is not None:
            ax.scatter(range(data.shape[0]), ctrl_data, label="Ctrl", color = 'black', s=10)

        if min_ is not None:
            ax.set_ylim(min_, max_)
        if add_line is not None:
            p.axhline(add_line, color="red")
        ax.set_title(title)
        ax.legend()
        ax.grid(True)
    #as we need to turn around the data to plot it we should also turn around the info of wells where data is absent
        
        if show:
            p.show()
        else:
            p.savefig(full_filename)
        p.close(fig)
            
        return