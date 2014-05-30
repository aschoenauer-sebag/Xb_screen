import string, getpass, datetime
from warnings import warn
from collections import defaultdict

from scipy.stats import nanmean, nanstd

if getpass.getuser()=='aschoenauer':
    import matplotlib as mpl
    mpl.use('Agg')
    show = False
elif getpass.getuser()=='lalil0u':
    import matplotlib as mpl
    show=True

from util import settings
from tracking.importPack.imp import importTargetedFromHDF5

from analyzer import *
from plateSetting import readPlateSetting

from util.plots import makeColorRamp, couleurs, basic_colors
from util.movies import makeMovieMultiChannels
import brewer2mpl
import matplotlib.pyplot as p

from plates.models import Well


class HTMLGenerator():
    def __init__(self, settings_file, nb_row = 8, nb_col =12):
        self.settings = settings.Settings(settings_file, globals())
        
        self.plate = self.settings.plate
        self.nb_row = nb_row
        self.nb_col = nb_col
        self.ap = ArrayPlotter(plotDir=self.settings.plot_dir, 
                               legendDir=self.settings.plot_dir, 
                               nb_row=self.nb_row, 
                               nb_col=self.nb_col)
        self.wp = WellPlotter(plotDir=self.settings.plot_dir)
        
        
    def targetedDataExtraction(self, plateL, featureL):
        newFrameLot = None
        dataFolder = self.settings.raw_result_dir  
        print "Looking for features ", featureL
        for plate in plateL:
            listW = filter(lambda x: '.hdf5' in x, os.listdir(os.path.join(dataFolder, plate, 'hdf5')))
            for filename in listW:
                well=filename[:-5]
                    
                filename = os.path.join(dataFolder, plate,"hdf5", filename)
                
                try:
                    featureL, frameLotC= importTargetedFromHDF5(filename, plate, well,featureL, secondary=self.settings.secondaryChannel)
                except ValueError:
                    sys.stderr.write(sys.exc_info())
                    sys.stderr.write("File {} containing data for plate {}, well {} does not contain all necessary data".format(filename, plate, well))

                if newFrameLot == None:
                    newFrameLot = frameLotC 
                else: newFrameLot.addFrameLot(frameLotC)
        return featureL, newFrameLot 
    
    def featureComparison(self, table, index):
        primary = table[:,index]; secondary =table[:,self.FEATURE_NUMBER+index]+0.00001
        return np.divide(primary,secondary)
    
    def featureHist(self, array, bins, binOfInterest):
        h1,_ = np.histogram(array, bins = bins)
        return h1[binOfInterest]/float(np.sum(h1))
        
    def formatData(self, frameLot, resD, featureL, featureChannels):
        featureL = list(featureL)
        for plate in frameLot.lstFrames:
            if plate not in resD:
                print plate, resD.keys()#just checkin
                warn("Plate {} not in result dictionary from CSV file extraction".format(plate))
                continue
            for well in frameLot.lstFrames[plate]:
                if int(well[:-3]) not in resD[plate]:
                    warn("Well {} not in result dictionary from CSV file extraction".format(well))
                    continue
                
                #well2=well[:-3]
                result={'cell_count':[], 'red_only':[], 'circularity':[]}
                argL = zip(filter(lambda x: x != 'roisize', featureL), featureChannels)
                result.update({'{}_ch{}'.format(arg[0], arg[1]+1):[] for arg in argL})
                
                for frame_nb in frameLot.lstFrames[plate][well]:
                    frame = frameLot.lstFrames[plate][well][frame_nb]
                    result["cell_count"].append(frame.centers.shape[0])
                    red_only_dist =self.featureComparison(frame.features, featureL.index("roisize"))
                    result["red_only"].append(self.featureHist(red_only_dist, bins = [0, 0.95, 1], binOfInterest = 1))

                    for arg in argL:
                        if arg[0]=='irregularity':
                            result['{}_ch{}'.format(arg[0], arg[1]+1)].append(frame.features[:,arg[1]*self.FEATURE_NUMBER+featureL.index(arg[0])]+1)
                            continue
                        result['{}_ch{}'.format(arg[0], arg[1]+1)].append(frame.features[:,arg[1]*self.FEATURE_NUMBER+featureL.index(arg[0])])
                    if "circularity_ch2" in result:
                    #circularity at well level
                        result['circularity'].append(self.featureHist(result["circularity_ch2"][-1], bins = [1, 1.5, 5], binOfInterest = 0))

                        
                result["initCellCount"]=result["cell_count"][0]
                result["endCellCount"]=result["cell_count"][-1]
                result["proliferation"]=result["cell_count"][-1]/float(result["cell_count"][0])

                result["initNucleusOnly"]=result["red_only"][0]
                result["endNucleusOnly"]=result["red_only"][-1]
                
                result['initCircularity']= self.featureHist(result["circularity_ch2"][0], bins = [1, 1.4, 5], binOfInterest = 0)
                result['endCircularity'] = self.featureHist(result["circularity_ch2"][-1], bins = [1, 1.4, 5], binOfInterest = 0)
                
                result['initCircMNucleus']= result['initCircularity']-result["initNucleusOnly"]
                result['endCircMNucleus']= result['endCircularity']-result["endNucleusOnly"]
                
                result["death"]=result['endCircMNucleus']/float(result['initCircMNucleus']+0.0001)
                resD[plate][int(well[:-3])].update(result)
                        
        return 1
    
    def addDummyExperiments(self, labtekId, resD):
        '''
        For experiments for which we don't have results, add them with values to zero
        '''
        nb_exp = self.nb_row * self.nb_col
        idL = range(nb_exp) if self.settings.startAtZero else range(1, nb_exp+1)#['%s--%03i' % (labtekId, (i+1)) for i in range(nb_exp)]
        for id in idL:
            if not resD.has_key(id):
                resD[id] = {'drug': 'XXX',
                            'concentration': '0um',
                            'initCellCount': np.nan,
                            'endCellCount': np.nan,
                            'proliferation': np.nan                            
                            }
            
        return  
    
    def saveResults(self, plate, resD):
        try:
            resCour = resD[plate]
        except KeyError:
            print plate, ' not in result dictionary.'
            return
        
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
                                    ax.scatter(range(len(resCour[well]['cell_count'])), val, 
                                        color=couleurs[paramIntVal.index(legendName)],
                                        label = legendName)
                                else:
                                    ax.scatter(range(len(resCour[well]['cell_count'])), val,
                                        color=couleurs[paramIntVal.index(legendName)])
                                ax.text(len(resCour[well]["cell_count"]), val[-1], 
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
            

    
    def generatePlatePlots(self, plate, resD):
        
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
            {'min': self.settings.density_plot_settings['min_circularity'], 'max': self.settings.density_plot_settings['max_circularity'], 'int_labels': False}
            ]
        
        plotDir = os.path.join(self.settings.plot_dir, plate)
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)
        try:
            missing_col = self.settings.missing_cols[plate]
        except KeyError:
            missing_col = []
        
        ###printing plate setup
        print "working on plate setup"
        filename = '%s--%s.png' % ('plate', plate)
        well_setup, texts = self.ap.plotPlate(resCour,self.params, filename, plotDir,missing_col = missing_col, title='{} {}'.format(plate, 'plate'))
        
        for pheno, setting in zip(['initCellCount', 'endCellCount', 'proliferation', 
                                   'initCircularity', 'endCircularity',  
                                   "initNucleusOnly", "endNucleusOnly",
                                   'initCircMNucleus', 'endCircMNucleus'], settingL):
            print "working on ", pheno
            filename = '%s--%s.png' % (pheno, plate)
            
            data = self.ap.prepareData(resCour, self.ap.getPhenoVal(pheno), missing_col)
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
#TODO investigate this normalization story
        return well_setup
    
    def generateWellPlots(self, plate, resD, well_setup):
        try:
            resCour = resD[plate]
        except KeyError:
            print plate, ' not in result dictionary.'
            return
        well_setup=well_setup.flatten()
        plotDir = os.path.join(self.settings.plot_dir, plate)
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        for well in resCour:
            #we want to plot cell count and circularity evolutions
            
            print "working on well ", well
            final_well_num = np.where(well_setup==well)[0][0]+1
            for pheno in self.settings.well_features:
                print "working on ", pheno
                filename = '{}_{}--W{:>05}.png'.format(pheno,plate, final_well_num)
                try:
                    data = self.wp.prepareData(resCour[well][pheno])
                except KeyError:
                    continue
                else:
                    if pheno =='circularity':
                        data2=self.wp.prepareData(resCour[well]['red_only'])
                    else:
                        data2=None
                    try:
                        min_=self.settings.well_plot_settings[pheno][0];max_=self.settings.well_plot_settings[pheno][1]
                    except KeyError:
                        max_=None; min_=None
                        
                    self.wp.plotEvolution(final_well_num, data, filename,pheno, plotDir=plotDir,
                    max_=max_, min_=min_, data2=data2,
                    title='Plate {}, well {}, evolution of {}'.format(plate, final_well_num, pheno))
        return
    
    def generateMovies(self, plate, well_setup):
        well_setup=well_setup.flatten()

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
                    makeMovieMultiChannels(imgDir=imgDir, outDir=self.settings.movie_dir, plate=plate, well=np.where(well_setup==well)[0][0]+1, 
                                           redo=self.settings.redoMovies)
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
        
        
    def __call__(self, plateL=None, featureL = None, featureChannels =None, saveDB=True, fillDBOnly=False):
        
    #Getting plate list
        if plateL is None:
            plateL = [self.settings.plate] if type(self.settings.plate)!=list else self.settings.plate
    #Getting features of interest
        if featureL is None:
            featureL = self.settings.featuresOfInterest
            featureChannels = self.settings.featureChannels
        assert(len(featureL)==len(featureChannels)), "Not same number of features of interest and channels to consider them"
        self.FEATURE_NUMBER = len(featureL)
        
    #Adding roisize to be able to assess number of cells that only have fluorescent nucleus
        if 'roisize' not in featureL:
            featureL.append('roisize')
            self.FEATURE_NUMBER +=1
            
        #first we get the plate setting, not counting empty wells. In particular it means that 
        #wells are counted according to their Zeiss number and not number on the plate.
        #After that step, resD contains only information regarding well treatments
        print ' *** reading plate settings and saving info for plate, well, condition and treatment in db for %s ***' % plateL
        
        #idL is a dictionary {plate:{well_Zeiss_num:well_id in db}} to renumber the wells in case
        #some wells were on the plate were not taken into account while doing plate setup on Zeiss software.
        #in resD, wells are numbered according to Zeiss numbering
        resD, self.params, idL = readPlateSetting(plateL, self.settings.confDir, self.nb_row, self.nb_col, 
                                             addPlateWellsToDB=saveDB,
                                             countEmpty = self.settings.countEmpty,
                                             startAtZero = self.settings.startAtZero,
                                             plateName = self.settings.name,
                                             dateFormat= self.settings.date_format,
                                             )

        if not fillDBOnly:
            print ' *** get result dictionary ***'
            featureL, frameLot = self.targetedDataExtraction(plateL, featureL)
            self.formatData(frameLot, resD, featureL, featureChannels)
            for plate in plateL:
                if True:
                    print ' *** generate density plots for %s: ***' % plate
                    #well_setup is an array representing the plate with the Zeiss well numbers
                    #on the plate. So well_setup.flatten() gives the following: on absolute position
                    #final well_number there is either a zero either the Zeiss well number.
                    
                    #This function is used to generate all 
                    well_setup = self.generatePlatePlots(plate, resD)
                    
                    print ' *** generate well plots ***'
                    self.generateWellPlots(plate, resD, well_setup)
    
                    print ' *** saving result dictionary ***'
                    self.saveResults(plate, resD)
    
                    if len(np.where(well_setup>0)[0])>1:
                        print ' *** changing well numbers in db, plate ', plate
                        self.changeDBWellNumbers(plate, well_setup, idL)
    
                    print ' *** generate movies ***'
                    self.generateMovies(plate, well_setup)
         
                else: 
                    print ' ERROR while working for %s' % plate
                    continue
        return
    
    
    
class ArrayPlotter():
    def __init__(self, plotDir, legendDir, nb_row, nb_col):
        self.plotDir = plotDir
        if not os.path.isdir(self.plotDir):
            os.makedirs(self.plotDir)
            
        self.legendDir = legendDir
        self.nb_row = nb_row
        self.nb_col = nb_col
        return
    
    class getPhenoVal(object):

        def __init__(self, phenoClass):
            self.phenoClass = phenoClass

        def __call__(self, valuesD):
            try:
                return valuesD[self.phenoClass]
            except KeyError:
                return -1


    def prepareData(self, resD, getLabel,missing_col, expL=None):
        if expL ==None:
            expL = resD.keys()
#TODO faire une correspondance entre les numeros des puits suivant si countEmpty ou pas

        nb_exp = len(expL)
        hitValueL = [0 for x in range(nb_exp)]
    
        for i in range(len(expL)):
            id = expL[i]
            if resD.has_key(id):
                label = getLabel(resD[id])
            else:
                label = -1
            hitValueL[i] = label
            
        try:
            hitData=np.reshape(hitValueL, (self.nb_row, self.nb_col))
        except ValueError:
            warn("Missing columns in Zeiss plate setup. We put zeros as columns "+str(missing_col))
            if missing_col !=[]:
                hitData=np.reshape(hitValueL, (self.nb_row, self.nb_col-len(missing_col)))
                miss = np.zeros(shape=(self.nb_row, len(missing_col)))
                miss.fill(-1)
                exp = np.hstack((miss, hitData))
                print np.min(exp[np.where(exp>-1)]), np.max(exp)
                return exp
            else:
                raise ValueError
        else:   
            print np.min(hitData[np.where(hitData>-1)]), np.max(hitData)
            print hitData
            return hitData    
  
    def plotPlate(self, res,params, filename, plotDir=None,missing_col=[], title='', show=False):
        
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
        texts = []
        nrow = self.nb_row
        ncol = self.nb_col
        fig = p.figure(figsize=(12,7))
        ax=fig.add_subplot(111)
        cmap = brewer2mpl.get_map('RdPu', 'Sequential', 9).mpl_colormap

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

        if len(res.keys())<nrow*ncol:
            if missing_col !=[]:
                exp=np.reshape(range(1, (ncol-len(missing_col))*nrow+1), (nrow, ncol-len(missing_col)))
                exp=(np.hstack((np.zeros(shape=(self.nb_row, len(missing_col))), exp))).astype(int)
            else:
                raise ValueError
        else:
            exp = range(1, nrow*ncol+1)
            exp=np.reshape(exp, (nrow, ncol))
        #exp=np.flipud(exp)
        for well in res.keys():
            a,b=np.where(exp==well)
            txt = res[well]['Name']
            txt+='\n'+res[well]['Xenobiotic']+' '+res[well]['Dose']
            txt+='\n'+res[well]['Medium']+' '+res[well]['Serum']
            ax.text(b+0.1, nrow-a-1+0.1, txt, fontsize=10)
            texts.append((b, nrow-a-1, txt))
        if show:
            p.show()
        else:
            p.savefig(full_filename)
        print exp
        return exp, texts
        
        
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
        nrow = self.nb_row
        ncol = self.nb_col
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
            ax.text(el[1]+0.5, el[0]+0.5, 'n', fontweight = 'bold', fontsize=20)
        if grid: ax.grid(True)
        if texts is not None:
            for el in filter(lambda x: (x[0], x[1]) not in zip(absent[1], absent[0]), texts):
                ax.text(el[0]+0.1, el[1]+0.1, el[2],fontsize=10)
        if show:
            p.show()
        else:
            p.savefig(full_filename)
        
        if legend_plot:
            self.plotNumLegend(cmap,max_, min_, title, filename=filename, legendDir=plotDir, int_labels=int_labels)
            
            
            
            
    def plotNumLegend(self, cmap,max_, min_, title, filename, legendDir=None, int_labels=False):
        filename = 'legend_{}'.format(filename)
        if legendDir is None:
            legendDir = self.legendDir
        nb_breaks = 5
        tickL = np.array([float(x)/(nb_breaks - 1) for x in range(nb_breaks)])
        
        if int_labels:
            colBreaksL = np.array([float(x) / (nb_breaks - 1) * (max_ - min_) for x in range(nb_breaks)], dtype=int)
        else:
            colBreaksL = np.array([float(x) / (nb_breaks - 1) * (max_ - min_) for x in range(nb_breaks)])

        full_filename = os.path.join(legendDir, filename)
        
        fig = p.figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, orientation='horizontal')
        cb.set_ticks(tickL)
        cb.set_ticklabels(["{0:.2f}".format(a) for a in colBreaksL])
        ax.set_title("Legend {}".format(title))
        ax.grid(True)
        p.savefig(full_filename)
   
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
        
        
    def plotEvolution(self, well_num, data, filename,pheno, plotDir=None, max_=None, min_=None,data2=None,
                  title='', ctrl_data=None, show=False):
        
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
        print title

        fig, ax = p.subplots(1)
        if data is not None:
            ax.scatter(range(data.shape[0]), data, label="{}".format(pheno), color = 'blue')
        if ctrl_data is not None:
            ax.scatter(range(data.shape[0]), ctrl_data, label="Ctrl", color = 'black')
        if data2 is not None:
            ax.scatter(range(data.shape[0]), data2, label="Fraction of cells simply marked", color = 'red')
        if min_ is not None:
            ax.set_ylim(min_, max_)
        ax.set_title(title)
        ax.legend()
        ax.grid(True)
    #as we need to turn around the data to plot it we should also turn around the info of wells where data is absent
        
        if show:
            p.show()
        else:
            p.savefig(full_filename)
            
        return