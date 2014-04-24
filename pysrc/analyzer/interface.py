import string, getpass, datetime
from warnings import warn

if getpass.getuser()=='aschoenauer':
    import matplotlib as mpl
    mpl.use('Agg')
elif getpass.getuser()=='lalil0u':
    import matplotlib as mpl

from util import settings
from tracking.importPack.imp import importTargetedFromHDF5

from analyzer import *
from plateSetting import readPlateSetting

from util.plots import makeColorRamp
from util.movies import makeMovieMultiChannels
import brewer2mpl
import matplotlib.pyplot as p


class HTMLGenerator():
    def __init__(self, settings_file, nb_row = 8, nb_col =12):
        self.settings = settings.Settings(settings_file, globals())
        
        self.plate = self.settings.plate
        self.nb_row = nb_row
        self.nb_col = nb_col
        self.ap = ArrayPlotter(plotDir=os.path.join(self.settings.plot_dir, 'plots'), 
                               legendDir=os.path.join(self.settings.plot_dir, 'plots'), 
                               nb_row=self.nb_row, 
                               nb_col=self.nb_col)
        self.wp = WellPlotter(plotDir=os.path.join(self.settings.plot_dir, 'plots'))
        
        
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
                    sys.stderr.write(sys.exc_info()[1])
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
        return h1[0]/float(np.sum(h1))
        
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
                result={'cell_count':[], 'red_only_dist':[], 'circularity':[]}
                argL = zip(filter(lambda x: x != 'roisize', featureL), featureChannels)
                result.update({'{}_ch{}'.format(arg[0], arg[1]+1):[] for arg in argL})
                
                for frame_nb in frameLot.lstFrames[plate][well]:
                    frame = frameLot.lstFrames[plate][well][frame_nb]
                    result["cell_count"].append(frame.centers.shape[0])
                    result["red_only_dist"].append(self.featureComparison(frame.features, featureL.index("roisize")))

                    for arg in argL:
                        if arg[0]=='irregularity':
                            result['{}_ch{}'.format(arg[0], arg[1]+1)].append(frame.features[:,arg[1]*self.FEATURE_NUMBER+featureL.index(arg[0])]+1)
                            continue
                        result['{}_ch{}'.format(arg[0], arg[1]+1)].append(frame.features[:,arg[1]*self.FEATURE_NUMBER+featureL.index(arg[0])])
                    if "circularity_ch2" in result:
                    #circularity at well level
                        result['circularity'].append(self.featureHist(result["circularity_ch2"][-1], bins = [1, 1.4, 5], binOfInterest = 0))

                        
                result["initCellCount"]=result["cell_count"][0]
                result["endCellCount"]=result["cell_count"][-1]
                result["proliferation"]=result["cell_count"][-1]/float(result["cell_count"][0])
                
                result['initCircularity']= self.featureHist(result["circularity_ch2"][0], bins = [1, 1.4, 5], binOfInterest = 0)
                result['endCircularity'] = self.featureHist(result["circularity_ch2"][-1], bins = [1, 1.4, 5], binOfInterest = 0)
                
                result["death"]=result['initCircularity']/float(result['endCircularity'])
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

    
    def generateDensityPlots(self, plate, resD):
        
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
            {'min': self.settings.density_plot_settings['min_death'], 'max': self.settings.density_plot_settings['max_death'], 'int_labels': False}
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
        self.ap.plotPlate(resCour,self.params, filename, plotDir, title='{} {}'.format(plate, 'plate'))
        
        for pheno, setting in zip(['initCellCount', 'endCellCount', 'proliferation', 'initCircularity', 'endCircularity', 'death'], settingL):
            print "working on ", pheno
            filename = '%s--%s.png' % (pheno, plate)
            
            data = self.ap.prepareData(resCour, self.ap.getPhenoVal(pheno), missing_col)
            absent = np.where(data==0)
            self.ap.plotArray(data, filename, plotDir=plotDir,
                title='{} {}'.format(plate, pheno),
                absent= absent,
                int_labels = setting['int_labels'],
                norm = True,
                min_val=setting['min'], max_val=setting['max'], 
                label_colors = False, #labeledPosD=labeledPosD,
                legend_plot = True, legend_filename='legend_%s' % filename)
#TODO investigate this normalization story
        return absent
    
    def generateAllPlots(self, plate, resD):
        try:
            resCour = resD[plate]
        except KeyError:
            print plate, ' not in result dictionary.'
            return

        plotDir = os.path.join(self.settings.plot_dir, plate)
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        for well in resCour:
            #we want to plot cell count and circularity evolutions
            
            print "working on well ", well
            
            for pheno in self.settings.well_features:
                print "working on ", pheno
                filename = '{}_{}--W{:>05}.png'.format(pheno,plate, well)
                try:
                    data = self.wp.prepareData(resCour[well][pheno])
                except KeyError:
                    continue
                else:
                    self.wp.plotEvolution(well, data, filename, plotDir=plotDir,
                    title='Plate {}, well {}, evolution of {}'.format(plate, well, pheno))
        return
    
    def generateMovies(self, plate, wellL=None):
        if wellL ==None:
            wellL = range(1, self.nb_col*self.nb_row+1) if not self.settings.startAtZero else range(self.nb_col*self.nb_row)
            
        for well in wellL:
            imgDir= os.path.join(self.settings.raw_data_dir, plate, 'W{:>05}'.format(well))
            try:
                l=os.listdir(imgDir)
            except OSError:
                print "No image for well ", well, 'plate ', plate
                continue
            else:
                if l==[]:
                    print "No image for well ", well, 'plate ', plate
                    continue
                else:
                    makeMovieMultiChannels(imgDir=imgDir, outDir=self.settings.movie_dir, plate=plate, well=well, redo=self.settings.redoMovies)
        return        
        
        
    def __call__(self, plateL=None, featureL = None, featureChannels =None):
        
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
        #first we get the plate setting, counting empty wells.
        #After that step, resD contains only information regarding well treatments
        print ' *** reading plate settings and saving info for plate, well, condition and treatment in db for %s ***' % plateL
        
        #la liste d'id sert dans le cas ou des puits sont vides pour renumeroter eventuellement les puits
        resD, self.params, idL = readPlateSetting(plateL, self.settings.confDir, self.nb_row, self.nb_col, 
                                             countEmpty = self.settings.countEmpty,
                                             startAtZero = self.settings.startAtZero,
                                             plateName = self.settings.name,
                                             dateFormat= self.settings.date_format,
                                             addPlateWellsToDB=True,
                                             addCondToDB = True, addTreatmentToDB=True)
        pdb.set_trace()
        print ' *** get result dictionary ***'
        featureL, frameLot = self.targetedDataExtraction(plateL, featureL)
        self.formatData(frameLot, resD, featureL, featureChannels)
        
        #for labtekId in labtekIdL:
        #   self.copyMovies(labtekId)
#todo REGLER LE PBL DES NUMEROS DE PUITS ET ORDRE DES PUITS SUR LA PLAQUE
        for plate in plateL:
            if True:
                print ' *** generate density plots for %s: ***' % plate
                absent = self.generateDensityPlots(plate, resD)
                
                print ' *** generate spot plots ***'
                self.generateAllPlots(plate, resD)

#                for control in self.settings.controlD.keys():
#                    print ' *** generate %s plots for %s' % (control, plate)
#                    #self.generateControlPlots(plate, control)
                print ' *** generate movies ***'
                self.generateMovies(plate)
#TODO : mettre les data directement dans le dossier pour les mettre sur l'interface ?                
                #self.copyMovies(plate)
#                
#                print ' *** adding plate and wells to databases, plate ' % plate
#                self.addToDB(plate, resD, absent)        
            else: 
                print ' ERROR while generating page for %s' % plate
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
                return 0


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
                label = 0
            hitValueL[i] = label
            
        try:
            hitData=np.reshape(hitValueL, (self.nb_row, self.nb_col))
        except ValueError:
            warn("Missing columns in Zeiss plate setup. We put zeros as columns "+str(missing_col))
            if missing_col !=[]:
                hitData=np.reshape(hitValueL, (self.nb_row, self.nb_col-len(missing_col)))
                exp = np.hstack((np.zeros(shape=(self.nb_row, len(missing_col))), hitData))
                return exp
            else:
                raise ValueError
        else:        
            return hitData    
  
    def plotPlate(self, res,params, filename, plotDir=None, title='', show=False):
        
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)

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

        exp = range(1, nrow*ncol+1)
        exp=np.reshape(exp, (nrow, ncol))
        #exp=np.flipud(exp)
        for well in res.keys():
            a,b=np.where(exp==well)
            txt = res[well]['Name']
            for param in params: txt+='\n'+res[well][param]
            ax.text(b+0.1, nrow-a-1+0.1, txt, fontsize=10)
        
        if show:
            p.show()
        else:
            p.savefig(full_filename)
        
        
    def plotArray(self, data, filename, plotDir=None, labeledPosD=None,
                  title='', label_colors = True, absent =None,
                  legend_plot = False, legend_filename=None, 
                  control_color='black', numeric_as_int = False,
                  norm = False, min_val=0, max_val=1,
                  grid=False,int_labels=False, show=False):
        
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
        print title, np.min(data), np.max(data)
        nrow = self.nb_row
        ncol = self.nb_col
        fig, ax = p.subplots(1)
        
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
            ax.text(el[1]+0.5, el[0]+0.5, 'n', fontweight = 'bold')
        if grid: ax.grid(True)
        
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
        
        
    def plotEvolution(self, well_num, data, filename, plotDir=None,
                  title='', ctrl_data=None, show=False):
        
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
        print title

        fig, ax = p.subplots(1)
        if data is not None:
            ax.scatter(range(data.shape[0]), data, label="Well {}".format(well_num), color = 'blue')
        if ctrl_data is not None:
            ax.scatter(range(data.shape[0]), ctrl_data, label="Ctrl", color = 'black')
        ax.set_title(title)
        ax.legend()
        ax.grid(True)
    #as we need to turn around the data to plot it we should also turn around the info of wells where data is absent
        
        if show:
            p.show()
        else:
            p.savefig(full_filename)
            
        return