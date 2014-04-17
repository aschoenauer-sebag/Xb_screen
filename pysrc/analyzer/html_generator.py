import pdb
#import classification
#import spotAnalyzer
from rpy2 import robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
RDevice = importr('grDevices')
import copy
#import plot_utilities
import numpy as np
import operator

# import colors from r (for heat maps)
#import colors
from util.plots import makeColorRamp
import os
import pickle
from optparse import OptionParser

from util.settings import Settings
import string
    
from scipy import stats
import shutil

import re

class HTMLGenerator(object):
    def __init__(self, settings_filename):
        
        self.settings =Settings(settings_filename, globals())
                
        self.html_dir = self.settings.html_dir
        self.baseDir = self.settings.result_dir
        self.ap = ArrayPlotter(plotDir=os.path.join(self.html_dir, 'plots'), 
                               legendDir=os.path.join(self.html_dir, 'plots'), 
                               nb_row=self.settings.nb_row, 
                               nb_col=self.settings.nb_col)
        
        return
        
    
    def readLabtekConfiguration(self):
        file = open(os.path.join(self.settings.confDir, self.settings.confFilename), 'r')
        labtekConfD = pickle.load(file)
        file.close()        
        return labtekConfD

    def getExpConditionDict(self, labtekId):
        labtekConfD = self.readLabtekConfiguration()        
        if labtekConfD.has_key(labtekId):
            spotConfD = labtekConfD[labtekId]
        else:
            labtekBase = labtekId.split('_')[0]
            spotConfD = labtekConfD[labtekBase]

        key_level1 = self.settings.key_level1
        key_level2 = self.settings.key_level2
        expConditionD = {}
        posL = spotConfD.keys()
        posL.sort()
        for pos in posL:
            id = '%s--%s' % (labtekId, pos)
            k1, k2 = spotConfD[pos][key_level1], spotConfD[pos][key_level2]
            if not expConditionD.has_key(k1):
                expConditionD[k1] = {}
            if not expConditionD[k1].has_key(k2):
                expConditionD[k1][k2] = []
            expConditionD[k1][k2].append(id)
        return expConditionD
        
    # plot proliferation
    def plotProliferationSingleExperiment(self, countL, fitCountL,  
                          filename, directory, 
                          scrambledResD = None,
                          title = 'Proliferation for experiment',  
                          xlab='time (frames)', ylab='Normalized cell count', 
                          plotType='png', plotHeight=None, plotWidth=None, plotdevice=None,
                          ylim=(0.0,5.0)):


        if plotHeight is None:
            if plotType == 'png': 
                plotHeight = 600
            elif plotType == 'pdf': 
                plotHeight = 6
        
        if plotWidth is None:
            if plotType == 'png':
                plotWidth = 900
            elif plotType == 'pdf':
                plotWidth = 9
        
        maxTime = min(len(countL), len(fitCountL), len(scrambledResD['proliferation_mean']))
        timeVec = range(maxTime)
        #timeVec = [18 + 0.5 * x for x in range(maxTime)]
        
        ps = plot_utilities.PlotSettings(filename=filename, directory=directory, plotType=plotType, 
                                         main=title,
                                         plotWidth=plotWidth, plotHeight=plotHeight)
        if plotdevice is None:
            rdev = ps.getDevice()
        else:
            rdev = plotdevice

        o = r.par(mgp=(1.6, 0.4, 0), mar=(3, 3, 2, 8))

        # plot relative counts and the fit
        r.plot(timeVec, countL[:maxTime],
               type='o',
               col='blue',
               main=title, ylim=ylim,
               xlab=xlab, ylab=ylab, pch=20, lwd=1, lty = 1, 
               cex=1.0, cex_lab=1.0, cex_main=1.0, cex_axis=0.8) 
        
        r.lines(timeVec, fitCountL[:maxTime],
                col='red',
                lty = 1,
                lwd=3)

        # plot scrambled
        upperL = [x + y for (x,y) in zip(scrambledResD['proliferation_mean'],scrambledResD['proliferation_stdev']) ]
        lowerL = [x - y for (x,y) in zip(scrambledResD['proliferation_mean'],scrambledResD['proliferation_stdev']) ]
    
        r.segments(timeVec, upperL[:maxTime],
                   timeVec, lowerL[:maxTime],
                   col = "lightblue")
        r.lines(timeVec, upperL[:maxTime], col = "lightblue", lwd = 1, lty = 2)
        r.lines(timeVec, lowerL[:maxTime], col = "lightblue", lwd = 1, lty = 2)
        r.lines(timeVec, scrambledResD['proliferation_mean'][:maxTime], col = "lightblue", lwd = 2, lty = 1)
    
        r.grid(col="darkgray")
        #r.par(o)

        if plotdevice is None:
            rdev.close()
        
        return

    # plots a bundle of time series
    def plotBundle(self, bundleD, full_filename, colorsD=None, bundlePointsD=None, 
        legendL=None, title=None, y_max=None, xlab="time in hours after transfection"):

        if y_max is None:
            y_max = 0.4
            
        if legendL is None:
            legendL = bundleD.keys()
            legendL.sort()
            
        if title is None:
            title = 'data'            

        bundleIdL = bundleD.keys()
        bundleIdL.sort()

        if colorsD is None:            
            colorsL = r.rainbow(len(bundleIdL))
            colorsD = dict(zip(bundleIdL, colorsL))
        
        colorsL = [colorsD[x] for x in bundleIdL]
        
        time_min = min([len(bundleD[x]) for x in bundleD.keys()])
        timeVec = [0.5 * x for x in range(time_min)]

        
        #try:
        if True:
            r.png(full_filename, width=900, height=600)
            oldPar = r.par(xpd = False, mar = [x + y for (x,y) in zip(r.par()['mar'], [0,0,0,12])])
            #oldPar = r.par(xpd = False, oma = [0,0,0,8])
        
            print 'plot %s' % full_filename
            r.plot(timeVec, timeVec,
                   type='n',
                   main=title, ylim=(0, y_max),
                   xlab=xlab, ylab="Relative Cell Counts",
                   pch=20, lwd=1, lty = 1, 
                   cex=1.0, cex_lab=1.2, cex_main=1.5)
        
        
            for bundleId in bundleIdL:
                
                if not bundlePointsD is None:
                    r.points(timeVec, bundlePointsD[bundleId][:time_min],
                             col=colorsD[bundleId], pch=20,
                             lwd=1)
                    r.lines(timeVec, bundlePointsD[bundleId][:time_min],
                            col=colorsD[bundleId],
                            lwd=1, lty = 1)

                r.lines(timeVec, bundleD[bundleId][:time_min],
                        col=colorsD[bundleId],
                        lwd=3, lty = 1)

            r.grid(col="darkgrey")

            r.par(xpd=True)
            r.legend(max(timeVec) * 1.1, y_max, legend=legendL, fill=colorsL, cex=1.0, bg= 'whitesmoke')
            r.par(oldPar)
            
            r.dev_off() 
        #except:
        else:
            r.dev_off()
            print full_filename + ' has not been printed.'
            

        return

    def generateAllPlots(self, labtekId):
        plotBaseDir = os.path.join(self.html_dir, 'plots', labtekId)

        labtekConfD = self.readLabtekConfiguration()
        if labtekConfD.has_key(labtekId):
            spotConfD = labtekConfD[labtekId]
        else:
            labtekBase = labtekId.split('_')[0]
            spotConfD = labtekConfD[labtekBase]
        

        # retrieve control information
        ctrlPosL = self.settings.controlD['negative']
        controlProliferationD = {}
        controlProliferationFitD = {}
       
        control_len = None
        for pos in ctrlPosL:

            try:
                id = '%s--%s' % (labtekId, pos)
                spotRes = spotAnalyzer.spotAnalyzer(labtekId, pos, baseDir=self.baseDir, 
                    spotConfD=spotConfD, prediction_suffix=self.settings.prediction_suffix)
                phenoSpotD, countL = spotRes.readPrediction()
    
                if control_len is None:
                    control_len = len(countL)
                if countL[0] > 0:    
                    controlProliferationD[pos] = [float(countL[i])/float(countL[0]) for i in range(len(countL))]
                else:
                    controlProliferationD[pos] = [0.0 for i in range(len(countL))]
                controlProliferationFitD[pos] = spotRes.fitCellCount(controlProliferationD[pos])
                control_len = min(control_len, len(countL))                    
            except:
                print 'skipping control %s' % id
                continue
        
        scrambledResD = {
        'proliferation_mean': [stats.mean([controlProliferationD[pos][i] for pos in controlProliferationD.keys()]) for i in range(control_len)],
        'proliferation_stdev': [stats.stdev([controlProliferationD[pos][i] for pos in controlProliferationD.keys()]) for i in range(control_len)]
        }


        # plotting ...
        posL = spotConfD.keys()
        posL.sort()
        for pos in posL:
            try:
                spotRes = spotAnalyzer.spotAnalyzer(labtekId, pos, baseDir=self.baseDir, 
                    spotConfD=spotConfD, prediction_suffix=self.settings.prediction_suffix)
                phenoSpotD, countL = spotRes.readPrediction()
    
    
                if countL[0] > 0:    
                    prol = [float(countL[i])/float(countL[0]) for i in range(len(countL))]
                else:
                    prol = [0.0 for i in range(len(countL))]
                prolFit = spotRes.fitCellCount(prol)
                
                # calculate the data for all pheno classes
                normD = spotRes.normalizePrediction(phenoSpotD, countL)
                fitPredD = spotRes.fitPrediction(normD)
    
                # ##########
                # HERE I AM
                # plotting
                #pheno_file = filter(lambda x: x[0:len('Pheno')] == 'Pheno', os.listdir(current_plot_dir))[0]
                #prol_file = filter(lambda x: x[0:len('Proliferation')] == 'Proliferation', os.listdir(current_plot_dir))[0]
    
                plotDir = os.path.join(plotBaseDir, pos)
                if not os.path.isdir(plotDir):
                    os.makedirs(plotDir)
                    
                full_filename = os.path.join(plotDir, 'PhenoKinetics--%s--%s.png' % (labtekId, pos))            
                self.plotBundle(fitPredD, full_filename, colorsD = self.settings.COLORD, 
                    bundlePointsD = normD,
                    title='Phenotypic kinetics for %s %s' % (labtekId, pos), y_max=0.8,
                    xlab='Time (Frames)')
        
                # proliferation
                filename = 'Proliferation--%s--%s.png' % (labtekId, pos)
                self.plotProliferationSingleExperiment(countL, prolFit, filename, plotDir, 
                    scrambledResD, 
                    title='Proliferation for %s %s' % (labtekId, pos))
                
            except:
                print 'a problem occurred while plotting results for %s %s' % (labtekId, pos)
                continue
        return

    def generateControlPlots(self, labtekId, control):
        plotDir = os.path.join(self.html_dir, 'plots', labtekId, 'controls')
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)
        
        ctrlPosL = self.settings.controlD[control]
            
        controlDataD = {}
        controlDataRawD = {}
        controlProliferationD = {}
        controlProliferationFitD = {}

    
        labtekConfD = self.readLabtekConfiguration()
        if labtekConfD.has_key(labtekId):
            spotConfD = labtekConfD[labtekId]
        else:
            labtekBase = labtekId.split('_')[0]
            spotConfD = labtekConfD[labtekBase]
            
        #spotConfD = LabtekConfiguration().getSpotConfForLabtek(labtekId)
        
        control_len = None
        for pos in ctrlPosL:

            try:
                id = '%s--%s' % (labtekId, pos)
                print id
                print self.baseDir
                spotRes = spotAnalyzer.spotAnalyzer(labtekId, pos, baseDir=self.baseDir, 
                    spotConfD=spotConfD, prediction_suffix=self.settings.prediction_suffix)
                    
                phenoSpotD, countL = spotRes.readPrediction()
    
                if control_len is None:
                    control_len = len(countL)
                if countL[0] > 0:    
                    controlProliferationD[pos] = [float(countL[i])/float(countL[0]) for i in range(len(countL))]
                else:
                    controlProliferationD[pos] = [0.0 for i in range(len(countL))]
                controlProliferationFitD[pos] = spotRes.fitCellCount(controlProliferationD[pos])
                control_len = min(control_len, len(countL))                    
    
                # calculate the data for all pheno classes
                normD = spotRes.normalizePrediction(phenoSpotD, countL)
                fitPredD = spotRes.fitPrediction(normD)
    
                for pheno in phenoSpotD.keys():
                    if not controlDataRawD.has_key(pheno):
                        controlDataRawD[pheno] = {}
                    if not controlDataD.has_key(pheno):
                        controlDataD[pheno] = {}
                    controlDataRawD[pheno][pos] = normD[pheno]
                    controlDataD[pheno][pos] = fitPredD[pheno]
            except:
                print 'no data for %s' % pos

        for pheno in ['Prometaphase', 'MetaphaseAlignment', 'Shape1', 'Shape3', 'Apoptosis']:
            try:
                full_filename = os.path.join(plotDir, control + '_' + pheno + '.png') 
                self.plotBundle(controlDataD[pheno], full_filename, colorsD = None, bundlePointsD = controlDataRawD[pheno],
                                title='%s for control %s' % (pheno, control), y_max=0.4)
    
                # proliferation
                full_filename = os.path.join(plotDir, control + '_proliferation.png') 
                self.plotBundle(controlProliferationFitD, full_filename, bundlePointsD = controlProliferationD, title='proliferation for control %s' % control, y_max=5.0)
            except:
                print 'no plot for %s' % control
        return

    def addDummyExperiments(self, labtekId, resD):
        nb_exp = self.settings.nb_row * self.settings.nb_col
        idL = ['%s--%03i' % (labtekId, (i+1)) for i in range(nb_exp)]
        for id in idL:
            if not resD.has_key(id):
                resD[id] = {'drug': 'XXX',
                            'concentration': '0um',
                            'initCellCount': 0,
                            'endCellCount': 0,
                            'proliferation': 0.0                            
                            }
            
        return
        
    def generateDensityPlots(self, labtekId, resD):

        self.addDummyExperiments(labtekId, resD)

        labeledPosD = self.settings.labeledPosD
    
        idL = filter(lambda x: x.split('--')[0] == labtekId, resD.keys())
        idL.sort()
        
        pattern = [(0,0,0), (0,0,0.5), (0,1,1), (1,1,1)]
        colVecL = colors.make_colors(pattern, 500)

        settingL = [
            {'min': self.settings.density_plot_settings['min_count'], 'max': self.settings.density_plot_settings['max_count'], 'colVecL': None},
            {'min': self.settings.density_plot_settings['min_count'], 'max': self.settings.density_plot_settings['max_count'], 'colVecL': None},
            {'min': self.settings.density_plot_settings['min_proliferation'], 'max': self.settings.density_plot_settings['max_proliferation'], 'colVecL': colVecL},
            ]
                    
        # init cell count
        for pheno, setting in zip(['initCellCount', 'endCellCount', 'proliferation'], settingL):
            plotDir = os.path.join(self.html_dir, 'plots', labtekId)
            if not os.path.isdir(plotDir):
                os.makedirs(plotDir)
            filename = '%s--%s' % (pheno, labtekId)
            data = self.ap.prepareData(resD, self.ap.getPhenoVal(pheno), idL)
            self.ap.plotArray(data, filename, plotDir=plotDir, colVecL=setting['colVecL'],
                min_val=setting['min'], max_val=setting['max'], 
                label_colors = False, labeledPosD=labeledPosD,
                legend_plot = True, legend_filename='legend_%s' % filename)
        
        return

    
    def __call__(self, labtekIdL=None):
        if labtekIdL is None:
            labtekIdL = self.settings.labtekL

        print ' *** get resD'
        file = open(self.settings.id_result_filename, 'r')
        resD = pickle.load(file)
        file.close()

        #for labtekId in labtekIdL:
        #   self.copyMovies(labtekId)

        for labtekId in labtekIdL:
            if True:
                print ' *** generate density plots for %s: ***' % labtekId
                self.generateDensityPlots(labtekId, resD)
        
                for control in self.settings.controlD.keys():
                    print ' *** generate %s plots for %s' % (control, labtekId)
                    #self.generateControlPlots(labtekId, control)
        
                print ' *** generate spot plots ***'
                #self.generateAllPlots(labtekId)
                
                #self.copyMovies(labtekId)
                
                print '*** assigning movie names ***'
                self.assignMovies(resD)
                
                print ' *** generate qc page for %s' % labtekId
                self.generateQCPage(labtekId, resD)        
            else: 
                print ' ERROR while generating page for %s' % labtekId
                continue
        return

    def makeNumeric(self, inStr):
        try:
            num_val = float(inStr.replace(',', '.'))
        except:
            num_str = ''
            for c in inStr:
                if c in '1234567890,.':
                    num_str += c
                else:
                    break
            if len(num_str) == 0:
                return 0.0
            num_val = string.atof(num_str.replace(',', '.'))

        return num_val

    def copyMovies(self, labtekId):
        
        movie_raw_dir = os.path.join(self.settings.raw_result_dir, labtekId, 'movies')
        print 'copying from ', movie_raw_dir
        movies = filter(lambda x: os.path.splitext(x)[-1]=='.avi', 
                        os.listdir(movie_raw_dir))                 
        movies_target_dir = os.path.join(self.html_dir, 'movies', labtekId)
        if not os.path.exists(movies_target_dir):
            os.makedirs(movies_target_dir)
        print '%s: %i movies to copy' % (labtekId, len(movies))
        for movie in movies:
            print 'copying %s to %s' % (os.path.join(movie_raw_dir, movie), movies_target_dir)
            shutil.copy2(os.path.join(movie_raw_dir, movie), movies_target_dir)
        return
        
#    def __oldcopyMovies(self, labtekId):
#        raw_result_dir
#        movies = filter(lambda x: os.path.splitext(x)[-1]=='.avi', 
#                        os.path.join(self.settings.movie_raw_dir, labtekId))
#        movies.sort()
#        movies_target_dir = os.path.join(self.html_dir, 'movies', labtekId)
#        if not os.path.exists(movies_target_dir):
#            os.makedirs(movies_target_dir)
#        for movie in movies:                        
#            shutil.copy2(os.path.join(movie_raw_dir, movie), movies_target_dir)
#        return
        
    def assignMovies(self, resD):

        lstLabteks = set([x.split('--')[0] for x in resD.keys()])
        
        for labtekId in lstLabteks:
            relativeMovieDir = os.path.join('.', 'movies', labtekId)
            absoluteMovieDir = os.path.join(self.html_dir, 'movies', labtekId)
            movieFileL = os.listdir(absoluteMovieDir)
    
            idL = filter(lambda x: x[0:len(labtekId)] == labtekId, resD.keys())
            idL.sort()
            nomovies = []
            for id in idL:
                lt, pos = id.split('--')
                
                if self.settings.multiple_positions:
                    regex_res = dict(zip(movieFileL, [self.settings.position_regex.search(x) for x in movieFileL]))
                    movieNames = filter(lambda x: not regex_res[x] is None
                                                  and int(regex_res[x].groupdict()['Well']) == int(pos), movieFileL)
                else:                        
                    movieNames = filter(lambda x: x[0:len(id)].lower() == id.lower(), movieFileL)

                fullNames = [os.path.join(relativeMovieDir, x) for x in movieNames]                
                resD[id]['movieNames'] = fullNames

        return
        
    def generateQCPage(self, labtekId, resD):
        #plotDir = os.path.join(self.html_dir, 'plots', labtekId, 'controls')
        
        plotDir = os.path.join(self.html_dir, 'plots', labtekId)
        relativePlotDir = os.path.join('.', 'plots', labtekId)
        
        idL = filter(lambda x: x[0:len(labtekId)] == labtekId, resD.keys())
        idL.sort()

        pheno = 'initCellCount'
        imageName = os.path.join(relativePlotDir, '%s--%s.png' % (pheno, labtekId))
        mapStr = {pheno: self.makeImageMapHTML(labtekId, resD, imageName, pheno, 0)}
        
        k = 1
        for pheno in ['endCellCount', 'proliferation']:
            imageName = os.path.join(relativePlotDir, '%s--%s.png' % (pheno, labtekId))
            mapStr[pheno] = self.makeImageMapHTML(labtekId, resD, imageName, pheno, k)
            k += 1
            
        mapD = {(os.path.join(relativePlotDir, '%s--%s.png') % ('initCellCount', labtekId)): mapStr['initCellCount'],
                (os.path.join(relativePlotDir, '%s--%s.png') % ('endCellCount', labtekId)): mapStr['endCellCount'],
                (os.path.join(relativePlotDir, '%s--%s.png') % ('proliferation', labtekId)): mapStr['proliferation']
            }
        #mapD = {(os.path.join(relativePlotDir, '%s--%s.png') % ('initCellCount', labtekId)): mapStr['initCellCount']}
        rowL = ['Overview', 'Negative']
        #rowL = ['Overview', 'Scrambled']        
        rowD = {}
#        rowD['Overview'] = {'header': ['Initial Cell Count', 'End Cell Count', 'Proliferation'],
#                            'images': [os.path.join(relativePlotDir, '%s--%s.png') % (pheno, labtekId)
#                                       for pheno in ['initCellCount', 'endCellCount', 'proliferation']]}
        rowD['Overview'] = {'header': ['Initial Cell Count', 'End Cell Count', 'Proliferation'],
                            'images': [os.path.join(relativePlotDir, '%s--%s.png') % (pheno, labtekId)
                                       for pheno in ['initCellCount', 'endCellCount', 'proliferation']],
                            'map': mapD}
#        rowD['Legend'] = {'header': None,
#                          'images': [os.path.join(relativePlotDir, 'legend_%s--%s.png') % (pheno, labtekId)
#                                     for pheno in ['initCellCount', 'proliferation']]}
        #relativePlotDir = os.path.join('..', labtekId, 'plots', 'controls')
        controlPlotDir = os.path.join(relativePlotDir, 'controls')
        rowD['Negative'] = {'header': None,
                             'images': [os.path.join(controlPlotDir, 'negative_%s.png' % pheno) for pheno in 
                                       ['proliferation', 'Apoptosis', 'Prometaphase', 'MetaphaseAlignment', 'Shape1', 'Shape3']],
                             'height': [500, 500, 500, 500, 500, 500]}
        if self.settings.controlD.has_key('empty'):
            rowD['Empty'] = {'header': None,
                             'images': [os.path.join(controlPlotDir, 'empty_%s.png' % pheno) for pheno in 
                                       ['proliferation', 'Apoptosis', 'Prometaphase', 'MetaphaseAlignment', 'Shape1', 'Shape3']],
                             'height': [500, 500, 500, 500, 500, 500]}
            rowL.append('Empty')
            
                                       
#        rowD['Incenp'] = {'header': None, #['Proliferation', 'Shape 3', 'Prometaphase', 'Apoptotic'],
#                          'images': [os.path.join(relativePlotDir, 'incenp_%s.png') % pheno
#                                     for pheno in ['proliferation', 'Shape3', 'Shape1', 'Prometaphase', 'MetaphaseAlignment', 'Apoptosis']] }
#        rowD['Scrambled'] = {'header': None, #['Proliferation', 'Shape', 'Mitotic', 'Apoptotic'],
#                             'images': [os.path.join(relativePlotDir, 'scrambled_%s.png') % pheno
#                                        for pheno in ['proliferation', 'Shape3', 'Shape1', 'Prometaphase', 'MetaphaseAlignment', 'Apoptosis']] }

        # generate one row per gene/siRNA
        if True:
            expConditionL = []
            expConditionD = self.getExpConditionDict(labtekId)
            keyL = expConditionD.keys()
            for key in keyL:
                lkL = expConditionD[key].keys()
                numL = [self.makeNumeric(lk) for lk in lkL]
                tempL = zip(lkL, numL)
                tempL.sort(key=operator.itemgetter(-1))
                lkSortL = [x[0] for x in tempL]
                for lk in lkSortL:
                    expConditionL.append((key, lk))
                    
                    
            #sirnaD = self.lc.getSirnaView(labtekId)
            #expConditionL.sort(key=operator.itemgetter(0))
            relativePlotDir = os.path.join('.', 'plots', labtekId)
            for hk, lk in expConditionL:
                
                if len(expConditionD[hk][lk]) != 1:
                    print hk, lk, expConditionD[hk][lk]
                    continue
                    
                id = expConditionD[hk][lk][0]
                lt_id, pos = id.split('--')
                current_plot_dir = os.path.join(self.html_dir, 'plots', labtekId, pos)
                if True:
                    pheno_file = filter(lambda x: x[0:len('Pheno')] == 'Pheno', os.listdir(current_plot_dir))[0]
                    prol_file = filter(lambda x: x[0:len('Proliferation')] == 'Proliferation', os.listdir(current_plot_dir))[0]
        
                    current_plot_dir = os.path.join(relativePlotDir, pos)
                    exp_cond_id = '%s (%s)' % (hk, lk)
                    
                    movieNames = resD[id]['movieNames']
                    moviePageName = self.getMoviePageName(id)
                    #moviePageName = os.path.join('.', labtekId, os.path.basename(movieName).replace('.avi', '.html'))
                    rowD[exp_cond_id] = {'header': None,
                                         'images': [os.path.join(current_plot_dir, pheno_file), os.path.join(current_plot_dir, prol_file)],
                                         'height': [500, 500],
                                         'header_link': moviePageName
                                         }
                    rowL.append(exp_cond_id)
                else:
                    print 'missing condition: ', hk, lk
                
        self.exportPlateHTML('Summary--%s.html' % labtekId, labtekId, rowD, rowL)
        return

    # HERE I AM
    #self.makeMovieLinkPage(id, movieNames, moviePageName, gene, sirna)
    def makeMovieLinkPage(self, id, movieNames, moviePageName, gene, sirna):
        labtekId, pos = id.split('--')
        movie_html_dir = os.path.join(self.html_dir, labtekId)
        if not os.path.isdir(movie_html_dir):
            os.makedirs(movie_html_dir)
            
        
        labtekId, pos = id.split('--')
        
        #filename = os.path.join(movie_html_dir, moviePageName)
        filename = os.path.join(self.html_dir, moviePageName)
        movie_html_file = open(filename, 'w')
        print 'generating ', filename

        
        tempStr = """
<html>        
  <head><title>Experiment: %s %s</title></head>
  <body>
    
    <table align="left" border="0" cellspacing="10" cellpadding="0">
""" % (labtekId, pos)

        for movieName in movieNames:
            tempStr += \
"""    <th align="center"> %s %s %s </th>
    <th align="center" ><EMBED ID="movie" SRC="%s" WIDTH="864" HEIGHT="670" CONTROLLER="TRUE" SCALE="tofit" AUTOSTART="true" LOOP=true pluginspage="http://www.apple.com/quicktime/download"></EMBED> </th> 
""" % (gene, sirna, os.path.basename(movieName), os.path.join('..', movieName))

        tempStr += \
"""    </table>
  </body>
</html>
""" 

        movie_html_file.write(tempStr)
        movie_html_file.close()
        return

    def makeSingleMovieLinkPage(self, labtekId, movieName, gene, sirna):
        movie_html_dir = os.path.join(self.html_dir, labtekId)
        if not os.path.isdir(movie_html_dir):
            os.makedirs(movie_html_dir)
            
        filename = os.path.join(movie_html_dir, os.path.basename(movieName).replace('.avi', '.html'))
        movie_html_file = open(filename, 'w')
        print 'generating ', filename

        tempStr = """
<html>        
  <head><title>Movie: %s</title></head>
  <body>
    
    <table align="left" border="0" cellspacing="10" cellpadding="0">
    <th align="center"> %s (%s) </th>
    <th align="center" ><EMBED ID="movie" SRC="%s" WIDTH="864" HEIGHT="670" CONTROLLER="TRUE" SCALE="tofit" AUTOSTART="true" LOOP=true pluginspage="http://www.apple.com/quicktime/download"></EMBED> </th> 
    </table>
  </body>
</html>
""" % (os.path.basename(movieName), gene, sirna, os.path.join('..', movieName))

        movie_html_file.write(tempStr)
        movie_html_file.close()
        return
        
    # make single 
    def makeImageMapHTML96(self, labtekId, resD, imageName, pheno='initCellCount'):        
        
        tempStr = \
"""                
<map name="map_0">
"""
        
        i=0
        for y in range(226, 2, -28):
            for x in range(24, 636, 51):                
                i += 1                
                pos = '%03i' % i
                id = '%s--%s' % (labtekId, pos)
                value = resD[id][pheno]
                movieNames = resD[id]['movieNames']
                moviePageName = self.getMoviePageName(id)
                #moviePageName = os.path.join('.', labtekId, os.path.basename(movieName).replace('.avi', '.html'))
                gene = resD[id][self.settings.gene_key]
                sirna = resD[id][self.settings.sirna_key]
                self.makeMovieLinkPage(id, movieNames, moviePageName, gene, sirna)
                #self.makeMovieLinkPage(labtekId, movieName, gene, sirna)
                
                tempStr += \
"""<area shape="rect" coords="%i, %i, %i, %i" alt="" target="_blank" 
    title="%i/96, pos: %s -- %s -- %s -- objects: %i" href="%s">
""" % (x, y-28, (x+51), y, i, pos, gene, sirna, value, moviePageName)
 
        tempStr += \
"""
</map>
"""       

        tempStr += \
"""
<td align="center"><img src="%s" border="0" usemap="#map_0"></td>
""" % imageName

        return tempStr
        
    def getMoviePageName(self, id):
        labtekId, pos = id.split('--')
        return os.path.join('.', labtekId, '%s.html' % id)

    # make single 
    def makeImageMapHTML(self, labtekId, resD, imageName, pheno='initCellCount', map_index=0):        
        
        tempStr = \
"""                
<map name="map_%i">
""" % map_index
        
        y_step = 20 # 776.0 / 32.0
        x_step = 20 # 276.0 / 12.0

        y_offset = 1
        x_offset = 22

        i=0
        nb_col = self.settings.nb_col
        nb_row = self.settings.nb_row
        
        #for y in range((nb_row-1), -1, -1):
        for y in range(0, nb_row): 
            for x in range(0, nb_col): 
            
                i += 1                
                pos = '%03i' % i
                id = '%s--%s' % (labtekId, pos)
                if resD.has_key(id):
                    value = resD[id][pheno]
                    #moviePageName = os.path.join('.', labtekId, '%s.html' % id)
                    moviePageName = self.getMoviePageName(id)
                    movieNames = resD[id]['movieNames']
                    #moviePageName = os.path.join('.', labtekId, os.path.basename(movieName).replace('.avi', '.html'))
                    gene = resD[id][self.settings.gene_key]
                    sirna = resD[id][self.settings.sirna_key]
                    print "%s: %i/%i, pos: %s -- %s -- %s -- %s -- objects: %i" % (pheno, i, 308, id, pos, sirna, gene, value)
                    self.makeMovieLinkPage(id, movieNames, moviePageName, gene, sirna)
                    tempStr += \
"""<area shape="rect" coords="%i, %i, %i, %i" alt="" target="_blank" 
    title="%i/%i, pos: %s -- %s -- %s -- value: %.2f" href="%s">
""" % ((x_offset + x*x_step), (y_offset + y * y_step), (x_offset + (x + 1) * x_step), (y_offset + (y + 1) * y_step), i, (nb_col*nb_row), pos, gene, sirna, value, moviePageName)

                else:
                    tempStr += \
"""<area shape="rect" coords="%i, %i, %i, %i" alt="" title="%i/384, pos: %s no data">
""" % ((x_offset + x*x_step), (y_offset + y * y_step), (x_offset + (x + 1) * x_step), (y_offset + (y + 1) * y_step), i, pos)
 
        tempStr += \
"""
</map>
"""       

        tempStr += \
"""
<td align="center"><img src="%s" border="0" usemap="#map_%i"></td>
""" % (imageName, map_index)

        return tempStr
#
#    # make single 
#    def OLDmakeImageMapHTML(self, labtekId, resD, imageName, pheno='initCellCount', map_index=0):        
#        
#        tempStr = \
#"""                
#<map name="map_%i">
#""" % map_index
#        
#        y_step = 20 # 776.0 / 32.0
#        x_step = 20 # 276.0 / 12.0
#
#        y_offset = 1
#        x_offset = 22
#
#        i=0
#        nb_col = self.settings.nb_col
#        nb_row = self.settings.nb_row
#        
#        for y in range((nb_row-1), -1, -1):
#            for x in range(0, nb_col): 
#            
#                i += 1                
#                pos = '%03i' % i
#                id = '%s--%s' % (labtekId, pos)
#                if resD.has_key(id):
#                    value = resD[id][pheno]
#                    #moviePageName = os.path.join('.', labtekId, '%s.html' % id)
#                    moviePageName = self.getMoviePageName(id)
#                    movieNames = resD[id]['movieNames']
#                    #moviePageName = os.path.join('.', labtekId, os.path.basename(movieName).replace('.avi', '.html'))
#                    gene = resD[id][self.settings.gene_key]
#                    sirna = resD[id][self.settings.sirna_key]
#                    print "%s: %i/%i, pos: %s -- %s -- %s -- %s -- objects: %i" % (pheno, i, 308, id, pos, sirna, gene, value)
#                    self.makeMovieLinkPage(id, movieNames, moviePageName, gene, sirna)
#                    tempStr += \
#"""<area shape="rect" coords="%i, %i, %i, %i" alt="" target="_blank" 
#    title="%i/%i, pos: %s -- %s -- %s -- value: %.2f" href="%s">
#""" % ((x_offset + x*x_step), (y_offset + y * y_step), (x_offset + (x + 1) * x_step), (y_offset + (y + 1) * y_step), i, (nb_col*nb_row), pos, gene, sirna, value, moviePageName)
#
#                else:
#                    tempStr += \
#"""<area shape="rect" coords="%i, %i, %i, %i" alt="" title="%i/384, pos: %s no data">
#""" % ((x_offset + x*x_step), (y_offset + y * y_step), (x_offset + (x + 1) * x_step), (y_offset + (y + 1) * y_step), i, pos)
# 
#        tempStr += \
#"""
#</map>
#"""       
#
#        tempStr += \
#"""
#<td align="center"><img src="%s" border="0" usemap="#map_%i"></td>
#""" % (imageName, map_index)
#
#        return tempStr
        
        
    def exportPlateHTML(self, filename, labtekId, rowD, rowL):
        file = open(os.path.join(self.html_dir, filename), 'w')
        
        # title and header
        titleString="""
<html>        
  <head><title>Time Lapse Analysis Primary Screen%s</title></head>
  <body>
    <h2>Basic Results for Plate %s </h2>
""" % (labtekId, labtekId)

        file.write(titleString)

        for rowTitle in rowL:
            
            labtekTableString="""
<br clear="all">
<hr width=6000>
<table align="left" border="0" cellspacing="10" cellpadding="0">
"""
            file.write(labtekTableString)

            labtekRowString = ""

            rowHeader = rowD[rowTitle]['header']

            if not rowHeader is None:

                # ----- HEADER OF EACH ROW
                # header: name of labteks
                labtekRowString += """
<tr>
<th align="center"> </th>
"""
                for colTitle in rowHeader:
                    labtekRowString += """
<th align="center"> %s </th>
""" % colTitle
                labtekRowString += """
</tr>
"""
                        
            
            # ----- FIRST COLUMN (ROWTITLE) 
            # write the title of the row 
            if not rowD[rowTitle].has_key('header_link'):
                labtekRowString += """
<tr> <td align="center" valign="middle"> <B> <font size=3> %s </font>   </B> </font> 
""" % rowTitle  
            else:
                # HERE I AM
                labtekRowString += """
<tr> <td align="center" valign="middle"> <B> <font size=3> <A HREF="%s" target="_blank"> %s </A>
</font>   </B> </font> 
""" % (rowD[rowTitle]['header_link'], rowTitle)  
            

                        
            # ----- OTHER COLUMNS
            imageL = rowD[rowTitle]['images']
            if rowD[rowTitle].has_key('height'):
                heightL = rowD[rowTitle]['height']
            else: 
                heightL = [None for x in imageL]
            for img_filename, height in zip(imageL, heightL):
                # no map: plot simply the image
                if not rowD[rowTitle].has_key('map'):
                    # if the height parameter is not given
                    if not height is None:
                        labtekRowString += """
<td align="center"><img src="%s" border="0" height=%i></td>
""" % (img_filename, height)
                    
                    # if the height parameter is given
                    else:
                        labtekRowString += """
<td align="center"><img src="%s" border="0"></td>
""" % img_filename

                # if there is a map, but not for this image
                elif not rowD[rowTitle]['map'].has_key(img_filename):#img_filename != rowD[rowTitle]['map'][0]:

                    # if there is a height parameter
                    if not height is None:
                        labtekRowString += """
<td align="center"><img src="%s" border="0" height=%i></td>
""" % (img_filename, height)

                    # if there is no height paramter
                    else:
                        labtekRowString += """
<td align="center"><img src="%s" border="0"></td>
""" % img_filename

                # if there is a map, and it is for this image
                else:
                    labtekRowString += rowD[rowTitle]['map'][img_filename]
                    
                    
                    
            labtekRowString += """
</tr>
"""
        
            labtekRowString += """
</table>
"""
            # ----- WRITE EVERYTHING TO THE FILE
            file.write(labtekRowString)

        endString = """
</body>
</html>
"""
        file.write(endString)
        file.close()
        
        return


# plots array data
class ArrayPlotter(object):
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
                return max(valuesD[self.phenoClass], 0)
            except KeyError:
                return 0


    def prepareData(self, resD, getLabel, expL=None):
        if expL ==None:
            expL = resD.keys()
            
        nb_exp = len(expL)
        hitValueL = [0 for x in range(nb_exp)]

        for i in range(len(expL)):
            id = expL[i]
            if resD.has_key(id):
                label = getLabel(resD[id])
            else:
                label = 0
            hitValueL[i] = label

        #hitData = rpy.reshape(hitValueL, (self.nb_col, self.nb_row))
        hitData=np.reshape(hitValueL, (self.nb_row, self.nb_col))        

        #hitData = r.matrix(hitValueL, nrow=self.nb_row, byrow=True)[self.nb_row:1,]
        #hitData = r("function(x, nb_row, nb_col) matrix(x, nrow=nb_row, byrow=TRUE)[nb_row:1,]")(hitValueL, self.nb_row, self.nb_col)
        #matrix(seq(1,40), nrow=5, byrow=TRUE)[,8:1]
        
        return hitData    


    # old
#    def plotLabtekArray(self, hitDataParam, filename, plotDir=None, labeledPosD=None,
#                     title='', type='png',
#                     label_colors = True, colVecL=None, legend_plot = False, legend_filename=None, 
#                     control_color='black', numeric_as_int = False,
#                     grid=True, max_val=0.3, min_val=0.0, cut_data = True):
#        if plotDir is None:
#            plotDir = self.plotDir
#        full_filename = os.path.join(plotDir, filename)
#    
#        #hitData = copy.deepcopy(hitDataParam)                
#        hitData = hitDataParam
#        #nrow = len(hitData)
#        #ncol = len(hitData[0])
#        nrow = self.nb_row
#        ncol = self.nb_col
#        
#        xL = range(1, self.nb_col + 1)
#        yL = range(1, self.nb_row + 1)
#        
#        
#        rdev = plot_utilities.RDevice(name = full_filename, title=title, plotType=type, 
#            width=22+20*ncol, height=22+20*nrow)
#
#        if label_colors:
#            # in this case, the data is considered to consist of labels
#            labelL = []
#            for i in range(nrow):
#                for j in range(ncol):
#                    if hitData[i][j] not in labelL:
#                        labelL.append(hitData[i][j]) 
#            
#            if colVecL is None:
#                colVecL = r.rainbow(len(labelL))
#            colBreaksL = range(len(colVecL) + 1)
#
#            if legend_plot:
#                self.plotLabelLegend(colVecL, plotType=type, filename=legend_filename)
#        else:
#            # in this case, the data is numeric.
#            if cut_data:
#                max_rel_cell_count = max_val
#                min_rel_cell_count = 0.0
#                for i in range(nrow):
#                    for j in range(ncol):
#                        hitData[i][j] = max(min_rel_cell_count, min(max_rel_cell_count, hitData[i][j])) 
#            else:
#                max_rel_cell_count = max([max(x) for x in hitData.tolist() ]) 
#                min_rel_cell_count = min([min(x) for x in hitData.tolist() ])
#            
#            if numeric_as_int:
#                nb_colors = max_rel_cell_count
#            else:
#                nb_colors = 500
#
#            if colVecL is None:                    
#                pattern = [(0,0,0),(0.7,0,0),(1,1,0),(1,1,1)]
#                colVecL = colors.make_colors(pattern, nb_colors)
#            colBreaksL = [1.0/ (len(colVecL) - 1) * x * (max_rel_cell_count - min_rel_cell_count) + min_rel_cell_count  for x in range(len(colVecL) + 1)]
#
#            if legend_plot:
#                self.plotNumLegend(colVecL, colBreaksL, 16, filename=legend_filename, type=type, int_labels=numeric_as_int, legendDir = plotDir)
#
#        axisSize = .7
#        r("par(mar=c(1.6,1.6,0.1,0.1))")
#        r.image(xL, yL, r.t(hitData), axes = False, ann=False, cex=1, col=colVecL, breaks=colBreaksL)
#        r.box()
#        if not labeledPosD is None:
#            for label in labeledPosD.keys():
#                posL = labeledPosD[label]
#                if len(posL) > 0:
#                    xlL = [(int(pos)-1) % self.nb_col + 1 for pos in posL]
#                    ylL = [(int(pos)-1) / self.nb_col + 1 for pos in posL]
#                    r.points(xlL, ylL, pch=label, col=control_color, cex=axisSize)
#        # grid
#        if grid:
#            for i in range(self.nb_col):
#                r.abline(v=i+.5, lty=3, lwd=1, col='grey')
#            for i in range(self.nb_row):
#                r.abline(h=i+.5, lty=3, lwd=1, col='grey')
#            
#        r.axis(1, at=xL, labels=[str(x) for x in xL], tick=False, line=-1.0, cex_axis=axisSize)
#        r.axis(2, at=yL, labels=[str(y) for y in yL], tick=False, line=-1.0, cex_axis=axisSize)
#        rdev.close()
#        
#        return

    def plotArray(self, hitDataParam, filename, plotDir=None, labeledPosD=None,
                  title='', type_='png',
                  label_colors = True, colVecL=None, legend_plot = False, legend_filename=None, 
                  control_color='black', numeric_as_int = False,
                  grid=True, max_val=0.3, min_val=0.0, cut_data = True):
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
    
        #hitData = copy.deepcopy(hitDataParam)                
        hitData = hitDataParam
        #nrow = len(hitData)
        #ncol = len(hitData[0])
        nrow_ = self.nb_row
        ncol = self.nb_col
        
        xL =robjects.IntVector( range(1, self.nb_col + 1))
        yL =robjects.IntVector(  range(1, self.nb_row + 1))

        
        rdev = RDevice.png(full_filename, #title=title, plotType=type, 
            width=22+20*ncol, height=22+20*nrow_)

        if label_colors:
            # in this case, the data is considered to consist of labels
            labelL = []
            for i in range(nrow_):
                for j in range(ncol):
                    if hitData[i][j] not in labelL:
                        labelL.append(hitData[i][j]) 
            
            if colVecL is None:
                colVecL = r.rainbow(len(labelL))
            colBreaksL = robjects.IntVector( range(len(colVecL) + 1))

            if legend_plot:
                self.plotLabelLegend(colVecL, filename=legend_filename)
        else:
            # in this case, the data is numeric.
            if cut_data:
                max_rel_cell_count = max_val
                min_rel_cell_count = 0.0
                for i in range(nrow_):
                    for j in range(ncol):
                        hitData[i,j] = max(min_rel_cell_count, min(max_rel_cell_count, hitData[i,j])) 
            else:
                max_rel_cell_count = max([max(x) for x in list(hitData) ]) 
                min_rel_cell_count = min([min(x) for x in list(hitData) ])
            
            if numeric_as_int:
                nb_colors = max_rel_cell_count
            else:
                nb_colors = 500

            if colVecL is None:                    
                pattern = [(0,0,0),(0.7,0,0),(1,1,0),(1,1,1)]
                colVecL = r["rainbow"](nb_colors)#colors.make_colors(pattern, nb_colors)
            colBreaksL = robjects.FloatVector( [1.0/ (len(colVecL) - 1) * x * (max_rel_cell_count - min_rel_cell_count) + min_rel_cell_count  
                                                for x in range(len(colVecL) + 1)])

            if legend_plot:
                self.plotNumLegend(colVecL, colBreaksL, 16, filename=legend_filename, int_labels=numeric_as_int, legendDir = plotDir)

        axisSize = .7
        if type(hitData[0][0])==np.int64:
            h = r['matrix'](robjects.IntVector(hitData.flatten()), nrow=hitData.shape[1])
        else:
            h = r['matrix'](robjects.FloatVector(hitData.flatten()), nrow=hitData.shape[1])

        r("par(mar=c(1.6,1.6,0.1,0.1))")
        r.image(xL, yL, h, axes = False, ann=False, cex=1, col=colVecL, breaks=colBreaksL)
        r.box()
        if not labeledPosD is None:
            for label in labeledPosD.keys():
                posL = labeledPosD[label]
                if len(posL) > 0:
                    xlL = [(int(pos)-1) % self.nb_col + 1 for pos in posL]
                    ylL = [self.nb_row - (int(pos)-1) / self.nb_col for pos in posL]
                    r.points(xlL, ylL, pch=label, col=control_color, cex=axisSize)
        # grid
        if grid:
            for i in range(self.nb_col):
                r.abline(v=i+.5, lty=3, lwd=1, col='grey')
            for i in range(self.nb_row):
                r.abline(h=i+.5, lty=3, lwd=1, col='grey')
            
        r("par(las=1)")
        r.axis(1, at=xL, labels=[str(x) for x in xL], tick=False, cex_axis=axisSize,
               line=-1)
        r.axis(2, at=yL, labels=[r.LETTERS[i] for i in range(len(yL)-1, -1, -1)], 
               tick=False, cex_axis=axisSize, hadj=True, line=-0.5)
        RDevice.dev_off()#rdev.close()
        
        return

    
    def plotNumLegend(self, colVecL, breakL, nb_breaks, filename=None, legendDir=None, int_labels=False):

        if filename is None:
            filename = 'legend_%i_%i_%i' % (len(colVecL), min(breakL), max(breakL))

        if legendDir is None:
            legendDir = self.legendDir
            
        full_filename = os.path.join(legendDir, filename)

        max_break = max(breakL)
        min_break = min(breakL)
        tickBreakL = robjects.FloatVector([float(x) / (nb_breaks - 1) * (max_break - min_break) - min_break for x in range(nb_breaks)])
        if int_labels:
            labels =robjects.IntVector( ['%i' % int(x) for x in tickBreakL])
        else:
            labels = robjects.FloatVector(['%.3f' % x for x in tickBreakL])
        
        rdev = RDevice.png(full_filename,# title='', plotType=type, 
            width=640, height=120)
        r("par(mar=c(3,1,0,1), las=2)")
        #legendA = rpy.reshape(breakL, (len(breakL), 1))
        legendA = r.matrix(breakL, byrow=False, ncol=1)
        r.image(legendA, 1, legendA, col=colVecL, axes=False, ann=False)
        r.box()
        r.axis(1, at=tickBreakL, labels=labels, tick=True, line=0, cex=0.8, cex_axis=0.8)
        
        RDevice.dev_off()
        #rdev.close()
    
        return

    # labelPhenoD = {0: 'no phenotype', 1: 'mitotic', 2: 'shape', 3: 'mitotic and shape'}
    def plotLabelLegend(self, colVecL, labelPhenoD, filename=None, legendDir=None, type='png'):

        if filename is None:
            filename = 'legend_label_%i' % len(labelPhenoD.keys())
        if legendDir is None:
            legendDir = self.legendDir
            
        full_filename = os.path.join(legendDir, filename)

        labelL = labelPhenoD.keys()
        labelL.sort()
        phenoL = [labelPhenoD[x] for x in labelL]
        
        rdev = RDevice.png(full_filename, #title='', plotType=type, 
            width=300, height=200)
        r("par(mar=c(0,0,0,0))")
        r.plot([1, 2, 3], type="n", xlab="", ylab="", main="", ann=False, axes=False)
        r.legend(x="center", legend = phenoL, fill=colVecL, bg = "white", bty="n", cex=0.7)
        RDevice.dev_off()#rdev.close()
    
    
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-s", "--settings_file", dest="settings_file",
                      help="name of the settings file")
    (options, args) = parser.parse_args() 
    htmlgen = HTMLGenerator(options.settings_file)
    htmlgen()
    
    #parser.add_option("-g", "--group", dest="group",
    #                  help="group")

#    parser.add_option("-l", "--labtek", dest="labtek",
#                      help="labtek")
#                      
#
#    (options, args) = parser.parse_args()        
#    group = options.group
#    labtek = options.labtek
#    
#    LOCAL_PLOT_DIR = os.path.join('/netshare/almf/EMBO-data/timelapse_result', group, 'html')
#    HTML_DIR = LOCAL_PLOT_DIR
#    RESULT_DIR = os.path.join('/netshare/almf/EMBO-data/timelapse_result', group)
#    
#    ht = HTMLGenerator(html_dir=HTML_DIR, baseDir=RESULT_DIR)
#    ht(labtek)
#    
    
