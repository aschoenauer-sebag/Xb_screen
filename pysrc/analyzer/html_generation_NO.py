import numpy, os
import copy
import pickle
import random

from settings import *
from plot_utilities import *
from util import *
import quality_control
import spotInfo

from spotInfo import spotDefinition
reload(spotInfo)

#class RDevice:
#    
#    def __init__(self, name=None, type="pdf", pointsize="10", title="",
#                 **optionD):
#        self.type = type
#        self.name = name
#        if not name is None:
#            name = ".".join((name, type))
#            if type == "pdf":
#                r.pdf(name, title=title, pointsize=pointsize, **optionD)
#            elif type == "ps":
#                r.postscript(name, title=title, **optionD)
#            elif type == "png":
#                r.png(name, **optionD)
#
#    def close(self):
#        if not self.name is None:
#            # ATTENTION! call 'dev.off()' to write all plot-data to file
#            r.dev_off()


    
        
# plots array data
class ArrayPlotter(object):
    def __init__(self, plotDir=os.path.join(LOCAL_PLOT_DIR, 'array'), legendDir = os.path.join(LOCAL_PLOT_DIR, 'legend')):
        self.plotDir = plotDir
        self.legendDir = legendDir
        return
    
    class getPhenoVal(object):

        def __init__(self, phenoClass):
            self.phenoClass = phenoClass

        def __call__(self, valuesD):
            return max(valuesD[self.phenoClass], 0)

    def prepareData(self, resD, getLabel, expL):

        if not len(expL) == 384:
            return None
        
        hitValueL = [0 for x in range(384)]

        for i in range(len(expL)):
            id = expL[i]
            if resD.has_key(id):
                label = getLabel(resD[id])
            else:
                label = 0
            hitValueL[i] = label

        hitData = numpy.reshape(hitValueL, (32, 12))
        
        return hitData    


    def plotArray(self, hitDataParam, filename, plotDir=None, labeledPosD=None,
                  title='', type='png',
                  label_colors = True, colVecL=None, legend_plot = False, legend_filename=None, control_color='black', numeric_as_int = False,
                  grid=True, max_val=0.3, min_val=0.0, cut_data = True):
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
    
        hitData = copy.deepcopy(hitDataParam)                
        nrow = len(hitData)
        ncol = len(hitData[0])
        
        xL = range(1, 33)
        yL = range(1, 13)
        
        rdev = RDevice(name = full_filename, title=title, type=type, 
                       width=640, height=250)

        if label_colors:
            # in this case, the data is considered to consist of labels
            labelL = []
            for i in range(nrow):
                for j in range(ncol):
                    if hitData[i][j] not in labelL:
                        labelL.append(hitData[i][j]) 
            
            if colVecL is None:
                colVecL = r.rainbow(len(labelL))
            colBreaksL = range(len(colVecL) + 1)

            if legend_plot:
                self.plotLabelLegend(colVecL, type=type, filename=legend_filename)
        else:
            # in this case, the data is numeric.
            if cut_data:
                max_rel_cell_count = max_val
                min_rel_cell_count = 0.0
                for i in range(nrow):
                    for j in range(ncol):
                        hitData[i][j] = max(min_rel_cell_count, min(max_rel_cell_count, hitData[i][j])) 
            else:
                max_rel_cell_count = max([max(x) for x in hitData.tolist() ]) 
                min_rel_cell_count = min([min(x) for x in hitData.tolist() ])
            
            if numeric_as_int:
                nb_colors = max_rel_cell_count
            else:
                nb_colors = 500
                
            pattern = [(0,0,0),(0.7,0,0),(1,1,0),(1,1,1)]
            colVecL = colors.make_colors(pattern, nb_colors)
            colBreaksL = [1.0/ (len(colVecL) - 1) * x * (max_rel_cell_count - min_rel_cell_count) + min_rel_cell_count  for x in range(len(colVecL) + 1)]

            if legend_plot:
                self.plotNumLegend(colVecL, colBreaksL, 16, filename=legend_filename, type=type, int_labels=numeric_as_int)

        axisSize = .5
        r("par(mar=c(1.6,1.6,0.1,0.1))")
        r.image(xL, yL, hitData, axes = False, ann=False, cex=1, col=colVecL, breaks=colBreaksL)        
        r.box()
        if not labeledPosD is None:
            for label in labeledPosD.keys():
                posL = labeledPosD[label]
                if len(posL) > 0:
                    xlL = [(int(x)-1) / 12 + 1 for x in posL]
                    ylL = [(int(x)-1) % 12 + 1 for x in posL]
                    r.points(xlL, ylL, pch=label, col=control_color, cex=axisSize)                

        # grid
        if grid:
            for i in range(12):
                r.abline(h=i+.5, lty=3, lwd=1, col='grey')
            for i in range(32):
                r.abline(v=i+.5, lty=3, lwd=1, col='grey')
            
        r.axis(1, at=xL, labels=[str(x) for x in xL], tick=False, line=-1.0, cex_axis=axisSize)
        r.axis(2, at=yL, labels=[str(y) for y in yL], tick=False, line=-1.0, cex_axis=axisSize)
        rdev.close()
        
        return

    
    def plotNumLegend(self, colVecL, breakL, nb_breaks, filename=None, legendDir=None, type='png', int_labels=False):

        if filename is None:
            filename = 'legend_%i_%i_%i' % (len(colVecL), min(breakL), max(breakL))

        if legendDir is None:
            legendDir = self.legendDir
            
        full_filename = os.path.join(legendDir, filename)

        max_break = max(breakL)
        min_break = min(breakL)
        tickBreakL = [float(x) / (nb_breaks - 1) * (max_break - min_break) - min_break for x in range(nb_breaks)]
        if int_labels:
            labels = ['%i' % int(x) for x in tickBreakL]
        else:
            labels = ['%.3f' % x for x in tickBreakL]
        
        rdev = RDevice(name = full_filename, title='', type=type, 
                       width=640, height=120)
        r("par(mar=c(3,1,0,1), las=2)")
        legendA = numpy.reshape(breakL, (len(breakL), 1))
        r.image(legendA, 1, legendA, col=colVecL, axes=False, ann=False)
        r.box()
        r.axis(1, at=tickBreakL, labels=labels, tick=True, line=0, cex=0.8, cex_axis=0.8)
        rdev.close()
    
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
        
        rdev = RDevice(name = full_filename, title='', type=type, 
                       width=300, height=200)
        r("par(mar=c(0,0,0,0))")
        r.plot([1, 2, 3], type="n", xlab="", ylab="", main="", ann=False, axes=False)
        r.legend(x="center", legend = phenoL, fill=colVecL, bg = "white", bty="n", cex=0.7)
        rdev.close()
    

class plotHitsOnLabteks(object):

    def __init__(self, metaDir=META_RESULT_DIR, qcDir=QC_DIR, spotqcD=None):
        self.metaDir = metaDir
        self.subsetFilename = os.path.join(self.metaDir, 'hit_id_subset.pickle')
        self.subsetShortFilename = os.path.join(self.metaDir, 'hit_id_subset_short.pickle')
        
        self.plotDir = '/netshare/sputnik/plots/hits_on_plate'
        #self.plotDir = '/data/twalter/Plots/hits_on_plate'
        #self.htmlDir = '/data/twalter/WWW'
        self.htmlDir = '/netshare/sputnik/www/localization'
        self.relativePlotDir = '../Plots/hits_on_plate'
        
        if not os.path.isdir(self.plotDir):
            os.mkdir(self.plotDir)
        if not os.path.isdir(self.htmlDir):
            os.makedirs(self.htmlDir)

        self.thresholdSet2D = {'maxSum_Shape':0.132, 'maxSum_MitosisPhenotype': 0.04}

        # import QC data
        self.qc = quality_control.QC()
        self.spotQC = quality_control.spotQC(self.qc)

        # spotqcD stores the information if a spot is accepted for evaluation or not
        if spotqcD is None:
            self.spotqcD = self.spotQC.readSpotQCFile()
        else:
            self.spotqcD = spotqcD


        # the list of valid labteks
        self.goodL = self.qc.getPassedLabtekList()    
        self.goodL.sort()        

        self.labtekBaseD = {}
        for labtekId in self.goodL:
            labtekBase = labtekId[0:len('LT0001')]
            if not self.labtekBaseD.has_key(labtekBase):
                self.labtekBaseD[labtekBase] = []
            self.labtekBaseD[labtekBase].append(labtekId)
                
        self.man = None

        return

    def getManualAnnotation(self):
        if self.man is None:
            self.man = manual_annotation.manualAnnotationObject()
        
        # read manual annotation
        manualAnnotationD = self.man.readAnnotationScore()
        
        # adaption for manual annotation
        for id in manualAnnotationD.keys():
            if not self.spotqcD[id]['accepted']:
                del(manualAnnotationD[id])
                
        idScoreD = {}
        for id in manualAnnotationD.keys():
            idScoreD[id] = manualAnnotationD[id]['score']

        return idScoreD
    
    def getShortSubsetFromHitList(self, resD, hitListD):
        # find all experiments from the hitListD structure
        subsetD = {}
        for pheno in hitListD.keys():
            for sirnaHitD, score in hitListD[pheno]:
                sirnaExpL = [x[0] for x in sirnaHitD['expL']]
                for exp in sirnaExpL:
                    lt, pos = exp.split('--')
                    ltBase = lt[0:len('LT0001')]
                    new_key = '%s--%s' % (ltBase, pos)
                    if not subsetD.has_key(new_key):
                        subsetD[new_key] = []
                    if pheno not in subsetD[new_key]:
                        subsetD[new_key].append(pheno)                    
        return subsetD

    def getSubsetFromHitList(self, hitListD):
        # find all experiments from the hitListD structure
        expL = []
        for pheno in hitListD.keys():
            for sirnaHitD, score in hitListD[pheno]:
                sirnaExpL = [x[0] for x in sirnaHitD['expL']]
                expL.extend(sirnaExpL)

        subsetD = {}
    
        for exp in expL:
            subsetD[exp] = resD[exp]

        return subsetD
    
    def showHitInformation(self, hitListD, sirnaId):
        for phenoClass in hitListD.keys():
            for sirnaHitD, score in hitListD[phenoClass]:
                if sirnaHitD['sirnaId'] == sirnaId:
                    print str(sirnaHitD)            
        return
    
    def dumpSubset(self, subsetD):
        file = open(self.subsetFilename, 'w')
        pickle.dump(subsetD, file)
        file.close()
        return

    def dumpSubsetShort(self, subsetD):
        file = open(self.subsetShortFilename, 'w')
        pickle.dump(subsetD, file)
        file.close()
        return

    def readSubset(self):
        file = open(self.subsetFilename, 'r')
        subsetD = pickle.load(file)
        file.close()
        return subsetD

    def readSubsetShort(self):
        file = open(self.subsetShortFilename, 'r')
        subsetD = pickle.load(file)
        file.close()
        return subsetD

    class getPhenoVal(object):

        def __init__(self, phenoClass):
            self.phenoClass = phenoClass

        def __call__(self, valuesD):
            return max(valuesD[self.phenoClass], 0)       
        
    def getLabel2D(self, valuesD):
        phenoClassL = self.thresholdSet2D.keys()
        phenoClassL.sort()

        labelD = dict(zip(phenoClassL, [2**x for x in range(len(phenoClassL))])) 
        label = 0
        for phenoClass in phenoClassL:
            if valuesD[phenoClass] >= self.thresholdSet2D[phenoClass]:
                label += labelD[phenoClass]
        #print '%s : %.4f\t%s : %.4f\t --> \tlabel : %i' % (phenoClassL[0], max(valuesD[phenoClassL[0]], 0), phenoClassL[1], max(valuesD[phenoClassL[1]], 0), label) 
        #print labelD
        return label
    
    def getLabelFromShortDescription(self, hitL):
        label = 0
        
        for hit in hitL:
            if hit == 'maxSum_MitosisPhenotype':
                label += 1
            if hit == 'maxSum_Shape':
                label += 2

        return label
        
    def getColorVec2D_old(self, valueL):        
        tempL = ['white', 'green', 'blue', 'cyan']
        colorVecL = [tempL[i] for i in valueL]
        return colorVecL            

    def getColorVec2D(self):        
        colorVecL = ['white', 'green', 'blue', 'cyan']
        return colorVecL            
                    
    def prepareHitData(self, resD, getLabel, labtekId):

        #expL = filter(lambda(x): x[0:len(labtekId)] == labtekId, resD.keys())
        posL = ['%03i' % x for x in range(1, 385)]
        expL = ['%s--%s' % (labtekId, pos) for pos in posL]
        expL.sort()
        
        hitValueL = [0 for x in range(384)]

        for id in expL:
            lt, pos = id.split('--')
            if resD.has_key(id):
                label = getLabel(resD[id])
            else:
                label = 0
            hitValueL[int(pos)-1] = label

        hitData = numpy.reshape(hitValueL, (32, 12))
        
        return hitData    

    def hitHisto(self, subsetD):
        expL = subsetD.keys()
        expL.sort()

        phenoL = []
        for exp in expL:
            for pheno in subsetD[exp]:
                if not pheno in phenoL:
                    phenoL.append(pheno)

        print phenoL
        hitValueD = {}
        for config in ['druggable', 'genome']:
            hitValueD[config] = {}
            for pheno in phenoL:
                hitValueD[config][pheno] = [0 for x in range(384)]
        
        for exp in expL:
            lt, pos = exp.split('--')
            config = self.getConfig(lt)
            
            for pheno in subsetD[exp]:
                hitValueD[config][pheno][int(pos) - 1] += 1

        hitDataD = {}
        for config in ['druggable', 'genome']:
            hitDataD[config] = {}
            for pheno in phenoL:
                hitDataD[config][pheno] = numpy.reshape(hitValueD[config][pheno], (32, 12))

        return hitDataD

    def getConfig(self, labtekId):
        if int(labtekId[2:6]) < 50:
            return 'druggable'
        else:
            return 'genome'

    def plotGenomeStatOnLabtek(self, resD):

        statD = {}
        # mean
        for pheno in ['maxSum_MitosisPhenotype', 'maxSum_Shape']:
            print ' ---- %s ----' % pheno
            statD[pheno] = {}

            # mean for druggable
            print '%s:\t%s\t%s' % (pheno, 'mean', 'druggable')
            druggableL = filter(lambda(x): int(x[2:6]) < 50, self.goodL)
            hitData = self.prepareDataAllInOne(resD, self.getPhenoVal(pheno), stats.mean, labtekL=druggableL)
            self.plotLocalDistribution('druggable', hitData, file_basename='mean_%s_' % pheno,
                                       title=pheno, label_colors=False, grid=False)
            statD[pheno]['mean_druggable'] = hitData
            
            # mean for genome wide
            print '%s:\t%s\t%s' % (pheno, 'mean', 'genome')
            genomeL = filter(lambda(x): int(x[2:6]) > 50, self.goodL)
            hitData = self.prepareDataAllInOne(resD, self.getPhenoVal(pheno), stats.mean, labtekL=genomeL)
            self.plotLocalDistribution('genome', hitData, file_basename='mean_%s_' % pheno,
                                       title=pheno, label_colors=False, grid=False)
            statD[pheno]['mean_genome'] = hitData

            # median for druggable
            print '%s:\t%s\t%s' % (pheno, 'median', 'druggable')
            hitData = self.prepareDataAllInOne(resD, self.getPhenoVal(pheno), stats.median, labtekL=druggableL)
            self.plotLocalDistribution('druggable', hitData, file_basename='median_%s_' % pheno,
                                       title=pheno, label_colors=False, grid=False)
            statD[pheno]['median_druggable'] = hitData
            
            # median for genome wide
            print '%s:\t%s\t%s' % (pheno, 'median', 'genome')
            hitData = self.prepareDataAllInOne(resD, self.getPhenoVal(pheno), stats.median, labtekL=genomeL)
            self.plotLocalDistribution('genome', hitData, file_basename='median_%s_' % pheno,
                                       title=pheno, label_colors=False, grid=False)
            statD[pheno]['median_genome'] = hitData

            # stdev for druggable
            print '%s:\t%s\t%s' % (pheno, 'stdev', 'druggable')
            druggableL = filter(lambda(x): int(x[2:6]) < 50, self.goodL)
            hitData = self.prepareDataAllInOne(resD, self.getPhenoVal(pheno), stats.stdev, labtekL=druggableL)
            self.plotLocalDistribution('druggable', hitData, file_basename='stdev_%s_' % pheno,
                                       title=pheno, label_colors=False, grid=False)
            statD[pheno]['stdev_druggable'] = hitData
            
            # stdev for genome wide
            print '%s:\t%s\t%s' % (pheno, 'stdev', 'genome')
            genomeL = filter(lambda(x): int(x[2:6]) > 50, self.goodL)
            hitData = self.prepareDataAllInOne(resD, self.getPhenoVal(pheno), stats.stdev, labtekL=genomeL)
            self.plotLocalDistribution('genome', hitData, file_basename='stdev_%s_' % pheno,
                                       title=pheno, label_colors=False, grid=False)
            statD[pheno]['stdev_genome'] = hitData
            
        return statD

    def plotGenomeStatOnLabtekInNumbers(self, statD):
        for pheno in statD.keys():
            print '%s:\t%s\t%s' % (pheno, 'mean', 'druggable')
            druggableL = filter(lambda(x): int(x[2:6]) < 50, self.goodL)
            hitData = statD[pheno]['mean_druggable']
            self.plotNumbersOnLabtek('druggable', hitData, file_basename='numbers_mean_%s_' % pheno,
                                      title=pheno)

            print '%s:\t%s\t%s' % (pheno, 'mean', 'genome')
            genomeL = filter(lambda(x): int(x[2:6]) < 50, self.goodL)
            hitData = statD[pheno]['mean_genome']
            self.plotNumbersOnLabtek('genome', hitData, file_basename='numbers_mean_%s_' % pheno,
                                      title=pheno)

            print '%s:\t%s\t%s' % (pheno, 'median', 'druggable')
            druggableL = filter(lambda(x): int(x[2:6]) < 50, self.goodL)
            hitData = statD[pheno]['median_druggable']
            self.plotNumbersOnLabtek('druggable', hitData, file_basename='numbers_median_%s_' % pheno,
                                      title=pheno)

            print '%s:\t%s\t%s' % (pheno, 'median', 'genome')
            genomeL = filter(lambda(x): int(x[2:6]) < 50, self.goodL)
            hitData = statD[pheno]['median_genome']
            self.plotNumbersOnLabtek('genome', hitData, file_basename='numbers_median_%s_' % pheno,
                                      title=pheno)

        return

    def printNeighborValues(self, statD, config='druggable', pheno='maxSum_MitosisPhenotype', spotType='eg5Neighbors'):
        # druggable
        reload(util)
        if config[0:len('drug')] == 'drug':
            spotDef = spotDefinition('LT0001')
            statEntry = 'mean_druggable'
        else:
            spotDef = spotDefinition('LT0060')
            statEntry = 'mean_genome'

        print ' --------------- '
        print 'configuration: %s' % config
        eg5Nb = spotDef.typeD[spotType]
        coordL = []
        for pos in eg5Nb:
            x, y = int(pos)/12 + 1, int(pos)%12
            coordL.append((x,y))
        coordL.sort(util.lexicoAsc)
        for x,y in coordL:
            val = statD[pheno][statEntry][x-1][y-1]
            print '(%i,%i):\t%.3f' % (x, y, val)

        return
    
    # calculates for the whole screen the calcVal for each of the positions.
    # getLabel defines which measure is calculated for each spot.
    def prepareDataAllInOne(self, resD, getLabel, calcVal, labtekL=None):
        if labtekL is None:
            labtekL = self.goodL

        posL = ['%03i' % x for x in range(1, 385)]
        expL = resD.keys()
        expL.sort()

        hitValueL = [0 for x in range(384)]

        for pos in posL:
            print 'processing %s' % pos
            subExpL = filter(lambda(x): x.split('--')[1] == pos and x.split('--')[0] in labtekL, expL)            
            valL = []
            for id in subExpL:
                valL.append(getLabel(resD[id]))

            hitValueL[int(pos)-1] = calcVal(valL)
            
        hitData = numpy.reshape(hitValueL, (32, 12))

        return hitData
    
    def prepareMedianData(self, resD, getLabel, labtekBase):

        labtekL = filter(lambda(x): x[0:len('LT0001')]==labtekBase, self.goodL)
        
        subResD = {}
        hitValueL = [0 for x in range(384)]

        posL = ['%03i' % x for x in range(1, 385)]
        possibleIdL = []
        for labtekId in labtekL:
            possibleIdL.extend(['%s--%s' % (labtekId, pos) for pos in posL])                             
        
        # get a valid id from resD
        for possibleId in possibleIdL:
            if resD.has_key(possibleId):                
                validId = possibleId
                break
        else:
            # in this case there is no valid id from resD
            hitData = numpy.reshape([0 for i in range(384)], (32, 12))
            return hitData
        
        # get the numeric fields from resD
        numericKeyL = filter(lambda(x): not x in ['geneName', 'sirnaId'], resD[validId].keys())

        for pos in posL:
            #for labtekId in labtekL:
            validIdL = []
            for labtekId in labtekL:               
                id = '%s--%s' % (labtekId, pos)
                if resD.has_key(id):
                    validIdL.append(id)
            if len(validIdL) == 0:
                #print 'no valid experiment for %s %s' % (labtekBase, pos)
                label=0
                hitValueL[int(pos)-1] = label
                continue

            firstId = validIdL[0]      
            subResD[pos] = {'geneName': resD[firstId]['geneName'],
                            'sirnaId': resD[firstId]['sirnaId'] }
            for numericKey in numericKeyL:
                vecL = [resD[id][numericKey] for id in validIdL]
                if len(vecL) < 2:
                    subResD[pos][numericKey] = 0.0
                if len(vecL) == 2:
                    subResD[pos][numericKey] = min(vecL)
                if len(vecL) >= 3:
                    subResD[pos][numericKey] = upperMedian(vecL)

            label = getLabel(subResD[pos])
            hitValueL[int(pos)-1] = label

        hitData = numpy.reshape(hitValueL, (32, 12))

        return hitData

    def __call__(self, resD=None, manualAnnotationD=None, scoreD=None):
        # read data
        if resD is None:
            gen = genomeAnalyzer()
            resD = gen.readIdFile(filename_addon='_maxDiff')

        # -----
        # single experiments: plot and generate the html file
        self.plotLocalDistributionAll(resD, manualAnnotationD=manualAnnotationD)
        self.exportHTMLSingleExperiments()

        # -----
        # median and hit detection: plot and generate html file
        self.plotMedianDistributionAll(resD, scoreD=scoreD)
        # subsetD generated from: subsetD = self.getShortSubsetFromHitList(resD, hitListD)
        subsetD = self.readSubsetShort()
        self.plotFinalHitDistributionAll(subsetD, scoreD=scoreD)
        self.exportHTMLMedian()

        # -----
        return    
    
    def getAdditionalMedianData(self, scoreD, labtekBase):
        additionalInfoD = {}
        
        hitSirnaL = filter(lambda(x): scoreD[x]['pheno'] == 1, scoreD.keys())
        posL = []
        for sirna in hitSirnaL:
            idL = scoreD[sirna]['idL']
            for id in idL:
                labtekId, pos = id.split('--')
                currentLabtekBase = labtekId[:-3]                
                if currentLabtekBase==labtekBase and (not pos in posL):
                    posL.append(pos)
        posL.sort()
        
        # manual hits
        infoSource = 'manual_hits'
        additionalInfoD[infoSource] = {}
        additionalInfoD[infoSource]['pch'] = 0
        additionalInfoD[infoSource]['cex'] = 1.2
        additionalInfoD[infoSource]['col'] = 'red'
        additionalInfoD[infoSource]['posL'] = posL
        additionalInfoD[infoSource]['lwd'] = 2

        # manual negative
        negSirnaL = filter(lambda(x): scoreD[x]['pheno'] == 0, scoreD.keys())
        posL = []
        for sirna in negSirnaL:
            idL = scoreD[sirna]['idL']
            for id in idL:
                labtekId, pos = id.split('--')
                currentLabtekBase = labtekId[:-3]                
                if (currentLabtekBase==labtekBase) and (not pos in posL):
                    posL.append(pos)
        posL.sort()
        
        # manual hits
        infoSource = 'manual_negative'
        additionalInfoD[infoSource] = {}
        additionalInfoD[infoSource]['pch'] = 0
        additionalInfoD[infoSource]['cex'] = 1.2
        additionalInfoD[infoSource]['col'] = 'gray'
        additionalInfoD[infoSource]['posL'] = posL
        additionalInfoD[infoSource]['lwd'] = 2

        return additionalInfoD
        
    # transforms manual annotation to additional data (ready for plotLocalDistribution)
    def getAdditionalData(self, manualAnnotationD, labtekId):
        additionalInfoD = {}
        initPosL = ['%s--%03i' % (labtekId, pos) for pos in range(1, 385)]

        # manual hits
        infoSource = 'manual_hits'
        additionalInfoD[infoSource] = {}
        posL = [y.split('--')[1] for y in filter(lambda(x): manualAnnotationD.has_key(x) and manualAnnotationD[x]['score'] > 0, initPosL)]
        additionalInfoD[infoSource]['pch'] = 0
        additionalInfoD[infoSource]['cex'] = 1.2
        additionalInfoD[infoSource]['col'] = 'red'
        additionalInfoD[infoSource]['posL'] = posL
        additionalInfoD[infoSource]['lwd'] = 2

        # manual hits
        infoSource = 'manual_negative'
        additionalInfoD[infoSource] = {}
        posL = [y.split('--')[1] for y in filter(lambda(x): manualAnnotationD.has_key(x) and manualAnnotationD[x]['score'] == 0, initPosL)]
        additionalInfoD[infoSource]['pch'] = 0
        additionalInfoD[infoSource]['cex'] = 1.2
        additionalInfoD[infoSource]['col'] = 'gray'
        additionalInfoD[infoSource]['posL'] = posL
        additionalInfoD[infoSource]['lwd'] = 2

        # manual hits
        infoSource = 'disabled'
        additionalInfoD[infoSource] = {}
        posL = [y.split('--')[1] for y in filter(lambda(x): not self.spotqcD[x]['accepted'], initPosL)]
        additionalInfoD[infoSource]['pch'] = 3
        additionalInfoD[infoSource]['cex'] = 1.0
        additionalInfoD[infoSource]['col'] = 'red'
        additionalInfoD[infoSource]['posL'] = posL
        additionalInfoD[infoSource]['lwd'] = 1
        
        return additionalInfoD
    
    def plotLocalDistributionAll(self, resD, manualAnnotationD=None):
        for labtekId in self.goodL:
            if not manualAnnotationD is None:
                additionalInfoD = self.getAdditionalData(manualAnnotationD, labtekId)
            else: 
                additionalInfoD = {}
            
            hitData = self.prepareHitData(resD, self.getLabel2D, labtekId)
            
            self.plotLocalDistribution(labtekId, hitData, additionalInfoD=additionalInfoD)
            print '%s has been plotted' % labtekId
        return

    def plotMedianDistributionAll(self, resD, scoreD=None):

        labtekBaseL = self.labtekBaseD.keys()
        labtekBaseL.sort()
        
        for labtekBase in labtekBaseL:

            if not scoreD is None:
                additionalInfoD = self.getAdditionalMedianData(scoreD, labtekBase)

            hitData =self.prepareMedianData(resD, self.getLabel2D, labtekBase)
            title = 'Hit localization for %s' % labtekBase
            self.plotLocalDistribution(labtekBase, hitData, file_basename='medianDistribution_',
                                       title=title, additionalInfoD=additionalInfoD)
            print '%s has been plotted' % labtekBase
        return

    def plotFinalHitDistributionAll(self, resD, scoreD=None):

        labtekBaseL = self.labtekBaseD.keys()
        labtekBaseL.sort()
        
        for labtekBase in labtekBaseL:
            if not scoreD is None:
                additionalInfoD = self.getAdditionalMedianData(scoreD, labtekBase)

            hitData =self.prepareHitData(resD, self.getLabelFromShortDescription, labtekBase)
            title = 'Hit localization for %s' % labtekBase
            self.plotLocalDistribution(labtekBase, hitData, file_basename='finalDistribution_',
                                       title=title, additionalInfoD=additionalInfoD)
            print '%s has been plotted' % labtekBase
        return

    def plotPhenoMedianDistributionAll(self, resD, pheno):

        labtekBaseL = self.labtekBaseD.keys()
        labtekBaseL.sort()
        
        for labtekBase in labtekBaseL:
            hitData =self.prepareMedianData(resD, self.getPhenoVal(pheno), labtekBase)
            title = '%s distribution for %s' % (pheno, labtekBase)
            self.plotLocalDistribution(labtekBase, hitData, file_basename='%s_distribution_' % pheno,
                                       title=title, label_colors=False, grid=False)
            print '%s has been plotted' % labtekBase
        return

    def plotPhenoDistributionAll(self, resD, pheno):
        for labtekId in self.goodL:
            hitData =self.prepareHitData(resD, self.getPhenoVal(pheno), labtekId)
            title = '%s distribution for %s' % (pheno, labtekId)
            self.plotLocalDistribution(labtekId, hitData, file_basename='%s_distribution_' % pheno,
                                       title=title, label_colors=False, grid=False)            
            print '%s has been plotted' % labtekId
        return

    def getValuesFromHitData(self, hitData):
        tempL = hitData.tolist()
        resL = []
        for vec in tempL:
            for data in vec:
                if data not in resL:
                    resL.append(data)
        resL.sort()
        return resL

    def getControlData(self, labtekId, controlD=None):
        if controlD is None:
            controlD = {'scrambled': 1,
                        'eg5': 2,
                        'incenp': 3,
                        'bcop': 4,
                        'marker': 5,
                        'empty': 6}
        spotDef = spotDefinition(labtekId)
        
        posL = ['%03i' % x for x in range(1, 385)]
        controlL = [0 for i in range(1, 385)]
        for pos in posL:
            for control in controlD.keys():
                if pos in spotDef.typeD[control]:
                    controlL[int(pos)-1] = controlD[control]

        controlData = numpy.reshape(controlL, (32, 12))

        return controlData            
        
    def plotNumbersOnLabtek(self, labtekId, hitDataParam, file_basename,
                            title=None, type='png', grid=True):
    
        if title is None:
            title = 'local distribution for %s' % labtekId
    
        hitData = copy.deepcopy(hitDataParam)
        
        full_filename = os.path.join(self.plotDir, '%s%s' % (file_basename, labtekId))
    
        xL = range(1, 33)
        yL = range(1, 13)
    
        if labtekId[0:len('dru')] == 'dru':
            labtekId = 'LT0001'
        if labtekId[0:len('gen')] == 'gen':
            labtekId = 'LT0100'
        spotDef = spotDefinition(labtekId)        

        controlData = self.getControlData(labtekId)
        control_max = max(max(controlData.tolist()))
        controlColVecL = ['white', 'grey', 'green', 'lightblue', 'red', 'lightblue', 'yellow']
        controlColVecL = controlColVecL[0:(control_max + 1)]
        print '%i: %s' % (control_max, str(controlColVecL))
        rdev = RDevice(name = full_filename, title=title, type=type, 
                       width=1280, height=500)                    
            
        axisSize = .5
        r("par(mar=c(1.6,1.6,0.1,0.1))")
        r.image(xL, yL, controlData, axes = False, ann=False, cex=1, col=controlColVecL, breaks=[x - 0.5 for x in range(control_max+2)])        
        r.box()
        for x in xL:
            labelL = ['%.3f' % hitData[x-1][y-1] for y in yL]
            r.text([x], yL, label=labelL, col='black', cex=axisSize)
    
        # grid
        if grid:
            for i in range(12):
                r.abline(h=i+.5, lty=3, lwd=1, col='grey')
            for i in range(32):
                r.abline(v=i+.5, lty=3, lwd=1, col='grey')
            
        r.axis(1, at=xL, labels=[str(x) for x in xL], tick=False, line=-1.0, cex_axis=axisSize)
        r.axis(2, at=yL, labels=[str(y) for y in yL], tick=False, line=-1.0, cex_axis=axisSize)
        rdev.close()
        
        return

    def plotHitHisto(self, hitHisto):
        for config in ['druggable', 'genome']:
            for pheno in hitHisto[config].keys():
                self.plotLocalDistribution(config, hitHisto[config][pheno], file_basename='hitHisto_%s_' % pheno,
                                           title='Histogram of hits in %s, configuration: %s' % (pheno, config),
                                           grid=False, label_colors=False, max_val = 15, legend_plot = True,
                                           control_color='white')
        return

    def plotLocalDistribution(self, labtekId, hitDataParam, file_basename='localDistribution_',
                              title=None, type='png', label_colors = True, grid=True, max_val=0.3,
                              legend_plot = False, control_color='black', additionalInfoD = {}):

        if title is None:
            title = 'local distribution for %s' % labtekId

        hitData = copy.deepcopy(hitDataParam)
        
        full_filename = os.path.join(self.plotDir, '%s%s' % (file_basename, labtekId))

        xL = range(1, 33)
        yL = range(1, 13)

        if labtekId[0:len('dru')] == 'dru':
            labtekId = 'LT0001'
        if labtekId[0:len('gen')] == 'gen':
            labtekId = 'LT0100'
        spotDef = spotDefinition(labtekId)        
        
        rdev = RDevice(name = full_filename, title=title, type=type, 
                       width=640, height=250)
        #data = numpy.reshape(range(1, 385), (32, 12))
        if label_colors:
            #valueL = self.getValuesFromHitData(hitData)
            colVecL = self.getColorVec2D()
            colBreaksL = [x - 1 for x in range(len(colVecL) + 1)]

            if legend_plot:
                self.plotLabelLegend(colVecL, type=type)
        else:
            max_rel_cell_count = max_val
            min_rel_cell_count = 0.0
            nrow = len(hitData)
            ncol = len(hitData[0])
            for i in range(nrow):
                for j in range(ncol):
                    hitData[i][j] = max(min_rel_cell_count, min(max_rel_cell_count, hitData[i][j])) 
            if max_val > 5:
                nb_colors = max_val
                int_labels=True
            else:
                nb_colors = 500
                int_labels = False
                
            pattern = [(0,0,0),(0.7,0,0),(1,1,0),(1,1,1)]
            colVecL = colors.make_colors(pattern, nb_colors)
            maxVal = len(colVecL) - 1
            colBreaksL = [1.0/ maxVal * max_rel_cell_count * x for x in range(len(colVecL) + 1)]

            if legend_plot:
                self.plotNumLegend(colVecL, colBreaksL, type=type, int_labels=int_labels)

        axisSize = .5
        r("par(mar=c(1.6,1.6,0.1,0.1))")
        r.image(xL, yL, hitData, axes = False, ann=False, cex=1, col=colVecL, breaks=colBreaksL)
        r.box()
        
        # plot the control letter
        for label in spotDef.labelD.keys():
            posL = spotDef.labelD[label]
            if len(posL) > 0:
                xlL = [(int(x)-1) / 12 + 1 for x in posL]
                ylL = [(int(x)-1) % 12 + 1 for x in posL]
                r.points(xlL, ylL, pch=label, col=control_color, cex=axisSize)

        # additional info like manually found hits or disabled spots
        for infoSource in additionalInfoD.keys():
            #additionalInfoD[infoSource]
            posL = additionalInfoD[infoSource]['posL']
            pch = additionalInfoD[infoSource]['pch']
            cex = additionalInfoD[infoSource]['cex']
            col = additionalInfoD[infoSource]['col']
            lwd = additionalInfoD[infoSource]['lwd']
            if len(posL) > 0:
                xlL = [(int(x)-1) / 12 + 1 for x in posL]
                ylL = [(int(x)-1) % 12 + 1 for x in posL]
                r.points(xlL, ylL, pch=pch, col=col, cex=cex, lwd=lwd)
                
            
        # grid
        if grid:
            for i in range(12):
                r.abline(h=i+.5, lty=3, lwd=1, col='grey')
            for i in range(32):
                r.abline(v=i+.5, lty=3, lwd=1, col='grey')
            
        r.axis(1, at=xL, labels=[str(x) for x in xL], tick=False, line=-1.0, cex_axis=axisSize)
        r.axis(2, at=yL, labels=[str(y) for y in yL], tick=False, line=-1.0, cex_axis=axisSize)
        rdev.close()
        
        return

    
    def plotNumLegend(self, colVecL, breakL, legendName=None, type='png', int_labels=False):

        if legendName is None:
            legendName = 'legend_%i_%i_%i' % (len(colVecL), min(breakL), max(breakL))
            
        full_filename = os.path.join(self.plotDir, legendName)

        #max_rel_cell_count = 0.3
        #min_rel_cell_count = 0.0
        #        
        #pattern = [(0,0,0),(0.7,0,0),(1,1,0),(1,1,1)]
        #colVecL = colors.make_colors(pattern, 500)
        #maxVal = len(colVecL) - 1
        #colBreaksL = [1.0/ maxVal * max_rel_cell_count * x for x in range(len(colVecL) + 1)]

        nb_breaks = 16
        max_break = max(breakL)
        min_break = min(breakL)
        tickBreakL = [float(x) / (nb_breaks - 1) * (max_break - min_break) - min_break for x in range(nb_breaks)]
        if int_labels:
            labels = ['%i' % int(x) for x in tickBreakL]
        else:
            labels = ['%.3f' % x for x in tickBreakL]
        
        rdev = RDevice(name = full_filename, title='', type=type, 
                       width=640, height=120)
        r("par(mar=c(3,1,0,1), las=2)")
        legendA = numpy.reshape(breakL, (len(breakL), 1))
        #legendA = eric.swapaxes(legendA, 1, 0)
        r.image(legendA, 1, legendA, col=colVecL, axes=False, ann=False)
        r.box()
        r.axis(1, at=tickBreakL, labels=labels, tick=True, line=0, cex=0.8, cex_axis=0.8)
        rdev.close()
    
        return
    
    def plotLabelLegend(self, colVecL, labelPhenoD=None, type='png'):
        if labelPhenoD is None:
            labelPhenoD = {0: 'no phenotype', 1: 'mitotic', 2: 'shape', 3: 'mitotic and shape'}

        full_filename = os.path.join(self.plotDir, 'legend_%i' % len(labelPhenoD.keys()))

        labelL = labelPhenoD.keys()
        labelL.sort()
        phenoL = [labelPhenoD[x] for x in labelL]
        
        rdev = RDevice(name = full_filename, title='', type=type, 
                       width=300, height=200)
        r("par(mar=c(0,0,0,0))")
        r.plot([1, 2, 3], type="n", xlab="", ylab="", main="", ann=False, axes=False)
        r.legend(x="center", legend = phenoL, fill=colVecL, bg = "white", bty="n", cex=0.7)
        rdev.close()

    def exportHTMLSingleExperiments(self):
        filename_base = 'localDistribution_'

        # construct ltBaseFilenameD
        ltBaseFilenameD = {}
        labtekBaseL = self.labtekBaseD.keys()
        labtekBaseL.sort()
        for labtekBase in labtekBaseL:
            labtekIdL = filter(lambda(x): x[0:len(labtekBase)] == labtekBase, self.goodL)
            fileNameL = [os.path.join(self.relativePlotDir,'%s%s.png' % (filename_base, labtekId)) for labtekId in labtekIdL]                        
            ltBaseFilenameD[labtekBase] = {'header': labtekIdL, 'filenames': fileNameL}
        
        self.exportHTML(filename='localdistribution.html', legendFilename = 'legend_4.png', labtek_header=True, ltBaseFilenameD=ltBaseFilenameD)
        return

    def exportHTMLMedian(self):
        filename_baseL = ['medianDistribution_', 'finalDistribution_'] 
    
        # construct ltBaseFilenameD
        ltBaseFilenameD = {}
        labtekBaseL = self.labtekBaseD.keys()
        labtekBaseL.sort()
        for labtekBase in labtekBaseL:
            labtekIdL = ['median', 'final hits']
            fileNameL = [os.path.join(self.relativePlotDir,'%s%s.png' % (filename_base, labtekBase)) for filename_base in filename_baseL]                        
            ltBaseFilenameD[labtekBase] = {'header': labtekIdL, 'filenames': fileNameL}
        
        self.exportHTML(filename='mediandistribution.html', legendFilename = 'legend_4.png', labtek_header=True, ltBaseFilenameD=ltBaseFilenameD)
        return

    def exportHTMLPhenoStat(self):
        filename_baseL = ['mean_maxSum_MitosisPhenotype_', 'mean_maxSum_Shape_']
        labtekBaseL = ['druggable', 'genome']
        ltBaseFilenameD = {}
        for labtekBase in labtekBaseL:
            filenameL = [os.path.join(self.relativePlotDir,'%s%s.png' % (x, labtekBase)) for x in filename_baseL]
            ltBaseFilenameD[labtekBase] = {'header' : ['Mitotic Phenotype', 'Shape Phenotype'],
                                           'filenames' : filenameL}
        self.exportHTML(filename='phenostat.html', legendFilename= 'legend_500_00_03.png', labtek_header = True, ltBaseFilenameD=ltBaseFilenameD)
        return

    def exportHTMLPhenoStatNumbers(self):
        filename_baseL = ['numbers_mean_maxSum_MitosisPhenotype_', 'numbers_mean_maxSum_Shape_']
        labtekBaseL = ['druggable', 'genome']
        ltBaseFilenameD = {}
        for labtekBase in labtekBaseL:
            filenameL = [os.path.join(self.relativePlotDir,'%s%s.png' % (x, labtekBase)) for x in filename_baseL]
            ltBaseFilenameD[labtekBase] = {'header' : ['Mitotic Phenotype', 'Shape Phenotype'],
                                           'filenames' : filenameL}
        self.exportHTML(filename='phenostatnumbers.html', legendFilename= None, labtek_header = True, ltBaseFilenameD=ltBaseFilenameD)
        return

    def exportHTMLHitHisto(self):
    
        filename_baseL = ['hitHisto_maxSum_MitosisPhenotype_', 'hitHisto_maxSum_Shape_']
        labtekBaseL = ['druggable', 'genome']
        ltBaseFilenameD = {}
        for labtekBase in labtekBaseL:
            filenameL = [os.path.join(self.relativePlotDir,'%s%s.png' % (x, labtekBase)) for x in filename_baseL]
            ltBaseFilenameD[labtekBase] = {'header' : ['Histogram: Mitotic Phenotype', 'Histogram: Shape Phenotype'],
                                           'filenames' : filenameL}
        self.exportHTML(filename='hithisto.html', legendFilename= 'legend_15_0_16.png', labtek_header = True, ltBaseFilenameD=ltBaseFilenameD)
        return
    
    # for all labteks:
    # pl = meta_analyzer.plotHitsOnLabteks()
    # pl.exportHTML()
    # for base labteks (median):
    # pl.exportHTML(pl.labtekD.keys(), filename_base='medianDistribution_', filename='mediandistribution.html', labtek_header=False)
    # and to join the information of median and definite hits:
    # labtekBaseL = self.labtekBaseD.keys()
    # rowLabtekD = dict(zip(labtekBaseL, zip(['medianDistribution_%s.png' % labtekBase for labtekBase in labtekBaseL],
    #                                       ['finalDistribution_%s.png' % labtekBase for labtekBase in labtekBaseL])))
    # pl.exportHTML(labtekBaseL, filename_base='', filename='localization_median_final_distribution.html', labtek_header=False, rowLabtekD = rowLabtekD)
    def exportHTML(self, filename='localdistribution.html',
                   legendFilename = 'legend_4.png', labtek_header=True, ltBaseFilenameD = None,
                   labtekL=None, filename_base='localDistribution_'):
        if labtekL is None:
            labtekL = self.goodL
            
        file = open(os.path.join(self.htmlDir, filename), 'w')

        # title and header
        titleString="""
<html>        
  <head><title>MitoCheck - Hit Localization on labteks</title></head>
  <body>
    <h2>Hit Localization on Labteks </h2>
<br>        
"""
        file.write(titleString)

        # legend
        if not legendFilename is None:
            legend_full_filename = os.path.join(self.relativePlotDir, legendFilename)
            legendString="""
<p align="left"><img src="%s" border="0"></td>
            """ % legend_full_filename
            file.write(legendString)
        else:
            legendString="""
<p align="left">
<br>
"""
            file.write(legendString)


        # output for all labtek
        if ltBaseFilenameD is None:
            labtekBaseL = self.labtekBaseD.keys()
        else:
            labtekBaseL = ltBaseFilenameD.keys()
        labtekBaseL.sort()

        labtekTableString="""
<table align="left" border="0" cellspacing="10" cellpadding="0">
"""
        file.write(labtekTableString)
        # one row per labtekBase ... for each labtekBase do
        for labtekBase in labtekBaseL:

            if ltBaseFilenameD is None:
                # in this case construct the row from labtekL and filename_base
                rowLabtekL = filter(lambda(x): x[0:len(labtekBase)] == labtekBase, labtekL)
                rowFileNamesL = [os.path.join(self.relativePlotDir,'%s%s.png' % (filename_base, labtekId)) for labtekId in rowLabtekL]
            else:                
                rowLabtekL = ltBaseFilenameD[labtekBase]['header']  
                rowFileNamesL = ltBaseFilenameD[labtekBase]['filenames']
                
            labtekRowString = ""

            # ----- HEADER OF EACH ROW
            # header: name of labteks
            if labtek_header:
                labtekRowString += """
<tr>
<th align="center"> </th>
"""
                for labtekId in rowLabtekL:
                    labtekRowString += """
<th align="center"> %s </th>
""" % labtekId
                labtekRowString += """
</tr>
"""

            # ----- FIRST COLUMN (LABTEKBASE) 
            # write the labtek base (like 'LT0001')
            labtekRowString += """
<tr> <td align="center" valign="middle"> <B> %s   </B> </font> 
""" % labtekBase        

            # ----- OTHER COLUMNS
            # write the labteks for this base:
            for img_filename in rowFileNamesL:
                #img_filename = rowFileNamesD[labtekId]
                #img_filename = os.path.join(self.relativePlotDir,'%s%s.png' % (filename_base, labtekId))
                labtekRowString += """
<td align="center"><img src="%s" border="0"></td>
""" % img_filename

            labtekRowString += """
</tr>
"""
        
            # ----- WRITE EVERYTHING TO THE FILE
            file.write(labtekRowString)

        
        # ----- CLOSE FILE            
        file.close()
        
        return


# ------------------------- #
# --- Neighbor analysis --- #
# ------------------------- #

class NeighborAnalysis(object):
    def __init__(self, metaDir=META_RESULT_DIR):
        self.primaryScreenL = ['LT%04i' % i for i in range(1, 50)] + \
                              ['LT%04i' % i for i in range(60, 160)] 
                            
        self.metaDir = metaDir
        self.subsetShortFilename = os.path.join(self.metaDir, 'hit_list.pickle')        
        return
    
    def readHitList(self, filename=None):
        if filename is None:
            filename = self.subsetShortFilename
        file = open(filename, 'r')
        hitListD = pickle.load(file)
        file.close()
        return hitListD
    
    def getPosition(self, x, y):
        pos = (x - 1) * 12 + (y - 1) + 1
        return pos
        
    def getCoordinates(self, pos): 
        x = (int(pos) - 1) / 12 + 1
        y = (int(pos) - 1) % 12 + 1
        return x, y

    def getNeighborList(self, pos):
        x, y = self.getCoordinates(pos)
        nbL = []
        for delta_x, delta_y in [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]:
            x_nb = x + delta_x
            y_nb = y + delta_y
            if x_nb in range(1, 33) and y_nb in range(1, 13):
                nbL.append(self.getPosition(x_nb, y_nb))
        return nbL

    def getNeighborListForPosList(self, posL):
        nbL = []
        for pos in posL:
            nbPosL = self.getNeighborList(pos)
            for nbPos in nbPosL:
                if nbPos not in nbL:
                    nbL.append(nbPos)
        return nbL
    
    def getSirnaList(self, hitListD):
        sirnaL = []
        for pheno in hitListD.keys():
            for infoD, score in hitListD[pheno]:
                if infoD['sirnaId'] not in sirnaL:
                    sirnaL.append(infoD['sirnaId'])
        return sirnaL

    def getBasePosList(self, hitListD, excludeLabtekL = []):
        expL = []
        for pheno in hitListD.keys():
            for infoD, score in hitListD[pheno]:
                for id, idScore in infoD['expL']:
                    lt, pos = id.split('--')
                    ltBase=lt[0:6]
                    baseId = ltBase + '--' + pos
                    if (not baseId in expL) and (not ltBase in excludeLabtekL):
                        expL.append(baseId)
        return expL

    def getSirnaBasePosDict(self, hitListD, excludeLabtekL = []):
        sirnaD = {}
        for pheno in hitListD.keys():
            for infoD, score in hitListD[pheno]:
                sirnaId = infoD['sirnaId']
                if not sirnaD.has_key(sirnaId):
                    sirnaD[sirnaId] = []
                for id, idScore in infoD['expL']:
                    lt, pos = id.split('--')
                    ltBase=lt[0:6]
                    baseId = ltBase + '--' + pos
                    if (not ltBase in excludeLabtekL) and (not baseId in sirnaD[sirnaId]):
                        sirnaD[sirnaId].append(baseId)
        return sirnaD

    def writeDictToFile(self, aD, filename):
        file = open(filename, 'w')
        for key in aD.keys():
            tempStr = '%s\t%s\n' % (str(key), str(aD[key]))
            file.write(tempStr)
        file.close()
        return
    
    def hitListNeighborAnalysis(self, hitListD):
        phenoL = hitListD.keys()
        idD = {}
        nbD = {}
        for pheno in phenoL:
            idD[pheno] = {}
            nbD[pheno] = {'neighbor': [], 'isolated': []}
            for infoD, score in hitListD[pheno]:
                for id, idScore in infoD['expL']:
                    lt, pos = id.split('--')
                    ltBase=lt[0:6]
                    if not idD[pheno].has_key(ltBase):
                        idD[pheno][ltBase] = []
                    baseId = ltBase + '--' + pos
                    if baseId not in idD[pheno][ltBase]:
                        idD[pheno][ltBase].append(baseId)
            for ltBase in idD[pheno].keys():
                tempL = idD[pheno][ltBase]
                for id in tempL:
                    labtekBase, pos = id.split('--')
                    if labtekBase not in ['LT0601', 'LT0602', 'LT0603']:
                        nbL = ['%s--%s' % (labtekBase, x) for x in self.getNeighborList(pos)]
                        for other in tempL:
                            if other in nbL:
                                nbD[pheno]['neighbor'].append(id)
                                break
                        else:
                            nbD[pheno]['isolated'].append(id)
                
        isolatedL = []
        neighborL = []
        for pheno in phenoL:
            isolatedL = union(nbD[pheno]['isolated'], isolatedL)
            neighborL = union(nbD[pheno]['neighbor'], neighborL)
        isolatedL = setminus(isolatedL, neighborL)
        hitL = union(isolatedL, neighborL)
        sirnaD = self.getSirnaBasePosDict(hitListD, ['LT0601', 'LT0602', 'LT0603'])

        print
        print ' **** '
        for pheno in phenoL:
            print '%s:\t%i' % (pheno, len(union(nbD[pheno]['isolated'], nbD[pheno]['neighbor'])))
        
        for sirnaId in sirnaD.keys():
            if len(sirnaD[sirnaId]) == 0:
                del(sirnaD[sirnaId])
                
        sirnaMoreThanOneSpotL = filter(lambda(x): len(sirnaD[x]) > 1, sirnaD.keys())
        #for sirnaId in sirnaMoreThanOneSpotL:
        #    print sirnaD[sirnaId]
        print 
        print ' **** '
        print 'Number of sirna (not on eg5 neighbor labteks): \t%i' % len(sirnaD.keys()) 
        print 'Number of sirna with more than one position: \t%i' % len(sirnaMoreThanOneSpotL)
        print 'Number of hits: \t%i' % len(hitL)
        
        print 'Number of isolated hits: \t%i \t(%f)' % (len(isolatedL), float(len(isolatedL))/ len(hitL))
        print 'Number of neighbor hits: \t%i \t(%f)' % (len(neighborL), float(len(neighborL))/ len(hitL))
        return nbD

    def getRandomSequence(self, nb_hits, nb_lt, nb_spots_per_labtek=384, controlL = []):
        nb_spots = nb_lt * nb_spots_per_labtek
        absolutePosL = range(1,(nb_spots + 1))
        for controlPos in controlL:
            for i in range(nb_lt):
                remPos = int(controlPos) + i * nb_spots_per_labtek
                absolutePosL.remove(remPos)
                
        randomAbsPosL = random.sample(absolutePosL, nb_hits)
        hitD = {}
        for randomAbsPos in randomAbsPosL:
            lt = randomAbsPos / nb_spots_per_labtek
            pos = randomAbsPos % nb_spots_per_labtek
            ltId = 'LT%04i' % lt
            if not hitD.has_key(ltId):
                hitD[ltId] = []
            hitD[ltId].append('LT%04i--%03i' % (lt, pos))
        return hitD

    def getRandomSequenceFromAllowedPositions(self, nb_hits, absolutePosL):
        randomAbsPosL = random.sample(absolutePosL, nb_hits)
        hitD = {}
        for randomAbsPos in randomAbsPosL:
            ltId, pos = randomAbsPos.split('--')
            if not hitD.has_key(ltId):
                hitD[ltId] = []
            hitD[ltId].append(randomAbsPos)
        return hitD
    
    def getIsolatedAndNeighborHits(self, randomHitD):

        isolatedL = []
        neighborL = []
        for ltBase in randomHitD.keys():
            tempL = randomHitD[ltBase]
            for id in tempL:
                labtekBase, pos = id.split('--')
                nbL = ['%s--%s' % (labtekBase, x) for x in self.getNeighborList(pos)]
                for other in tempL:
                    if other in nbL:
                        neighborL.append(id)
                        break
                else:
                    isolatedL.append(id)
                
        return isolatedL, neighborL
    
    def simulateLocalHitDistribution(self, nbD, controlL=None):
        if controlL is None:
            controlL = ['004', '049', '052', '290', '301', '338', '349', '015', '026', '063',
                        '074', '304', '315', '352', '363', '333', '336', '381', '384', '001']
        phenoL = nbD.keys()
        random_nbD = {}
        for pheno in phenoL:
            random_nbD[pheno] = {}
            nb_hits = len(union(nbD[pheno]['isolated'], nbD[pheno]['neighbor']))
            print '%s:\t%i' % (pheno, nb_hits)
            randomHitD = self.getRandomSequence(nb_hits, nb_lt=150, nb_spots_per_labtek=384, controlL = controlL)
            isolatedL, neighborL = self.getIsolatedAndNeighborHits(randomHitD)
            random_nbD[pheno]['isolated'] = isolatedL
            random_nbD[pheno]['neighbor'] = neighborL

        isolatedL = []
        neighborL = []
        for pheno in phenoL:
            isolatedL = union(isolatedL, random_nbD[pheno]['isolated'])
            neighborL = union(neighborL, random_nbD[pheno]['neighbor'])

        hitL = union(isolatedL, neighborL)

        print ' **** RANDOM DISTRIBUTION **** '
        print 'Number of hits: \t%i' % len(hitL)        
        print 'Number of isolated hits: \t%i \t(%f)' % (len(isolatedL), float(len(isolatedL))/ len(hitL))
        print 'Number of neighbor hits: \t%i \t(%f)' % (len(neighborL), float(len(neighborL))/ len(hitL))
        
        return len(isolatedL), len(neighborL)

    def nRandomRun(self, nbD, nb_runs):
        isoL = []
        nbL = []
        for i in range(nb_runs):            
            iso, nb = self.simulateLocalHitDistribution(nbD)
            isoL.append(iso)
            nbL.append(nb)
        return isoL, nbL

    
    def getBasePosFromSirnaList(self, sirnaL, sirnaD):
        baseIdL = []
        for sirna in sirnaL:
            idL = sirnaD[sirna]['idL']
            for id in idL:
                labtekId, pos = id.split('--')
                baseId = '%s--%s' % (labtekId[:-3], pos)
                if baseId not in baseIdL:
                    baseIdL.append(baseId)
        return baseIdL    
    
    def getIsolatedAndNeighborsForReferenceList(self, testIdL, referenceIdL):
        isoL = []
        nbL = []
        confirmedL = []
        labtekNbD = {}
        
        for baseId in testIdL:
            labtekBase, pos = baseId.split('--')

            if not labtekNbD.has_key(labtekBase):
                labtekNbD[labtekBase] = {'nbL': [], 'isoL': [], 'confirmedL': []}    
            
            # first we check if the hit was automatically found. 
            # if this is the case we do not consider it any further.
            if baseId in referenceIdL:
                confirmedL.append(baseId)
                labtekNbD[labtekBase]['confirmedL'].append(baseId)
                continue
            
            nbPosL = ['%s--%s' % (labtekBase, nbPos) for nbPos in self.getNeighborList(pos)]
            neighbor_flag = False
            for nbPos in nbPosL:
                if nbPos in referenceIdL:
                    neighbor_flag = True
                    continue
            if neighbor_flag:
                # in this case the hit is a neighbor of an automatically found hit
                nbL.append(baseId)
                labtekNbD[labtekBase]['nbL'].append(baseId)
            else: 
                labtekNbD[labtekBase]['isoL'].append(baseId)
                isoL.append(baseId)
        
        return isoL, nbL, confirmedL, labtekNbD

    def getIsolatedAndNeighborsIntrinsic(self, hitPosL):
        isoL = []
        nbL = []
        
        for baseId in hitPosL:
            labtekBase, pos = baseId.split('--')

            nbPosL = ['%s--%s' % (labtekBase, nbPos) for nbPos in self.getNeighborList(pos)]
            neighbor_flag = False
            for nbPos in nbPosL:
                if nbPos in referenceIdL:
                    neighbor_flag = True
                    continue
            if neighbor_flag:
                # in this case the hit is a neighbor of an automatically found hit
                nbL.append(baseId)
            else: 
                isoL.append(baseId)

        return nbL, isoL
                
    def analyseManualAnnotation(self, scoreD, nbD, sirnaD=None):
        incenpPosL = ['001', '004', '049', '052']
        incenpNbL = ['%03i' % x for x in self.getNeighborListForPosList(incenpPosL)]
        
        hitSirnaL = filter(lambda(x): scoreD[x]['pheno'] == 1, scoreD.keys())
        if sirnaD is None:
            spotQC = quality_control.spotQC()
            sirnaD = spotQC.readSirnaView()
        baseIdL = self.getBasePosFromSirnaList(hitSirnaL, sirnaD)

        isolatedL = []
        neighborL = []
        for pheno in nbD.keys():
            isolatedL = union(isolatedL, nbD[pheno]['isolated'])
            neighborL = union(neighborL, nbD[pheno]['neighbor'])
        automaticHitL = union(isolatedL, neighborL)
            
        manIsoL, manNbL, confirmedL, labtekNbD =  self.getIsolatedAndNeighborsForReferenceList(baseIdL, automaticHitL)

        manNotIncenpNbL = []
        manIncenpNbL = []        
        for id in baseIdL:
            labtekId, pos = id.split('--')
            if pos in incenpNbL:
                manIncenpNbL.append(id)
            else:
                manNotIncenpNbL.append(id)
                
        print 'manual hits confirming an automatic hit: %i' % len(confirmedL)
        print 'isolated manual hits: %i' % (len(manIsoL))
        print 'manual hits next to an automatic hit: %i' % (len(manNbL))
        print 'manual hits found in incenp neighborhood: %i' % len(manIncenpNbL)
        print 'total number of siRNA hits: %i' % (len(hitSirnaL))        
        print 'total number of hit positions: %i' % len(baseIdL)
        print 'total number of not manually found position hits: %i' % (len(baseIdL) - len(confirmedL))
        
        print
        
        # randomized distribution
        allowedIdL = self.getAllowedPositions()
        print 'number of allowed spots: %i' % len(allowedIdL)
        allowedIdL = setminus(allowedIdL, automaticHitL)
        print 'number of allowed spots (without automatic hits): %i' % len(allowedIdL)
        
        randomAbsPosL = random.sample(allowedIdL, len(manNbL) + len(manIsoL))
        randIsoL, randNbL, randConfirmedL, randNbD =  self.getIsolatedAndNeighborsForReferenceList(randomAbsPosL, automaticHitL)

        randNotIncenpNbL = []
        randIncenpNbL = []
        for id in randomAbsPosL:
            labtekId, pos = id.split('--')
            if pos in incenpNbL:
                randIncenpNbL.append(id)
            else:
                randNotIncenpNbL.append(id)                

        print
        print ' ****************** '
        print ' RANDOMIZED '
        print 'random hits confirming an automatic hit: %i' % len(randConfirmedL)
        print 'isolated random hits: %i' % (len(randIsoL))
        print 'random hits next to an automatic hit: %i' % (len(randNbL))
        print 'manual hits found in incenp neighborhood: %i' % len(randIncenpNbL)
        print 'total number of random hits %i' % (len(randomAbsPosL))
        print
        print
        
        # the same for both, manual and automatic hits
        allHitIdL = union(automaticHitL, baseIdL)
        
        
        return labtekNbD
    
    def getAllowedPositions(self):
        self.lc = spotInfo.LabtekConfiguration()
        idL = self.lc.getIdListFromTypeList(self.primaryScreenL, ['experiment', 'bcopNeighbors', 'incenpNeighbors', 'eg5Neighbors'])
        return idL
    
    def __call__(self, hitListD=None, scoreD=None):
        
        if hitListD is None:
            hitListD = self.readHitList(filename)
        
        if scoreD is None:
            man = manual_annotation.manualAnnotationObject()
            manualAnnotationD = man.readAnnotationScore()
            scoreD = man.getSirnaScore(manualAnnotationD, penetranceThreshold=1)
        
        nbD = self.hitListNeighborAnalysis(hitListD)
        iso, nb = self.simulateLocalHitDistribution(nbD)
        
        print 
        print 
        print ' *** '
        print 'WITH MANUAL ANNOTATION'
        labtekNbD = self.analyseManualAnnotation(scoreD, nbD)
        
        return
    