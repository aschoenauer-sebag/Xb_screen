
import classification
import os
import spotAnalyzer
from rpy import r
import copy
import plot_utilities
import numpy 
import operator

# import colors from r (for heat maps)
import colors

from optparse import OptionParser


#GROUP = 'group1'
#LABTEK = 'LTValidmitosis7x3course2009_01'
RESULT_DIR = '/netshare/almf/EMBO-data/timelapse_result/group1'
LOCAL_PLOT_DIR = '/netshare/almf/EMBO-data/timelapse_result/group1/html'
HTML_DIR = LOCAL_PLOT_DIR
NB_ROW = 8
NB_COL = 12

class LabtekConfiguration(object):
    def __init__(self):
        self.labtekL = ['LTValidMitosis7x3course2009']
        self.negCtrPosL = ['002', '008', '035', '051', '058', '073', '094']
        self.emptyL = ['093', '095']
        self.incenpL = ['047', '069']
        self.pirD = classification.plateIdReader().read()
        return
    
    def getSpotConfForLabtek(self, labtekBase=None):
        if labtekBase is None:
            labtekBase = self.labtekL[0]
        labtekBase = labtekBase.split('_')[0][2:]
        spotConfD = {}
        for pos in range(1, 97):
            infoD = self.pirD[labtekBase]['%03i' % pos]
            if pos in self.negCtrPosL:
                infoD['type'] = 'scrambled'
            else:
                infoD['type'] = 'experiment'  
            spotConfD['%03i' % pos] = infoD
        return spotConfD
    
    def getSirnaView(self, labtekId, spotConfD=None):
        labtekBase = labtekId.split('_')[0]
        if spotConfD is None:
            spotConfD = self.getSpotConfForLabtek(labtekBase)

        sirnaD = {}
        for pos in spotConfD.keys():
            id = '%s--%s' % (labtekId, pos)
            sirna = spotConfD[pos]['sirnaId']
            gene = spotConfD[pos]['geneName']
            
            if not sirnaD.has_key(sirna):
                sirnaD[sirna] = { 'idL': [], 'gene': gene}
            sirnaD[sirna]['idL'].append(id)
            
        return sirnaD
        
class BasicDescriptionPlate(object):
    def __init__(self, baseDir=RESULT_DIR):
        self.baseDir = baseDir
        self.labtekL = filter(lambda x: x[0:2] == 'LT', os.listdir(self.baseDir))
        return
    
    def __call__(self, labtekL=None):
        
        if labtekL is None:
            labtekL = self.labtekL
            
        print 'BASE_DIR: ', self.baseDir
        resD = {}
        for labtek in labtekL:
            spotConfD = LabtekConfiguration().getSpotConfForLabtek(labtek)
            posL = spotConfD.keys()
            for pos in posL:
                sirna = spotConfD[pos]['sirnaId']
                gene = spotConfD[pos]['geneName']
                try:
                    spotRes = spotAnalyzer.spotAnalyzer(labtek, pos, 
                        baseDir=self.baseDir, 
                        spotConfD=spotConfD, prediction_suffix='_prediction_track.dat')
                    phenoSpotD, countL = spotRes.readPrediction()
                    id = '%s--%s' % (labtek, pos)
                    resD[id] = {'initCellCount': countL[0], 'endCellCount': countL[-1], 
                                'proliferation': float(countL[-1]) / float(countL[0]),
                                'sirna': sirna, 'gene': gene}
                except:
                    resD[id] = {'initCellCount': 0, 'endCellCount': 0, 
                                'proliferation': 0.0, 
                                'sirna':sirna, 'gene': gene}
        return resD
    
    #LTValidMitosis7x3course2009_01/results/068
    
    
class HTMLGenerator(object):
    def __init__(self, html_dir=HTML_DIR, baseDir=RESULT_DIR):
        
        self.html_dir = html_dir
        self.baseDir = baseDir
        self.ap = ArrayPlotter()
        self.lc = LabtekConfiguration()
        self.ap = ArrayPlotter(plotDir=os.path.join(self.html_dir, 'plots'), 
                               legendDir=os.path.join(self.html_dir, 'plots'))
        return
        

    # plots in general a bundle of time series
    def plotBundle(self, bundleD, full_filename, colorsD=None, bundlePointsD=None, legendL=None, title=None, y_max=None):

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

        try:
            r.png(full_filename, width=800, height=600)
            oldPar = r.par(xpd = True, mar = [x + y for (x,y) in zip(r.par()['mar'], [0,0,0,6])])
        
            print 'plot %s' % full_filename
            r.plot(timeVec, timeVec,
                   type='n',
                   main=title, ylim=(0, y_max),
                   xlab="time in hours after transfection", ylab="Relative Cell Counts",
                   pch=20, lwd=1, lty = 1, 
                   cex=1.0, cex_lab=1.2, cex_main=1.5)
        
        
            for bundleId in bundleIdL:
                
                if not bundlePointsD is None:
                    r.points(timeVec, bundlePointsD[bundleId],
                             col=colorsD[bundleId], pch=20,
                             lwd=1)
                    r.lines(timeVec, bundlePointsD[bundleId],
                            col=colorsD[bundleId],
                            lwd=1, lty = 1)

                r.lines(timeVec, bundleD[bundleId],
                        col=colorsD[bundleId],
                        lwd=3, lty = 1)

            r.legend(max(timeVec) * 1.1, y_max, legend=legendL, fill=colorsL, cex=1.0, bg= 'whitesmoke')
            r.par(oldPar)
            r.grid(col="darkgrey")
        
            r.dev_off() 
        except:
            r.dev_off()
            print full_filename + ' has not been printed.'
            

        return

    # plot the controls according to settings file.
    def generateControlPlots(self, labtekId, control='scrambled'):
        plotDir = os.path.join(self.baseDir, labtekId, 'plots', 'controls')
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)
        
        if control == 'scrambled':
            ctrlPosL = self.lc.negCtrPosL
        elif control== 'empty':
            ctrlPosL = self.lc.emptyL            
        elif control=='incenp':
            ctrlPosL = self.lc.incenpL
        else:
            return
            
        controlDataD = {}
        controlDataRawD = {}
        controlProliferationD = {}
        controlProliferationFitD = {}

        spotConfD = LabtekConfiguration().getSpotConfForLabtek(labtekId)
        
        control_len = None
        for pos in ctrlPosL:

            try:
                id = '%s--%s' % (labtekId, pos)
                spotRes = spotAnalyzer.spotAnalyzer(labtekId, pos, baseDir=self.baseDir, 
                    spotConfD=spotConfD, prediction_suffix='_prediction_track.dat')
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

    def generateDensityPlots(self, labtekId, resD):
        
        labeledPosD = {'s': self.lc.negCtrPosL,
                       '-': self.lc.emptyL}
    
        idL = filter(lambda x: x.split('--')[0] == labtekId, resD.keys())
        idL.sort()
        
        pattern = [(0,0,0), (0,0,0.5), (0,1,1), (1,1,1)]
        colVecL = colors.make_colors(pattern, 500)
        settingL = [
                    {'min': 15, 'max': 800, 'colVecL': None},
                    {'min': 15, 'max': 800, 'colVecL': None},
                    {'min': 0.0, 'max': 4.0, 'colVecL': colVecL}
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

    
    def __call__(self, labtekId):
        print ' *** BasicDescriptionPlate'
        bdp = BasicDescriptionPlate(baseDir=self.baseDir)
        print ' *** get resD'
        resD = bdp([labtekId])
        print ' *** generate density plots'
        self.generateDensityPlots(labtekId, resD)
        print ' *** generate scrambled plots'
        self.generateControlPlots(labtekId, 'scrambled')
        print ' *** generate empty plots'
        self.generateControlPlots(labtekId, 'empty')
        print ' *** generate incenp plots'
        self.generateControlPlots(labtekId, 'incenp')
        print ' *** generate qc page'
        self.generateQCPage(labtekId, resD)        
        return
        
    def generateQCPage(self, labtekId, resD):
        plotDir = os.path.join(self.html_dir, labtekId)
        relativePlotDir = os.path.join('./plots', labtekId)
        
        pheno = 'initCellCount'
        relativeMovieDir = os.path.join('..', labtekId, 'movies')
        absoluteMovieDir = os.path.join(self.html_dir, '..', labtekId, 'movies')
        idL = filter(lambda x: x[0:len(labtekId)] == labtekId, resD.keys())
        movieFileL = os.listdir(absoluteMovieDir)
        print absoluteMovieDir
        print relativeMovieDir
        print idL
        print 'movieFileL: ', movieFileL
        for id in idL:
            try:
                movieName = filter(lambda x: x[0:len(id)].lower() == id.lower(), movieFileL)[0]
            except:
                print 'no movie for ', id
                print filter(lambda x: x[0:len(id)] == id, movieFileL)
                
                resD[id]['movieName'] = 'XXX'
                continue
            resD[id]['movieName'] = os.path.join(relativeMovieDir, movieName)
        imageName = os.path.join(relativePlotDir, '%s--%s.png' % (pheno, labtekId))
        mapStr = {pheno: self.makeImageMapHTML(labtekId, resD, imageName, pheno)}
        
        for pheno in ['endCellCount', 'proliferation']:
            imageName = os.path.join(relativePlotDir, '%s--%s.png' % (pheno, labtekId))
            mapStr[pheno] = self.makeImageMapHTML(labtekId, resD, imageName, pheno)
            
        mapD = {(os.path.join(relativePlotDir, '%s--%s.png') % ('initCellCount', labtekId)): mapStr['initCellCount'],
                (os.path.join(relativePlotDir, '%s--%s.png') % ('endCellCount', labtekId)): mapStr['endCellCount'],
                (os.path.join(relativePlotDir, '%s--%s.png') % ('proliferation', labtekId)): mapStr['proliferation']
            }
        #mapD = {(os.path.join(relativePlotDir, '%s--%s.png') % ('initCellCount', labtekId)): mapStr['initCellCount']}
        rowL = ['Overview', 'Scrambled', 'Incenp']        
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
        relativePlotDir = os.path.join('..', labtekId, 'plots', 'controls')
        rowD['Scrambled'] = {'header': None,
                             'images': [os.path.join(relativePlotDir, 'scrambled_%s.png' % pheno) for pheno in 
                                       ['proliferation', 'Apoptosis', 'Prometaphase', 'MetaphaseAlignment', 'Shape1', 'Shape3']],
                             'height': [500, 500, 500, 500, 500, 500]}
        rowD['Empty'] = {'header': None,
                         'images': [os.path.join(relativePlotDir, 'empty_%s.png' % pheno) for pheno in 
                                   ['proliferation', 'Apoptosis', 'Prometaphase', 'MetaphaseAlignment', 'Shape1', 'Shape3']],
                         'height': [500, 500, 500, 500, 500, 500]}

        rowD['Incenp'] = {'header': None,
                         'images': [os.path.join(relativePlotDir, 'incenp_%s.png' % pheno) for pheno in 
                                   ['proliferation', 'Apoptosis', 'Prometaphase', 'MetaphaseAlignment', 'Shape1', 'Shape3']],
                         'height': [500, 500, 500, 500, 500, 500]}
                                       
#        rowD['Incenp'] = {'header': None, #['Proliferation', 'Shape 3', 'Prometaphase', 'Apoptotic'],
#                          'images': [os.path.join(relativePlotDir, 'incenp_%s.png') % pheno
#                                     for pheno in ['proliferation', 'Shape3', 'Shape1', 'Prometaphase', 'MetaphaseAlignment', 'Apoptosis']] }
#        rowD['Scrambled'] = {'header': None, #['Proliferation', 'Shape', 'Mitotic', 'Apoptotic'],
#                             'images': [os.path.join(relativePlotDir, 'scrambled_%s.png') % pheno
#                                        for pheno in ['proliferation', 'Shape3', 'Shape1', 'Prometaphase', 'MetaphaseAlignment', 'Apoptosis']] }

        # generate one row per gene/siRNA
        sirnaD = self.lc.getSirnaView(labtekId)
        
            
        sirnaL = sirnaD.keys()
        geneL = [sirnaD[sirna]['gene'] for sirna in sirnaL]
        expCondL = zip(sirnaL, geneL)
        expCondL.sort(key=operator.itemgetter(-1))
        relativePlotDir = os.path.join('..', labtekId, 'plots')
        for sirna, gene in expCondL:
            if len(sirnaD[sirna]['idL']) != 1:
                print sirna, gene, sirnaD[sirna]['idL']
                continue
                
            id = sirnaD[sirna]['idL'][0]
            lt_id, pos = id.split('--')
            current_plot_dir = os.path.join(self.html_dir, '..', labtekId, 'plots' , pos)
            try:
                pheno_file = filter(lambda x: x[0:len('PhenoPlotAnalysis')] == 'PhenoPlotAnalysis', os.listdir(current_plot_dir))[0]
                prol_file = filter(lambda x: x[0:len('ProliferationAnalysis')] == 'ProliferationAnalysis', os.listdir(current_plot_dir))[0]
    
                current_plot_dir = os.path.join(relativePlotDir, pos)         
                exp_cond_id = '%s (%s)' % (gene, sirna)
                
                movieName = resD[id]['movieName']
                moviePageName = os.path.join('.', labtekId, os.path.basename(movieName).replace('.avi', '.html'))
                rowD[exp_cond_id] = {'header': None,
                                     'images': [os.path.join(current_plot_dir, pheno_file), os.path.join(current_plot_dir, prol_file)],
                                     'height': [500, 500],
                                     'header_link': moviePageName
                                     }
                rowL.append(exp_cond_id)
            except:
                print 'missing condition: ', exp_cond_id
                
        self.exportPlateHTML('Summary--%s.html' % labtekId, labtekId, rowD, rowL)
        return

    def makeMovieLinkPage(self, labtekId, movieName, gene, sirna):
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
    def makeImageMapHTML(self, labtekId, resD, imageName, pheno='initCellCount'):        
        
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
                movieName = resD[id]['movieName']
                moviePageName = os.path.join('.', labtekId, os.path.basename(movieName).replace('.avi', '.html'))
                gene = resD[id]['gene']
                sirna = resD[id]['sirna']
                self.makeMovieLinkPage(labtekId, movieName, gene, sirna)
                
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
<tr> <td align="center" valign="middle"> <B> <font size=5> %s </font>   </B> </font> 
""" % rowTitle  
            else:

                labtekRowString += """
<tr> <td align="center" valign="middle"> <B> <font size=5> <A HREF="%s" target="_blank"> %s </A>
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
    def __init__(self, plotDir = os.path.join(LOCAL_PLOT_DIR, 'array'), legendDir = os.path.join(LOCAL_PLOT_DIR, 'legend'),
                 nb_row = NB_ROW, nb_col=NB_COL):
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
            return max(valuesD[self.phenoClass], 0)

    def prepareData(self, resD, getLabel, expL):
        
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
        hitData = r.matrix(hitValueL, nrow=self.nb_row, byrow=True)
        
        return hitData    


    def plotArray(self, hitDataParam, filename, plotDir=None, labeledPosD=None,
                  title='', type='png',
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
        nrow = self.nb_row
        ncol = self.nb_col
        
        xL = range(1, self.nb_col + 1)
        yL = range(1, self.nb_row + 1)
        
        rdev = plot_utilities.RDevice(name = full_filename, title=title, plotType=type, 
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
                self.plotLabelLegend(colVecL, plotType=type, filename=legend_filename)
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

            if colVecL is None:                    
                pattern = [(0,0,0),(0.7,0,0),(1,1,0),(1,1,1)]
                colVecL = colors.make_colors(pattern, nb_colors)
            colBreaksL = [1.0/ (len(colVecL) - 1) * x * (max_rel_cell_count - min_rel_cell_count) + min_rel_cell_count  for x in range(len(colVecL) + 1)]

            if legend_plot:
                self.plotNumLegend(colVecL, colBreaksL, 16, filename=legend_filename, type=type, int_labels=numeric_as_int, legendDir = plotDir)

        axisSize = .8
        r("par(mar=c(1.6,1.6,0.1,0.1))")
        r.image(xL, yL, r.t(hitData), axes = False, ann=False, cex=1, col=colVecL, breaks=colBreaksL)
        r.box()
        if not labeledPosD is None:
            for label in labeledPosD.keys():
                posL = labeledPosD[label]
                if len(posL) > 0:
                    xlL = [(int(x)-1) % self.nb_col + 1 for x in posL]
                    ylL = [(int(x)-1) / self.nb_col + 1 for x in posL]
                    r.points(xlL, ylL, pch=label, col=control_color, cex=axisSize)
                    print
                    print xlL
                    print ylL
                    print
        # grid
        if grid:
            for i in range(self.nb_col):
                r.abline(h=i+.5, lty=3, lwd=1, col='grey')
            for i in range(self.nb_row):
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
        
        rdev = plot_utilities.RDevice(name = full_filename, title='', plotType=type, 
            width=640, height=120)
        r("par(mar=c(3,1,0,1), las=2)")
        #legendA = rpy.reshape(breakL, (len(breakL), 1))
        legendA = r.matrix(breakL, byrow=False, ncol=1)
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
        
        rdev = plot_utilities.RDevice(name = full_filename, title='', plotType=type, 
            width=300, height=200)
        r("par(mar=c(0,0,0,0))")
        r.plot([1, 2, 3], type="n", xlab="", ylab="", main="", ann=False, axes=False)
        r.legend(x="center", legend = phenoL, fill=colVecL, bg = "white", bty="n", cex=0.7)
        rdev.close()
    
    
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-g", "--group", dest="group",
                      help="group")
    parser.add_option("-l", "--labtek", dest="labtek",
                      help="labtek")
                      

    (options, args) = parser.parse_args()        
    group = options.group
    labtek = options.labtek
    
    LOCAL_PLOT_DIR = os.path.join('/netshare/almf/EMBO-data/timelapse_result', group, 'html')
    HTML_DIR = LOCAL_PLOT_DIR
    RESULT_DIR = os.path.join('/netshare/almf/EMBO-data/timelapse_result', group)
    
    ht = HTMLGenerator(html_dir=HTML_DIR, baseDir=RESULT_DIR)
    ht(labtek)
    
    