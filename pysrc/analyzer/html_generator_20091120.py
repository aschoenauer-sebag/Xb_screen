from util import *
import os
import colors
from rpy import r
import plot_utilities
import copy
import re
import time

r.library("RColorBrewer")

RESULT_FILENAME = '/home/twalter/data/JoanaDruggableResult.csv'
HTML_DIR = '/netshare/mitofiler/Thomas/Joana_HTML/html'
IMAGE_BASE_DIR = '/netshare/mitofiler/Thomas/Joana'
#PICKLERESULTDIR = '/netshare/mitofiler/Thomas/Joana_HTML'

LOCAL_PLOT_DIR = os.path.join(HTML_DIR, 'plots')
NB_ROW = 32
NB_COL = 12


def get_color_brewer_pattern(brew_str):
    hex_pattern = r.brewer_pal(9, brew_str)
    #def hex2int(hexStr): return([int(hexStr[1:3], 16), int(hexStr[3:5], 16), int(hexStr[5:], 16)])
    pattern = [(int(hexStr[1:3], 16) / 255.0, int(hexStr[3:5], 16) / 255.0, int(hexStr[5:], 16) / 255.0) for hexStr in hex_pattern]
    return pattern
    
def convert_numeric_table_to_res(filename=RESULT_FILENAME):
    tabD = readTableFromFile(filename, sep=';', header=True)
    lt_idL = ['%s_%s' % (x[2:], y) for x,y in zip(tabD['Labtek'], tabD['Date']) ]
    idL = ['%s--%03i' % (x, int(y)) for x, y in zip(lt_idL, tabD['Spot'])]
    
    ratioL = [float(x.replace(',', '.')) for x in  tabD['Ratio']]
    ratioCorrL = [x if (x==x) else 0.0 for x in ratioL]
    numberL = [int(x) for x in tabD['Number']]
    faL = [float(x.replace(',', '.')) for x in  tabD['FractionIndex']]
    faCorrL = [ x if (x==x) else 0.0 for x in faL]
    
    sirnaL = tabD['siRNAID']
    geneL = tabD['Gene']
    
    resD = {}
    for i in range(len(idL)):
        id = idL[i]
        resD[id] = {'gene': geneL[i],
                    'sirna': sirnaL[i],
                    'pos': int(tabD['Spot'][i]),
                    'number': numberL[i],
                    'ratio': ratioCorrL[i],
                    'fractionIndex': faCorrL[i]
                    }

    return resD
    
def add_new_values(resD, result_filename):
    tabD = readTableFromFile(result_filename, header=True)
    idL = tabD['ID']
    val_keyL = tabD.keys()
    val_keyL.remove('ID')
    
    #missingIdL = filter(lambda x: x not in idL, resD.keys())
    newIdL = []
    for i in range(len(idL)):
        id = idL[i]
        
        if resD.has_key(id):
            for key in val_keyL:
                if tabD[key][i] == 'NA':
                    resD[id][key] = 0.0
                else:
                    resD[id][key] = float(tabD[key][i])
        else:
            newIdL.append(id)
            
    if len(newIdL) > 0:
        print 'There were %i ids that were not in the initial resD structure' % len(newIdL)
#    if len(missingIdL) > 0:
#        print 'There were %i ids that were in the initial resD structure and for which there is no new data' % len(missingIdL)
        
    return #{'new': newIdL, 'missing': missingIdL}
    
def add_raw_data_info(resD, imageBaseDir=IMAGE_BASE_DIR):
    idL = resD.keys()
    idL.sort()
    
    labtekL = set([id.split('--')[0] for id in idL])
    for labtek in labtekL:
        print labtek
        subsetIdL = filter(lambda x: x.split('--')[0] == labtek, idL)
        imageDir = os.path.join(imageBaseDir, labtek, 'data')
        imageNameL = filter(lambda x: x.split('.')[-1][0:3].lower() == 'jpg', os.listdir(imageDir))
        #posImageD = dict(zip([x.split('--')[1][1:] for x in imageNameL], imageNameL))
        for id in subsetIdL:
            pos = id.split('--')[-1]
            matchingImageNameL = filter(lambda x: not re.search('W0*%s' % pos, x) is None, imageNameL)
            img1 = filter(lambda x: not re.search('T0*0', x) is None, matchingImageNameL)[0]
            img2 = filter(lambda x: not re.search('T0*1', x) is None, matchingImageNameL)[0]
            resD[id]['img1'] = img1
            resD[id]['img2'] = img2
        
    return 

def add_type_info(resD):
    controlPosD = {'scrambled': ['015', '026', '063', '074', '304', '315', '352'],
                   'copb' : ['363'],
                   'empty' : ['290', '301', '338', '349'],
                   'incenp' : ['004', '049', '052'],
                   'marker' : ['001']
                   }
    posL = ['%03i' % (i + 1) for i in range(384)]
    posD = {}
    for pos in posL:
        posD[pos] = 'experiment'
        for control in controlPosD.keys():
            if pos in controlPosD[control]:
                posD[pos] = control
                
    for id in resD.keys():
        labtekId, pos = id.split('--')
        resD[id]['type'] = posD[pos]
             
    return
    
class HTMLGenerator(object):
    def __init__(self, html_dir=HTML_DIR):
        
        self.html_dir = html_dir
        self.ap = ArrayPlotter()
        self.ap = ArrayPlotter(plotDir=os.path.join(self.html_dir, 'plots'), 
                               legendDir=os.path.join(self.html_dir, 'plots'))
        return
        

    def generateScreenPage(self, resD):
        
        
        phenoL = ['number', 'ratio', 'fractionIndex', 'fractionIndexPlateThresh']
        
        file = open(os.path.join(self.html_dir, 'index.html'), 'w')
        platePageL = filter(lambda x: x[:len('Summary')] == 'Summary', os.listdir(self.html_dir))
        platePageL.sort()
        
        nb_cols = 6
        nb_rows = len(platePageL) / 6
        
        geneL = set([resD[x]['gene'] for x in resD.keys()])
        sirnaL = set([resD[x]['sirna'] for x in resD.keys()])
        expL = filter(lambda x: resD[x]['type'] == 'experiment', resD.keys())
         
        # start of the page
        screen_page_string = """
<html>        
  <head><title>Screen Plate Overview</title></head>
  <body>    
    <h2> Summary </h2>
    
    Number of processed labteks: %i
    <br>
    Number of spots (including controls): %i
    <br>
    Number of experiments: %i 
    <br>
    Number of siRNAs: %i
    <br>
    Number of genes: %i
    <br>
    Plate Overview: number, ratio, fraction index
    <br>
    images converted to 8bit-jpg
    <br>
    <h2> Screen Plot Overview </h2>
    <table align="left" border="0" cellspacing="4" cellpadding="0">
""" % (len(platePageL), (len(platePageL) * 384), len(expL), len(sirnaL), len(geneL) )
        # labtek entries
        i = 0
        
#'images': [os.path.join(relativePlotDir, '%s--%s.png') % (pheno, labtekId)
#                                    for pheno in phenoL]
        
        screen_page_string += """         
    <tr>
""" 

        for header in ['LabtekId'] + phenoL:

            screen_page_string += """         
        <th align="center" valign="middle"> <font size=3> %s
        </td>
""" % header

        screen_page_string += """         
    </tr>
""" 
            
        for plate_page in platePageL:
            
            i += 1
            labtekId = os.path.splitext(os.path.basename(plate_page))[0].split('--')[-1]
            screen_page_string += """
    <tr> <td align="left" valign="middle"> <font size=2> <A HREF="%s" target="_blank"> %s (%i/%i) </A> </td>
""" % (plate_page, labtekId, i,len(platePageL))

            relativePlotDir = os.path.join('./plots', labtekId)
            for pheno in phenoL:
                screen_page_string += """         
        <td align="left" valign="middle"> <img src=%s height=60 border="0" vspace="0" style="image-orientation: 90deg;">
        </td>
""" % (os.path.join(relativePlotDir, 'Rotate--%s--%s.png') % (pheno, labtekId))


            screen_page_string += """         
    </tr>
""" 
        
        # end of the page
        screen_page_string += """
    </table>
  </body>
</html>
"""
        file.write(screen_page_string)
        file.close()
        return
        
    def __call__(self, resD, labtekL=None):
        if labtekL is None:
            labtekL = list(set([x.split('--')[0] for x in resD.keys()]))

        labtekL.sort()
        i = 0
        startTime = time.time()
        badLabtekL = []
        for labtek in labtekL:
            i += 1
            try:
                print
                print 'generating plots for %s' % labtek
                self.generateDensityPlots(labtek, resD)
                print 'generating plate page for %s' % labtek
                self.generateQCPage(labtek, resD)
                diffTime = time.time() - startTime
                print 'finished %i / %i (%02i:%02i:%02i)' % (i, len(labtekL), (diffTime/3600), ((diffTime%3600)/60), (diffTime%60))
                startTime = time.time()
            except:
                print ' *** PROBLEM WITH LABTEK: %s' % labtek
                badLabtekL.append(labtek)
        
        print 'problems with the following labteks: ' 
        print badLabtekL
        
        print
        print 'generating screen page' 
        self.generateScreenPage(resD)
        
        print 'DONE!'
        
        return
        
    def generateDensityPlots(self, labtekId, resD):
        
        labeledPosD = {'s': ['015', '026', '063', '074', '304', '315', '352'],
                       'b' : ['363'],
                       'e' : ['290', '301', '338', '349'],
                       'i' : ['004', '049', '052'],
                       'm' : ['001']
                       }      
    
        idL = filter(lambda x: x.split('--')[0] == labtekId, resD.keys())
        idL.sort()
        
        #pattern = [(0,0,0), (0,0,1), (0,0.5,1), (1,1,1)]
#        pattern = [(0.00, 0.10, 0.30)
#                   (0.00, 0.22, 0.61),
#                   (0.1, 0.31, 0.74), 
#                   (0.3, 0.48, 0.84),
#                   (0.6, 0.64, 0.91),
#                   (0.94, 0.95, 1.00)]
        pattern = get_color_brewer_pattern("YlGnBu")
        col2 = colors.make_colors(pattern, 500)
        
#        pattern = [(0, 0, 0), 
#                   (0.0, 0.26, 0.1), 
#                   (0.0, 0.43, 0.17), 
#                   (0.1, 0.54, 0.27), 
#                   (0.2, 0.67, 0.36), 
#                   (0.3, 0.77, 0.46), 
#                   (0.4, 0.85, 0.60), 
#                   (0.6, 0.91, 0.75), 
#                   (0.8, 0.96, 0.88), 
#                   (0.97, 0.98, 1.0)]        
#        pattern = [(0,0,0), (0,0.5,0), (0.3,1,0), (1,1,1)]
        pattern = get_color_brewer_pattern("YlGn")
        col1 = colors.make_colors(pattern, 500)

        #paletteName = 'RdPu'
        #colRampFunc = r.colorRampPalette(brewer.pal(9, paletteName))
        #colVec = colRampFunc(384)

        settingL = [
                    {'min': 0.0, 'max': 80.0, 'colVecL': col1},
                    {'min': 0.0, 'max': 100.0, 'colVecL': col2},
                    {'min': 0, 'max': 800, 'colVecL': None},
                    {'min': 0, 'max': 384, 'colVecL': None},
                    {'min': 0.0, 'max': 100.0, 'colVecL': col2},                   
                    ]
                    
        # init cell count
        for pheno, setting in zip(['ratio', 'fractionIndex', 'number', 'pos', 'fractionIndexPlateThresh'], settingL):
            plotDir = os.path.join(self.html_dir, 'plots', labtekId)
            if not os.path.isdir(plotDir):
                os.makedirs(plotDir)
            filename = '%s--%s' % (pheno, labtekId)
            data = self.ap.prepareData(resD, self.ap.getPhenoVal(pheno), idL)
            self.ap.plotArray(data, filename, plotDir=plotDir, colVecL=setting['colVecL'],
                min_val=setting['min'], max_val=setting['max'], 
                label_colors = False, labeledPosD=labeledPosD,
                legend_plot = True, legend_filename='legend_%s' % filename)

            filename = 'Rotate--%s--%s' % (pheno, labtekId)
            data = self.ap.prepareData(resD, self.ap.getPhenoVal(pheno), idL, flip_axis=True)
            self.ap.plotArray(data, filename, plotDir=plotDir, colVecL=setting['colVecL'],
                min_val=setting['min'], max_val=setting['max'], 
                label_colors = False, labeledPosD={},
                legend_plot = False, flip_axis=True, axis=False)
        
        return

            
    def generateQCPage(self, labtekId, resD):
        plotDir = os.path.join(self.html_dir, labtekId)
        relativePlotDir = os.path.join('./plots', labtekId)
        
        phenoL = ['number', 'ratio', 'fractionIndex', 'fractionIndexPlateThresh']       
        mapStr = {} 
        for i in range(len(phenoL)):
            pheno = phenoL[i]
            imageName = os.path.join(relativePlotDir, '%s--%s.png' % (pheno, labtekId))
            mapStr[pheno] = self.makeImageMapHTML(labtekId, resD, imageName, pheno, i)
            
        mapD = {(os.path.join(relativePlotDir, '%s--%s.png') % ('number', labtekId)): mapStr['number'],
                (os.path.join(relativePlotDir, '%s--%s.png') % ('ratio', labtekId)): mapStr['ratio'],
                (os.path.join(relativePlotDir, '%s--%s.png') % ('fractionIndex', labtekId)): mapStr['fractionIndex'],
                (os.path.join(relativePlotDir, '%s--%s.png') % ('fractionIndexPlateThresh', labtekId)): mapStr['fractionIndexPlateThresh']
                }

        #rowL = ['Overview', 'Legend']        
        rowD = {}        
        rowD['Plate'] = {'header': ['Cell Count', 'Ratio', 'Fraction Index', 'Fraction Index (Plate Thresh)'],
                         'images': [os.path.join(relativePlotDir, '%s--%s.png') % (pheno, labtekId)
                                    for pheno in phenoL],
                         'map': mapD}
        rowD['Legend'] = {'header': None,
                          'images': [os.path.join(relativePlotDir, 'legend_%s--%s.png') % (pheno, labtekId)
                                     for pheno in phenoL],
                          'height': [80 for x in phenoL]}
                
        self.exportPlateViewHTML('Summary--%s.html' % labtekId, labtekId, rowD)
        
        return

    def makeRawDataLinkPage(self, labtekId, img1, img2, gene, sirna, id, infoD=None):
        rawdata_html_dir = os.path.join(self.html_dir, labtekId)
        if not os.path.isdir(rawdata_html_dir):
            os.makedirs(rawdata_html_dir)
            
        filename = os.path.join(rawdata_html_dir, '%s.html' % id)
        html_file = open(filename, 'w')
        #print 'generating ', filename

        relative_img_dir = os.path.join('..', '..', labtekId, 'data')
                
        infoStr1 = "%s<br>%s<br>%s<br>time:0" % (gene, sirna, id) 
        infoStr2 = "%s<br>%s<br>%s<br>time:1" % (gene, sirna, id)         
                        
        if not infoD is None:
            for key, val in infoD.iteritems():
                infoStr1 += '<br>%s:%.2f' % (key, val)
            for key, val in infoD.iteritems():
                infoStr2 += '<br>%s:%.2f' % (key, val)
                
        # HERE I AM
        tempStr = """
<html>        
  <head><title>Raw images: %s %s %s</title></head>
  <body>    
    <table align="left" border="0" cellspacing="10" cellpadding="0">
    <tr align="left">
        <td> %s </td>
        <td> <img src="%s" border="0" width=800 vspace=10> </td>
    </tr>
    <tr align="left">
        <td> %s </td>
        <td> <img src="%s" border="0" width=800 vspace=10> </td>
    </tr>
    </table>
  </body>
</html>
""" % (id, gene, sirna, infoStr1, os.path.join(relative_img_dir, img1), infoStr2, os.path.join(relative_img_dir, img2))

        html_file.write(tempStr)
        html_file.close()
        return
        
    # make single 
    def makeImageMapHTML(self, labtekId, resD, imageName, pheno, map_index=0):        
        
        tempStr = \
"""                
<map name="map_%i">
""" % map_index
        
        y_step = 20.0 
        x_step = 20.0 
        
        y_offset = 1
        x_offset = 22
        
        i=0
        for y in range(31, -1, -1): 
            for x in range(0, 12):
                
                i += 1                
                pos = '%03i' % i
                id = '%s--%s' % (labtekId, pos)                
                value = resD[id][pheno]
                img1 = resD[id]['img1']
                img2 = resD[id]['img2']
                pagename = os.path.join('.', labtekId, '%s.html' % id)
                #movieName = resD[id]['movieName']
                #moviePageName = os.path.join('.', labtekId, os.path.basename(movieName).replace('.avi', '.html'))
                gene = resD[id]['gene']
                sirna = resD[id]['sirna']
                
                infoD = {'number': resD[id]['number'], 
                         'ratio': resD[id]['ratio'], 
                         'fractionIndex': resD[id]['fractionIndex']}
                self.makeRawDataLinkPage(labtekId, img1, img2, gene, sirna, id, infoD)
                
                tempStr += \
"""<area shape="rect" coords="%i, %i, %i, %i" alt="" target="_blank" 
    title="%i/384, pos: %s -- %s -- %s -- value: %.2f" href="%s">
""" % ((x_offset + x * x_step), (y_offset + y * y_step), (x_offset + (x + 1) * x_step), (y_offset + (y + 1) * y_step), i, pos, gene, sirna, value, pagename)
 
        tempStr += \
"""
</map>
"""       

        tempStr += \
"""
<td align="center"><img src="%s" border="0" usemap="#map_%i"></td>
""" % (imageName, map_index)

        return tempStr
        

    def exportPlateViewHTML(self, filename, labtekId, entryD):
        
        file = open(os.path.join(self.html_dir, filename), 'w')
        
        # title and header
        tempStr="""
<html>        
  <head><title>Plate View: %s</title></head>
  <body>
    <h2>Basic Results for Plate %s </h2>
    <table align="left" valign="bottom" border="0" cellspacing="10" cellpadding="0">
""" % (labtekId, labtekId)

        file.write(tempStr)

        labtekRowString = ""

        # write headers for plates
        if entryD['Plate'].has_key('header'):
            header = entryD['Plate']['header']
        
            if not header is None:

                # ----- HEADER OF EACH ROW
                # header: name of labteks
                labtekRowString += """
    <tr>
"""
                for colTitle in header:
                    labtekRowString += """
    <th align="center"> %s </th>
""" % colTitle
                labtekRowString += """
    <th align="center"> Legends </th>
"""
                labtekRowString += """
    </tr>
    <tr valign="bottom">
"""
    
        # ----- plate density plots
        imageL = entryD['Plate']['images']
        if entryD['Plate'].has_key('height'):
            heightL = entryD['Plate']['height']
        else: 
            heightL = [None for x in imageL]
        for img_filename, height in zip(imageL, heightL):
            # no map: plot simply the image
            if not entryD['Plate'].has_key('map'):
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
            elif not entryD['Plate']['map'].has_key(img_filename):#img_filename != rowD[rowTitle]['map'][0]:

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
                labtekRowString += entryD['Plate']['map'][img_filename]
                    
        
        # --- legends
        imageL = entryD['Legend']['images']
        if entryD['Legend'].has_key('height'):
            heightL = entryD['Legend']['height']
        else: 
            heightL = [None for x in imageL]
        labtekRowString += """
    <td align="center">
"""
        for img_filename, height in zip(imageL, heightL):

            labtekRowString += """
    <img src="%s" border="0" height=%i> 
    <p>
""" % (img_filename, height)
        labtekRowString += """
    </td>
"""


#            labtekRowString += """
#    <td align="center">
#    <table align="left" border="0" cellspacing="10" cellpadding="0">
#    <tr align="center"><img src="%s" border="0" height=%i></tr>    
#    </table>
#    </td>
#""" % (img_filename, height)

        
        # --- end of row        
        labtekRowString += """
</tr>
"""
        
        # --- end of table
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
        
        
#        rowD['Plate'] = {'header': ['Cell Count', 'Ratio', 'Fraction Index'],
#                         'images': [os.path.join(relativePlotDir, '%s--%s.png') % (pheno, labtekId)
#                                    for pheno in phenoL],
#                         'map': mapD}
#        rowD['Legend'] = {'header': None,
#                          'images': [os.path.join(relativePlotDir, 'legend_%s--%s.png') % (pheno, labtekId)
#                                     for pheno in phenoL],
#                          'height': [70, 70, 70]}


        return
        
        
    def exportPlateHTML(self, filename, labtekId, rowD, rowL):
        file = open(os.path.join(self.html_dir, filename), 'w')
        
        # title and header
        titleString="""
<html>        
  <head><title>Plate View: %s</title></head>
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
                # HERE I AM
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

    def prepareData(self, resD, getLabel, expL, flip_axis=False):
        
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
        if flip_axis:
            hitData = r.matrix(hitValueL, nrow=self.nb_col, byrow=False)        
        else:
            hitData = r.matrix(hitValueL, nrow=self.nb_row, byrow=True)        
        #hitData = r("function(x,y){matrix(x, byrow=TRUE, nrow=y)[y:1,]}")(hitValueL, self.nb_row)

        return hitData    


    def plotArray(self, hitDataParam, filename, plotDir=None, labeledPosD=None,
                  title='', type='png',
                  label_colors = True, colVecL=None, legend_plot = False, legend_filename=None, 
                  control_color='black', numeric_as_int = False,
                  grid=True, max_val=0.3, min_val=0.0, cut_data = True, axis=True, flip_axis=False):
        if plotDir is None:
            plotDir = self.plotDir
        full_filename = os.path.join(plotDir, filename)
    
        #hitData = copy.deepcopy(hitDataParam)                
        hitData = hitDataParam
        #nrow = len(hitData)
        #ncol = len(hitData[0])
        
        if flip_axis:
            nrow = self.nb_col
            ncol = self.nb_row
        else:
            nrow = self.nb_row
            ncol = self.nb_col
        
        xL = range(1, ncol + 1)
        yL = range(1, nrow + 1)
        
        rdev = plot_utilities.RDevice(name = full_filename, title=title, plotType=type, 
                                      width=22+20*ncol, height=22+20*nrow)

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
                min_rel_cell_count = min_val
                for i in range(nrow):
                    for j in range(ncol):
                        hitData[i][j] = max(min_rel_cell_count, min(max_rel_cell_count, hitData[i][j])) 
            else:
                max_rel_cell_count = max([max(x) for x in hitData.tolist() ]) 
                min_rel_cell_count = min([min(x) for x in hitData.tolist() ])

            # we set the maximal/minimal values
            max_rel_cell_count = max_val
            min_rel_cell_count = min_val
            
            if numeric_as_int:
                nb_colors = max_rel_cell_count
            else:
                nb_colors = 500

            if colVecL is None:                    
                pattern = [(0,0,0),(0.7,0,0),(1,1,0),(1,1,1)]
                colVecL = colors.make_colors(pattern, nb_colors)
            colBreaksL = [1.0/ len(colVecL) * x * (max_rel_cell_count - min_rel_cell_count) + min_rel_cell_count  for x in range(len(colVecL) + 1)]

            if legend_plot:
                self.plotNumLegend(colVecL, colBreaksL, 11, filename=legend_filename, type=type, int_labels=numeric_as_int, legendDir = plotDir)

        axisSize = .75
        if axis:
            r("par(mar=c(1.6,1.6,0.1,0.1))")
        else:
            r("par(mar=c(0.1,0.1,0.1,0.1))")
                    
        r.image(xL, yL, r.t(hitData), axes = False, ann=False, cex=1, col=colVecL, breaks=colBreaksL)
        r.box()
        if not labeledPosD is None:
            for label in labeledPosD.keys():
                posL = labeledPosD[label]
                if len(posL) > 0:
                    xlL = [(int(x)-1) % self.nb_col + 1 for x in posL]
                    ylL = [(int(x)-1) / self.nb_col + 1 for x in posL]
                    r.points(xlL, ylL, pch=label, col=control_color, cex=axisSize)
        # grid
        if grid:
            for i in range(self.nb_row):
                r.abline(h=i+.5, lty=3, lwd=1, col='grey')
            for i in range(self.nb_col):
                r.abline(v=i+.5, lty=3, lwd=1, col='grey')

        if axis:
            r.axis(1, at=xL, labels=[str(x) for x in xL], tick=False, line=-1.0, cex_axis=axisSize)
            #r.axis(2, at=yL, labels=[str(y) for y in yL].reverse(), tick=False, line=-1.0, cex_axis=axisSize)
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
        tickBreakL = [min_break + float(x) / (nb_breaks - 1) * (max_break - min_break)  for x in range(nb_breaks)]
        #print filename, tickBreakL
        if int_labels:
            labels = ['%i' % int(x) for x in tickBreakL]
        else:
            labels = ['%.1f' % x for x in tickBreakL]
        
        rdev = plot_utilities.RDevice(name = full_filename, title='', plotType=type, 
            width=320, height=80)
        r("par(mar=c(4,1,0,1), las=2)")
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


from optparse import OptionParser
import pickle

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-p", "--pickle_file", dest="pickle_file",
                      help="pickle file containing result table")
    parser.add_option("-n", "--new_values", dest="new_val_file",
                      help="text file containing new values")
    
    
    (options, args) = parser.parse_args()

    if not (options.pickle_file is None):
        file = open(options.pickle_file, 'r')
        resD = pickle.load(file)
        file.close()

        # TEMPORARY: include new values.
        add_type_info(resD)
    
    if not (options.new_val_file is None):
        add_new_values(resD, options.new_val_file)

    hg = HTMLGenerator()
    hg(resD)

    
        
