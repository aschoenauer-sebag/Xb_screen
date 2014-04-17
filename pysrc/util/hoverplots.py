"""Example of how to use wx tooltips on a matplotlib figure window.
Adapted from http://osdir.com/ml/python.matplotlib.devel/2006-09/msg00048.html"""

import matplotlib as mpl
mpl.use('WXAgg')
mpl.interactive(False)
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from pylab import get_current_fig_manager as gcfm
import wx
import numpy as np
import random
from trajPack.exploiting_clustering import filmCharacterization, testing,\
    exploitingKMeans
    
'''
    To produce plots on which you can hover and get the name of the experiments. All need to be supplied is 
    the raw data, the number of clusters and the testing edges (0.05 usually), the function hoveringPlot calls
    the function that performs the clustering and statistical tests. 
'''


class onHover(object):
    def __init__(self, m, ctrlStatus,labls, colors, Td=False):
        self.figure = pl.figure()
        self.axis = self.figure.add_subplot(111) if not Td else self.figure.add_subplot(111, projection ='3d')
#ATTENTION NE MARCHE PAS EN TROIS D
        # create a long tooltip with newline to get around wx bug (in v2.6.3.3)
        # where newlines aren't recognized on subsequent self.tooltip.SetTip() calls
        self.tooltip = wx.ToolTip(tip='Pas de points ici\n')
        gcfm().canvas.SetToolTip(self.tooltip)
        self.tooltip.Enable(False)
        self.tooltip.SetDelay(0)
        self.figure.canvas.mpl_connect('motion_notify_event', self._onMotion)

        self.dataX = m[:,0]
        self.dataY = m[:,1]
        self.dataZ = m[:,2]
        self.text=labls
        if not Td: 
            self.axis.scatter(self.dataX, self.dataY, color=colors, label='myplot') 
        else:  
            self.axis.scatter(self.dataX, self.dataY,self.dataZ, color=colors, label='myplot')
            
    def _onMotion(self, event):
        collisionFound = False
        #ind=
        if event.xdata != None and event.ydata != None: # mouse is inside the axes
            for i in xrange(len(self.dataX)):
                radius = 0.15
                if abs(event.xdata - self.dataX[i]) < radius and abs(event.ydata - self.dataY[i]) < radius:
                    top = tip='%s,%s' % (self.text[i][0], self.text[i][1])
                    self.tooltip.SetTip(tip) 
                    self.tooltip.Enable(True)
                    collisionFound = True
                    break
        if not collisionFound:
            self.tooltip.Enable(False)
            
def hoveringPlot(k,seuilp,seuilc,data, length, who, ctrlStatus, genes, sirna,\
                  show=True, ctrl_only=False):
    labels, percentages=exploitingKMeans(data, length, ctrlStatus, k, who, genes, sirna, False)
#    percentages = filmCharacterization(labels, length, ctrlStatus, None, None, False, False, False)
    
    print "Using statistical tests for colors"
    ctest, cptest, ptest, pptest = testing(labels, ctrlStatus, length, who)
    if show:
        cperc=[]; ccols=[]; pperc=[]; pcols=[]
        for i,el in enumerate(who):
#            if el==
#            elif 
#            elif ctrlStatus[i]==0: ccols.append('black')
#            else: pcols.append('black')
            if ctrlStatus[i]==0:
                cperc.append(percentages[i])
#                if el==('LT0068_45--ex2005_07_15--sp2005_06_13--tt18--c2', '00074_01'):
#                    ccols.append('red')
#                elif el==('LT0068_45--ex2005_07_15--sp2005_06_13--tt18--c2', '00315_01'):
#                    ccols.append('green')
#                else: ccols.append("black") 
                if ctest[el[0]]<seuilc and cptest[el[0]]<seuilc: ccols.append('pink')
                else: ccols.append("red")
            else:
                pperc.append(percentages[i])
                if ptest[el]<seuilp and pptest[el]<seuilp:pcols.append('black')
                else:pcols.append('green')
                
        ww=zip(who, genes)
        if ctrl_only:
            localc=list(ccols)
            localp=np.array(cperc)
            localw=list(np.array(ww)[np.where(np.array(ctrlStatus)==0)])
            ex=onHover(localp, None, localw, localc)
            pl.show()
        else:
            for k in range(11):
                localc=list(ccols)
                localc.extend(pcols[1000*k:1000*(k+1)])
                localp=np.vstack((cperc, pperc[1000*k:1000*(k+1)]))
                localw=list(np.array(ww)[np.where(np.array(ctrlStatus)==0)])
                localw.extend(np.array(ww)[np.where(np.array(ctrlStatus)==1)][1000*k:1000*(k+1)])
                ex=onHover(localp, None, localw, localc)
                pl.show()
        
    return percentages, ctest, cptest, ptest




#example = wxToolTipExample()
#pl.show()
