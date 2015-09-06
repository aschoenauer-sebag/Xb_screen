import numpy as np
import matplotlib.pyplot as p
import matplotlib as mpl
import os
from operator import itemgetter
import cPickle as pickle
from phenotypes.drug_screen_utils import CLASSES

basic_colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ]
markers =('o', 'v','*', '^', '<','8', '>',  's', 'p',  'h', 'H', 'D', 'd', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
linestyles=['-', '--', '-.', ':']
couleurs=[]
couleurs.append("#9e0142")
couleurs.append("#5e4fa2")

couleurs.append("#d53e4f")
couleurs.append("#3288bd")

couleurs.append("#f46d43")
couleurs.append("#66c2a5")

couleurs.append("#fdae61")
couleurs.append("#abdda4")

couleurs.append("#fee08b")
couleurs.append("#e6f598")

couleurs.append('#a50026')
couleurs.append('#d73027')

couleurs.append("#9e0142")
couleurs.append("#5e4fa2")
couleurs.append("#9e0142")
couleurs.append("#5e4fa2")

def barplot_genelist(gg):
    
    res_=[(el, gg[el]/2452.0) for el in np.array(gg.keys())[np.argsort(gg.values())]]
    f=p.figure()
    ax=f.add_subplot(111)
    ax.grid(True)
    ax.bar(range(len(res_)), [el[1] for el in res_], alpha=0.7, color='blue')
    ax.set_xticks(np.array(range(len(res_)))+0.5)
    ax.set_xticklabels([el[0] for el in res_], rotation='vertical')
    
    ax.set_xlabel('Hit category')
    ax.set_ylabel('Percentage of gene in each category')
    ax.set_title('Distribution of genes selected for target inference')
    
    p.show()

def plotTransportCost(cost):
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.93,0.18,0.02]
    f=p.figure()
    ax=f.add_subplot(111)
    ax.matshow(cost, cmap=mpl.cm.YlOrRd)
    ax.set_xticks(range(len(CLASSES)))
    ax.set_xticklabels(CLASSES, rotation='vertical')
    
    ax.set_yticks(range(len(CLASSES)))
    ax.set_yticklabels(CLASSES)
    axcb = f.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False) 
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=mpl.cm.YlOrRd,orientation='horizontal', ticks=[])
    #cb.ax.set_xticklabels(range(1,7))
    p.show()



def plotBar(counts, second_counts, legends=['Very low (Uncertain)', 'Low (Uncertain)', 'Medium (Supportive)', 'High (Supportive)']):
    values=[counts[0][el] for el in legends]
    f=p.figure()
    ax=f.add_subplot(121)
    ax.bar(range(len(values)), values, alpha=0.7, color='blue', label="Scilifelab")
    for i,el in enumerate(legends):
        ax.text(i, counts[0][el]+2, el)
    ax.set_xticklabels([])
    
    if second_counts !=None:
        values=[second_counts[0][el] for el in legends]
    ax.bar(range(len(values)), values, alpha=0.7, color='green', label="Intersection w/ Mitocheck as found by MotIW")
    ax.legend()
        
    width=0.8
    values=sorted(counts[1].iteritems(), key=itemgetter(1))
    
    ax=f.add_subplot(122)
    ax.bar(range(len(counts[1])), [el[1] for el in values], width, color='blue', alpha=0.7)
    ax.set_xticklabels([])
    for i,el in enumerate(values):
        if el[1]>20:
            ax.text(i, el[1]-2, el[0], rotation='vertical')
        else:
            ax.text(i, len(el[0])-2, el[0], rotation='vertical')
    if second_counts is not None:
        values=[second_counts[1][el[0]] for el in values]
        ax.bar(range(len(values)), values, width, color='green', alpha=0.7)
    
    p.show()


def plotLabelDistributions(distributions, nb_experiments, nb_dist=8):
    f,axes=p.subplots(1,nb_dist, sharey=True)
    blist=[]
    for k in range(nb_dist):
        n,bins,patches=axes[k].hist(distributions[:nb_experiments,k], bins=50, color='#B20000', normed=True,alpha=1, label='Exp')
        blist.append(bins); axes[k].set_ylim(0,0.4); axes[k].grid(True)
    for k in range(nb_dist):
        n,bins,patches=axes[k].hist(distributions[nb_experiments:,k], bins=blist[k], color='#13C528', normed=True,alpha=0.6,label='Ctrl')
    axes[0].legend()
    p.show()

def plotBarPlot(means, ylabel, xlabels, xticklabels, title,name,target,  stds=None,folder=None, fig=None,ax=None, sh=True, save=True):    
    '''
    Xticklabels = sequence of labels for x ticks
    '''
    rects=[]
    ind = np.arange(means.shape[1])  # the x locations for the groups[0, 1, 2, 3, 4]
    width = 0.4 # the width of the bars
    space = 4.6#approx width*len(means)
    ind = np.array([space* el for el in ind])
    if ax is None or fig is None:
        fig = p.figure(figsize=(14,8))
        ax = fig.add_subplot(111)

    for k in range(len(means)):
        if stds is None or (stds is not None and stds[k] is None):
            rects.append(ax.bar(ind+width*k, means[k], width, color=couleurs[k]))
        elif stds[k] is not None:
            rects.append(ax.bar(ind+width*k, means[k], width, color=couleurs[k], yerr=stds[k]))
            
    # add some labels
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks([2+space*el for el in np.arange(means.shape[1])])
    ax.set_xticklabels( xticklabels )
    ax.set_ylim(0,6)
    ax.tick_params(axis='both', which='major', labelsize=5)
    fig.legend( rects, xlabels)
    p.grid(True)
    p.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    if sh:
        p.show()
    elif save:
        p.savefig(os.path.join(folder, '{}_barplot{}.png'.format(target, name)))
        
    return

def makeColorRamp(N, hex_output=False):
    basic_colors = couleurs[:6]
    if N<1:
        return []
    if N==1:
        return [basic_colors[0]]

    xvals = np.linspace(0, len(basic_colors)-1, N)

    single_channel = {}
    for c in range(3):
        xp = range(len(basic_colors))
        yp = [int(x[(1 + 2*c):(3+2*c)], base=16) for x in basic_colors]

        single_channel[c] = [x / 256.0 for x in np.interp(xvals, xp, yp)]

    if hex_output:
#            colvec = ['#' + hex(numpy.int32(min(16**4 * single_channel[0][i], 16**6 - 1) +
#                                            min(16**2 * single_channel[1][i], 16**4 - 1) +
#                                            min(single_channel[2][i], 16**2 -1) )).upper()[2:]
#                      for i in range(N)]
        colvec = ['#' + hex(np.int32(
                                        (((256 * single_channel[0][i]  ) +
                                          single_channel[1][i]) * 256 +
                                          single_channel[2][i]) * 256
                                          ))[2:]
                  for i in range(N)]

    else:
        colvec = zip(single_channel[0], single_channel[1], single_channel[2])

    return colvec

class AnnoteFinder:
    """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.
    
  Register this function like this:
    
  scatter(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

    def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
        self.data = zip(xdata, ydata, annotes)
        if xtol is None:
            xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
        if ytol is None:
            ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
        self.xtol = xtol
        self.ytol = ytol
        if axis is None:
            self.axis = pylab.gca()
        else:
            self.axis= axis
        self.drawnAnnotations = {}
        self.links = []

    def distance(self, x1, x2, y1, y2):
        """
        return the distance between two points
        """
        return math.hypot(x1 - x2, y1 - y2)

    def __call__(self, event):
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            if self.axis is None or self.axis==event.inaxes:
                annotes = []
                for x,y,a in self.data:
                    if  clickX-self.xtol < x < clickX+self.xtol and  clickY-self.ytol < y < clickY+self.ytol :
                        annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
                if annotes:
                    annotes.sort()
                    distance, x, y, annote = annotes[0]
                    self.drawAnnote(event.inaxes, x, y, annote)
                    for l in self.links:
                        l.drawSpecificAnnote(annote)

    def drawAnnote(self, axis, x, y, annote):
        """
        Draw the annotation on the plot
        """
        if (x,y) in self.drawnAnnotations:
            markers = self.drawnAnnotations[(x,y)]
            for m in markers:
                m.set_visible(not m.get_visible())
            self.axis.figure.canvas.draw()
        else:
            t = axis.text(x,y, "(%3.2f, %3.2f) - %s"%(x,y,annote), )
            m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
            self.drawnAnnotations[(x,y)] =(t,m)
            self.axis.figure.canvas.draw()

    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
        for x,y,a in annotesToDraw:
            self.drawAnnote(self.axis, x, y, a)