import numpy as np
import matplotlib.pyplot as p
import os

basic_colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ]
markers =('o', 'v','*', '^', '<','8', '>',  's', 'p',  'h', 'H', 'D', 'd', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
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
#f46d43
#fdae61
#fee08b
#ffffbf
#d9ef8b
#a6d96a
#66bd63
#1a9850
#006837

#couleurs.extend(basic_colors)


#couleurs.append("#A50026")
#couleurs.append("#D73027")
#couleurs.append("#4575B4")
#couleurs.append("#ABDDA4")
#couleurs.append("#F46D43")
#couleurs.append("#FDAE61")
#couleurs.append("#A50026")
#couleurs.append("#D73027")

def plotBarPlot(means, ylabel, xlabels, xticklabels, title,name,target,  stds=None,folder=None, sh=True):    
    '''
    Xticklabels = sequence of labels for x ticks
    '''
    rects=[]
    ind = np.arange(means.shape[1])  # the x locations for the groups[0, 1, 2, 3, 4]
    width = 0.4 # the width of the bars
    space = 4.6#approx width*len(means)
    ind = np.array([space* el for el in ind])
    
    fig = p.figure(figsize=(14,8))
    ax = fig.add_subplot(111)
    for k in range(len(means)):
        if stds is None or (stds is not None and stds[k] is None):
            rects.append(ax.bar(ind+width*k, means[k], width, color=couleurs[k]))
        elif stds[k] is not None:
            rects.append(ax.bar(ind+width*k, means[k], width, color=couleurs[k], yerr=stds[k][:,0]))
            
    # add some labels
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks([2+space*el for el in np.arange(means.shape[1])])
    ax.set_xticklabels( xticklabels )
    ax.set_ylim(0,6)
    fig.legend( rects, xlabels)
    p.grid(True)
    p.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    if sh:
        p.show()
    else:
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