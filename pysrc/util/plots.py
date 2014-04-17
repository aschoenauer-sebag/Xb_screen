import numpy as np

basic_colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ]

def makeColorRamp(N, hex_output=False):

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