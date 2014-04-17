from tracking.trajPack import featuresHisto

DIVERGENCES=['transportation','etransportation', 'hellinger',  'total_variation']#'kullback_leibler',
WEIGHTS = [[1,1,1,1,1,1], [0,1,1,1,1,1], [1,0, 0, 0,0,0],
                        [0,1, 0, 0,0,0], [0,0, 1, 0,0,0],[0,0, 0, 1,0,0],[0,0, 0, 0,1,0],[0,0, 0, 0,0,1]]
ddWEIGHTS = [[1,1], [0,1], [1,0]]
WEIGHTNAMES = ['all', 'histograms only', 'numeric only']; WEIGHTNAMES.extend(featuresHisto)