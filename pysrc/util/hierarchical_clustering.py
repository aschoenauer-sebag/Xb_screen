### hierarchical_clustering.py
#Copyright 2005-2012 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#################
### Produces and hierarchically clustered heatmap
#################
import getpass
if getpass.getuser()=='aschoenauer':
    import matplotlib as mpl
    mpl.use('Agg')
    show = False
elif getpass.getuser()=='lalil0u':
    import matplotlib as mpl
    show=True

import matplotlib.pyplot as pylab
import scipy
from optparse import OptionParser
import cPickle as pickle
from sklearn.cluster import MiniBatchKMeans
import scipy.cluster.hierarchy as sch
import fastcluster
import scipy.spatial.distance as dist
import numpy
import numpy as np
import string
import time
import sys, os, pdb
import getopt
from collections import Counter

from util.listFileManagement import EnsemblEntrezTrad, multipleGeneListsToFile
from tracking.trajPack import featuresNumeriques, featuresSaved

def replotHeatmap(folder, data_filename, indices, outputfile,action='hierarchical', level=0.4,trad=False,
                  num_clusters=8, labels_filename='labelsKM_k{}_halfM.pkl', pcaed_filename='halfM_max_05_pcaed.pkl'):
    
    #DONC FAIRE TOUT DANS LA MEME FONCTION
    
    features=list(featuresNumeriques); features.append('mean persistence'); features.append('mean straight')
    f=open(os.path.join(folder,data_filename))
    r,  who,ctrlStatus, length, genes, sirna, time_length=pickle.load(f); f.close()
    
    r=np.hstack((r[:,:len(featuresNumeriques)], r[:,featuresSaved.index('mean persistence'), np.newaxis], r[:, featuresSaved.index('mean straight'), np.newaxis]))
    
    r=np.hstack((r, np.array(time_length)[:,np.newaxis])); fL=list(features);fL.append('time length')
    
    mean_=np.mean(r,0)
    for k,feature in enumerate(features):
        print feature, mean_[k]
    nr=(r-mean_)/np.std(r,0)

    #nr=np.hstack((nr, indices[:,np.newaxis])); fL.append('hierarchical cluster')
    small_nr=None
    
    if action=='hierarchical':
        num_clusters=len(np.bincount(indices)); begin_=1
        print 'Going for hierarchical clustering (ward, euclidean) with {} clusters'.format(num_clusters-1)

    elif action=='kmeans':
        print 'Going for mini batch k-means with {} clusters'.format(num_clusters); begin_=0
        try:
            f=open(os.path.join(folder, labels_filename.format(num_clusters)), 'r')
            labels,percentages, who, length=pickle.load(f); f.close()
        except OSError:
            print 'File Error ', os.path.join(folder, labels_filename.format(num_clusters))
        else:
            indices=labels

#        f=open(os.path.join(folder, pcaed_filename), 'r')
#        _,pcaed_data=pickle.load(f); f.close()
#        
#        model = MiniBatchKMeans(n_clusters=num_clusters, batch_size = 2000, init='k-means++',n_init=1000,max_iter=1000, max_no_improvement=100, compute_labels = True)
#        indices = model.fit(pcaed_data[:,:7])
#        indices=indices.labels_

    for k in range(begin_, num_clusters):
        where_=np.where(np.array(indices)==k)[0]
        np.random.shuffle(where_)
        small_nr = np.vstack((small_nr, nr[where_[:1000]])) if small_nr is not None else nr[where_[:1000]]
    print small_nr.shape

    print 'Showing trajectory clusters'
    heatmap(small_nr.T,fL, range(small_nr.shape[0]), 'ward', None, 'euclidean', None, 
            color_gradient='red_white_blue', filename=outputfile+'TRAJ', trad=False, save=False)
    if action=='kmeans':
        heatmap(percentages[:1723], genes[:1723],range(begin_, num_clusters), None, 'ward', None, 'euclidean', 
            color_gradient='red_white_blue', filename=outputfile+'MOV', trad=False, save=False, level=level)
        heatmap(percentages[1723:], genes[1723:],range(begin_, num_clusters), None, 'ward', None, 'euclidean', 
            color_gradient='red_white_blue', filename=outputfile+'CTRL', trad=False, save=False, level=level)

    
    return

################# Perform the hierarchical clustering #################

def heatmap(x, row_header, column_header, row_method,
            column_method, row_metric, column_metric,
            color_gradient, filename, log=False, trad=False, level=0.4, save=True):
    
    print "\nPerforming hiearchical clustering using %s for columns and %s for rows" % (column_metric,row_metric),
    if numpy.any(numpy.isnan(x)):
        sys.stderr.write("WARNING, there are NaN values in the data. Hence distances with data elements that have NaN values will have value NaN, which might perturb the hierarchical clustering.")
        
    """
    This below code is based in large part on the protype methods:
    http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
    http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre
    
    Possibilities for methods: single, complete, average, centroid, median, ward
    
    Possibilities for metrics: 'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 
    'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching', 
    'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule

    x is an m by n ndarray, m observations, n genes, or m rows,n columns
    
    WARNING WARNING
    This is a modified version to work with "big data" (starting with m=50,000). Indeed, the previous version actually stores
    the distance matrix in the memory which makes it crash. Here, we use the package fastcluster (see http://danifold.net/fastcluster.html)
    in its memory-efficient implementation.
    The parameter method must be one of 'single', 'centroid', 'median', 'ward'
    
    """
    print level
        
    #for export
    if numpy.any(~numpy.array([type(s)==str for s in row_header])):
        row_header=[str(el) for el in row_header]
    if numpy.any(~numpy.array([type(s)==str for s in column_header])):
        column_header=[str(el) for el in column_header]
        
    ### Define the color gradient to use based on the provided name
    n = len(x[0]); m = len(x)
    if color_gradient == 'red_white_blue':
        cmap=pylab.cm.bwr
    if color_gradient == 'red_black_sky':
        cmap=RedBlackSkyBlue()
    if color_gradient == 'red_black_blue':
        cmap=RedBlackBlue()
    if color_gradient == 'red_black_green':
        cmap=RedBlackGreen()
    if color_gradient == 'yellow_black_blue':
        cmap=YellowBlackBlue()
    if color_gradient == 'seismic':
        cmap=pylab.cm.seismic
    if color_gradient == 'green_white_purple':
        cmap=pylab.cm.PiYG_r
    if color_gradient == 'coolwarm':
        cmap=pylab.cm.coolwarm

    ### Scale the max and min colors so that 0 is white/black
    vmin=numpy.nanmin(x)
    vmax=numpy.nanmax(x)
    vmax = max([vmax,abs(vmin)])
    #vmin = vmax*-1
#    if log:
#        norm = mpl.colors.LogNorm(vmin, vmax) ### adjust the max and min to scale these colors
#    elif normalization:
#        norm = mpl.colors.Normalize(10**(-70), 1)
#    else:
    if numpy.any(x<0):
        norm = mpl.colors.Normalize(-2,2)
    else:
        norm = mpl.colors.Normalize(0,50)
    ### Scale the Matplotlib window size
    default_window_hight = 8.5
    default_window_width = 12
    fig = pylab.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    color_bar_w = 0.015 ### Sufficient size to show
        
    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    #if row_method != None: w1 = 
    [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.22,0.2,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
    width_between_ax1_axr = 0.004
    height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix
    
    # axr, placement of row side colorbar
    [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
    axr_x = ax1_x + ax1_w + width_between_ax1_axr
    axr_y = ax1_y; axr_h = ax1_h
    width_between_axr_axm = 0.004

    # axc, placement of column side colorbar
    [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
    axc_x = axr_x + axr_w + width_between_axr_axm
    axc_y = ax1_y + ax1_h + height_between_ax1_axc
    height_between_axc_ax2 = 0.004

    # axm, placement of heatmap for the data matrix
    [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
    axm_x = axr_x + axr_w + width_between_axr_axm
    axm_y = ax1_y; axm_h = ax1_h
    axm_w = axc_w

    # ax2, placement of dendrogram 2, on the top of the heatmap
    [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls hight of the dendrogram
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
    ax2_w = axc_w

    # axcb - placement of the color legend
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.09]

    # Compute and plot top dendrogram
    if column_method != None:
        start_time = time.time()
#        d2 = dist.pdist(x.T)
#        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
        
        Y2 = fastcluster.linkage_vector(x.T, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z2 = sch.dendrogram(Y2)
        ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
        ax2.set_xticks([]) ### Hides ticks
        ax2.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print 'Column clustering completed in %s seconds' % time_diff
    else:
        ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
        
    # Compute and plot left dendrogram.
    if row_method != None:
        start_time = time.time()
#        d1 = dist.pdist(x)
#        D1 = dist.squareform(d1)  # full matrix
        ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True) # frame_on may be False
        
        Y1 = fastcluster.linkage_vector(x, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z1 = sch.dendrogram(Y1, orientation='right')
        ind1 = sch.fcluster(Y1,level*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
        ax1.set_xticks([]) ### Hides ticks
        ax1.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print 'Row clustering completed in %s seconds' % time_diff
    else:
        ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
    if save:
        print 'Saving flat clusters in', 'Flat_clusters_{}_{}.pkl'.format(filename, level) 
        f=open('Flat_clusters_{}_{}.pkl'.format(filename,level), 'w')
        pickle.dump([ind1, ind2],f); f.close()
    
    if trad:
        if len(row_header)>100:
            genes=list(row_header)
            clustering = numpy.array(ind1)
        elif len(column_header)>100:
            genes=list(column_header)
            clustering=numpy.array(ind2)
        else:
            print 'Tell which of column and row is the gene list'
            pdb.set_trace()
        #il faut d'abord traduire de SYMBOL en ENSEMBL
        trad = EnsemblEntrezTrad('../data/mapping_2014/mitocheck_siRNAs_target_genes_Ens75.txt')
        trad['ctrl']='None'
        
        result=[Counter([trad[genes[k]] for k in numpy.where(clustering==cluster)[0]]).keys() for cluster in range(1,numpy.max(clustering)+1)]
        for geneList in result:
            for i,gene in enumerate(geneList):
                if '/' in gene:
                    geneList[i]=gene.split('/')[0]
                    geneList.append(gene.split('/')[1])
                
        #ensuite on va enregistrer les genes des differents clusters dans differents fichiers
        #background par defaut c'est genes_list.txt
        print "Nb of cluster found", numpy.max(clustering)
        multipleGeneListsToFile(result, ['Cluster {}'.format(k+1) for k in range(numpy.max(clustering))], 'gene_cluster_{}_{}.txt'.format(column_method, filename))
    
    # Plot distance matrix.
    axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
    xt = x
    if column_method != None:
        idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
        xt = xt[:,idx2]
        ind2 = ind2[:,idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
    if row_method != None:
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        xt = xt[idx1,:]   # xt is transformed x
        ind1 = ind1[idx1,:] ### reorder the flat cluster to match the order of the leaves the dendrogram
    ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
    im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
    axm.set_xticks([]) ### Hides x-ticks
    axm.set_yticks([])

    # Add text
    new_row_header=[]
    new_column_header=[]
    for i in range(x.shape[0]):
        if row_method != None:
            if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i, '  {}'.format(row_header[idx1[i]]))
            new_row_header.append(row_header[idx1[i]])
        else:
            if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i, ' {}'.format(row_header[i])) ### When not clustering rows
            new_row_header.append(row_header[i])
    for i in range(x.shape[1]):
        if column_method != None:
            if len(column_header)<100:
                axm.text(i, -0.9, '{}'.format(column_header[idx2[i]]), rotation=270, verticalalignment="top") # rotation could also be degrees
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            if len(column_header)<100:
                axm.text(i, -0.9, '{}'.format(column_header[i]), rotation=270, verticalalignment="top")
            new_column_header.append(column_header[i])

    # Plot colside colors
    # axc --> axes for column side colorbar
    if column_method != None:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        dc = numpy.array(ind2, dtype=int)
        dc.shape = (1,len(ind2)) 
        im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
        axc.set_xticks([]) ### Hides ticks
        axc.set_yticks([])
    
    # Plot rowside colors
    # axr --> axes for row side colorbar
    if row_method != None:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        dr = numpy.array(ind1, dtype=int)
        dr.shape = (len(ind1),1)
        #print ind1, len(ind1)
        cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
        axr.set_xticks([]) ### Hides ticks
        axr.set_yticks([])

    # Plot color legend
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap,norm=norm, orientation='horizontal')
    axcb.set_title("colorkey")
    
    filename = 'Clustering-%s-hierarchical_%s_%s.pdf' % (filename[:10],column_metric,column_method)
    cb.set_label("Range")
    exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2)

    ### Render the graphic
    if len(row_header)>50 or len(column_header)>50:
        pylab.rcParams['font.size'] = 5
    else:
        pylab.rcParams['font.size'] = 8

    pylab.savefig(filename)
    print 'Exporting:',filename
    
    if show:
        pylab.show()
    if trad:
        return result

def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax=x.max()
    vmin=x.min()
    if vmax<0 and vmin<0: direction = 'negative'
    elif vmax>0 and vmin>0: direction = 'positive'
    else: direction = 'both'
    if direction == 'both':
        vmax = max([vmax,abs(vmin)])
        vmin = -1*vmax
        return vmax,vmin
    else:
        return vmax,vmin
    
################# Export the flat cluster data #################

def exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """
    
    filename = string.replace(filename,'.pdf','.txt')
    export_text = open(filename,'w')
    column_header = string.join(['UID','row_clusters-flat']+new_column_header,' ')+'\n' ### format column-names for export
    export_text.write(column_header)
    column_clusters = string.join(['column_clusters-flat','']+ map(str, ind2),' ')+'\n' ### format column-flat-clusters for export
    export_text.write(column_clusters)
    
    ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
    new_row_header = new_row_header[::-1]
    xt = xt[::-1]
    
    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_text.write(string.join([new_row_header[i],str(ind1[i])]+map(str, row),' ')+'\n')
        i+=1
    export_text.close()
    
    ### Export as CDT file
    filename = string.replace(filename,'.txt','.cdt')
    export_cdt = open(filename,'w')
    column_header = string.join(['UNIQID','NAME','GWEIGHT']+new_column_header,'\t')+'\n' ### format column-names for export
    export_cdt.write(column_header)
    eweight = string.join(['EWEIGHT','','']+ ['1']*len(new_column_header),'\t')+'\n' ### format column-flat-clusters for export
    export_cdt.write(eweight)
    
    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_cdt.write(string.join([new_row_header[i]]*2+['1']+map(str, row),'\t')+'\n')
        i+=1
    export_cdt.close()

################# Create Custom Color Gradients #################
#http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html

def RedBlackSkyBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.9),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackGreen():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'blue': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'green':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def YellowBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.8),
                       (0.5, 0.1, 0.0),
                       (1.0, 1.0, 1.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
    ### modulate between blue and cyan using the last y var in the first green tuple
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

################# General data import methods #################

def importData(filename):
    start_time = time.time()
    matrix=[]
    row_header=[]
    first_row=True

    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]
        
    for line in open(filename,'rU').xreadlines():         
        t = string.split(line[:-1],'\t') ### remove end-of-line character - file is tab-delimited
        if first_row:
            column_header = t[1:]
            first_row=False
        else:
            if ' ' not in t and '' not in t: ### Occurs for rows with missing data
                s = map(float,t[1:])
                if (abs(max(s)-min(s)))>0:
                    matrix.append(s)
                    row_header.append(t[0])
            
    time_diff = str(round(time.time()-start_time,1))
    try:
        print '\n%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
    except Exception:
        print 'No data in input file.'; force_error
    return numpy.array(matrix), column_header, row_header
  
if __name__ == '__main__':
    
    parser = OptionParser(usage="usage: %prog [options]")    
    parser.add_option('--level', type=float, default=0.4)
    parser.add_option('--outputname', type=str, default='halfM_median_05')
    
    (options, args) = parser.parse_args()
    
    f=open('../resultData/features_on_films/{}_pcaed.pkl'.format(options.outputname))
    _, narr=pickle.load(f)
    f.close()
    narr=narr[:,:7]
    print 'data loaded'
    r=heatmap(narr.T, range(narr.shape[1]),range(narr.shape[0]), 
              'ward', 'ward', 'euclidean', 'euclidean', 'red_black_sky', 'traj_HCsiRNA_iter5_clustering',trad=False,
              level=options.level)

    
    
#    ################  Default Methods ################
#    row_method = 'average'
#    column_method = 'single'
#    row_metric = 'cityblock' #cosine
#    column_metric = 'euclidean'
#    color_gradient = 'red_white_blue'
#    
#    """ Running with cosine or other distance metrics can often produce negative Z scores
#        during clustering, so adjustments to the clustering may be required.
#        
#    see: http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
#    see: http://docs.scipy.org/doc/scipy/reference/spatial.distance.htm  
#    color_gradient = red_white_blue|red_black_sky|red_black_blue|red_black_green|yellow_black_blue|green_white_purple'
#    """
#    ################  Comand-line arguments ################
#    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
#        print "Warning! Please designate a tab-delimited input expression file in the command-line"
#        print "Example: python hierarchical_clustering.py --i /Users/me/logfolds.txt"
#        sys.exit()
#    else:
#        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','row_header','column_method',
#                                                    'row_metric','column_metric','color_gradient'])
#        for opt, arg in options:
#            if opt == '--i': filename=arg
#            elif opt == '--row_header': row_header=arg
#            elif opt == '--column_method': column_method=arg
#            elif opt == '--row_metric': row_metric=arg
#            elif opt == '--column_metric': column_metric=arg
#            elif opt == '--color_gradient': color_gradient=arg
#            else:
#                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
#            
#    matrix, column_header, row_header = importData(filename)
#
#    if len(matrix)>0:
#        try:
#            heatmap(matrix, row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, filename)
#        except Exception:
#            print 'Error using %s ... trying euclidean instead' % row_metric
#            row_metric = 'euclidean'
#            try:
#                heatmap(matrix, row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, filename)
#            except IOError:
#                print 'Error with clustering encountered'
