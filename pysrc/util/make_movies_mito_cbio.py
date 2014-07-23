import os
import cPickle as pickle
import shutil, pickle
from optparse import OptionParser

import vigra
import numpy as np

import pdb

BASE_DIR = '/share/data20T/mitocheck/compressed_data'
TRAJECTORY_DIR = '/share/data20T/mitocheck/tracking_results'

FEATURES = [            
            'ball number1','ball number2',
            'norm convex hull area',
            'corrected straightness index',
            'diffusion coefficient',
            'effective space length', 'effective speed',
            'entropy1', 'entropy2',#rayons 0, 14 ci-dessus
            'largest move',
            'mean squared displacement',
            'movement type',
            'signed turning angle', 
            ]

RADIUS = 11

f=open('/cbio/donnees/aschoenauer/workspace2/Xb_screen/resultData/features_on_films/distExp_123etctrl_minmax_50.pkl', 'r')
minMax = pickle.load(f); f.close()

FEATURE_RANGE = dict(zip(FEATURES, minMax))

#FEATURE_RANGE =  {'entropy1': (0.15, 0.75), 
#                  'mean squared displacement': (0.0, 50.0), 
#                  'entropy2': (0.0, 0.4), 
#                  'corrected straightness index': (0.0, 4.0), 
#                  'effective space length': (0.0, 60.0), 
#                  'norm convex hull area': (0.0, 30.0), 
#                  'movement type': (0.0, 0.8), 
#                  'ball number2': (0.1, 0.6), 
#                  'signed turning angle': (-3.141592653589793, 3.141592653589793), 
#                  'ball number1': (0.09, 2.5), 
#                  'largest move': (0.0, 15.0), 
#                  'effective speed': (0.0, 9.0), 
#                  'diffusion coefficient': (0.0, 10.0)
#                 }


# colors 
class ColorMap(object):
    def __init__(self):
        # divergent color maps
        self.div_basic_colors_intense = ["#E41A1C",
                                         "#377EB8",
                                         "#4DAF4A",
                                         "#984EA3",
                                         "#FF7F00" ]
        self.div_basic_colors_soft = ["#7FC97F",
                                      "#BEAED4",
                                      "#FDC086",
                                      "#FFFF99",
                                      "#386CB0" ]


    def getRGBValues(self, hexvec):
        single_channel = {}
        for c in range(3):
            single_channel[c] = [int(x[(1 + 2*c):(3+2*c)], base=16) / 256.0 for x in hexvec]
        rgbvals = zip(single_channel[0], single_channel[1], single_channel[2])
        return rgbvals

    def makeDivergentColorRamp(self, N, intense=True, hex_output=False):
        if intense:
            basic_colors = self.div_basic_colors_intense
        if not intense:
            basic_colors = self.div_basic_colors_soft

        cr = self.makeColorRamp(N, basic_colors, hex_output)
        return cr

    def makeColorRamp(self, N, basic_colors, hex_output=False):

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
#            colvec = ['#' + hex(np.int32(min(16**4 * single_channel[0][i], 16**6 - 1) +
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

    # cr = cm.makeColorRamp(256, ["#FF0000", "#FFFF00", "#FFFFFF"])
    def getColorFromMap(self, value, color_map, minval, maxval):
        if value > maxval:
            value = maxval
        if value < minval:
            value = minval
            
        N = len(color_map)
        value_range = np.arange(minval, maxval, float(maxval-minval)/N)
        index = sum(value > value_range) - 1
        if index < 0:
            index = 0
        if index > N-1:
            index = N-1
        
        return color_map[index]
    

class MovieMaker(object):
    #def __init__(self):
    
    def __init__(self, base_dir=BASE_DIR):
        self.base_dir=base_dir
        return
    
    def get_feature_ranges(self, base_dir):
        feat_max = []
        feat_min = []
        feat_mean = []
        feat_std = []
        for lt in os.listdir(base_dir):
            lt_dir = os.path.join(base_dir, lt)
            feature_files = filter(lambda x: x[0:6]=='Thomas', 
                                   os.listdir(lt_dir))
            
            for ff in feature_files:
                fp = open(os.path.join(lt_dir, ff), 'r')
                feat = pickle.load(fp)
                fp.close()
                feat_max.append(np.max(feat[0], axis=0))
                feat_min.append(np.min(feat[0], axis=0))
                feat_mean.append(np.mean(feat[0], axis=0))
                feat_std.append(np.std(feat[0], axis=0))
                
        max_mat = np.array(feat_max)
        min_mat = np.array(feat_min)
        mean_mat = np.array(feat_mean)
        std_mat = np.array(feat_std)
        
        print max_mat.shape
        print min_mat.shape
        print np.max(max_mat, axis=0)
        print np.min(min_mat, axis=0)
        
        maxvec = np.max(max_mat, axis=0).tolist()
        minvec = np.min(min_mat, axis=0).tolist()
        meanvec = np.mean(mean_mat, axis=0).tolist()
        stdvec = np.mean(std_mat, axis=0).tolist()
        
        res = dict(zip(FEATURES, zip(minvec, maxvec, meanvec, stdvec)))
        
        print res
                
        return
    
    def get_img_dir(self, id, base_dir=None):
        if base_dir is None:
            base_dir = self.base_dir
        
        lt, pos = id.split('--')
        labteks_found = filter(lambda x: x[:len(lt)]==lt, os.listdir(base_dir))
        if len(labteks_found) == 0:
            raise ValueError('%s not found in %s' % (lt, base_dir))
        if len(labteks_found) > 1:
            raise ValueError('multiple labteks found')
        labtek_dir = os.path.join(base_dir, labteks_found[0])
        pos_found = filter(lambda x: x[:len(pos)]==pos, os.listdir(labtek_dir))
        if len(pos_found) == 0:
            raise ValueError('%s not found in %s' % (pos, labtek_dir))
        if len(pos_found) > 1:
            raise ValueError('multiple positions found')
        img_dir = os.path.join(labtek_dir, pos_found[0])
        
        return img_dir
    
    def make_circle_offset(self, r):
        xvec = []
        yvec = []
        for x in range(-r, r):
            for y in range(-r, r):
                if x**2 + y**2 <= r:
                    xvec.append(x)
                    yvec.append(y)
        res = [np.array(xvec), np.array(yvec)]
        return res
    
    def make_movie(self, id, out_path, tempDir=None, sirna=None, gene=None, outDir=None,
                   trajectory_dir=None, feature_movie_dir=None, feature=None):
        
        in_path = self.get_img_dir(id)
        ltId, pos = id.split('--')
        
        # temp directory
        if tempDir is None:
            tempDir = os.path.join(out_path, 'temp')
        if not os.path.isdir(tempDir):
            os.makedirs(tempDir)

        if trajectory_dir is None:
            trajectory_dir = TRAJECTORY_DIR
                
        # read feature data
        if not feature_movie_dir is None and not feature is None:
            
            co = self.make_circle_offset(RADIUS)      
            
            cm = ColorMap()
            cr = cm.makeColorRamp(256, ["#FFFF00", "#FF0000"])
            
            # centers_PLT0023_11--ex2005_06_03--sp2005_04_19--tt17--c3_w00210_01.pkl 
            # traj_noF_densities_w00210_01.hdf5.pkl
            ending = '%05i_01.pkl' % int(pos)
            feature_index = FEATURES.index(feature)
                        
            # find the feature_files
            labteks_found = filter(lambda x: x[:len(ltId)]==ltId, os.listdir(trajectory_dir))
            if len(labteks_found) == 0:
                raise ValueError('%s not found in %s' % (ltId, self.base_dir))
            if len(labteks_found) > 1:
                raise ValueError('multiple labteks found')
            labtek_dir = os.path.join(trajectory_dir, labteks_found[0])
            
            feature_files = filter(lambda x: 'hist_tabFeatures' in x and x[-len(ending):]==ending, 
                                   os.listdir(labtek_dir))
            if len(feature_files)==0:
                print 'file not found ... skipping feature projection'
            else:
                fp = open(os.path.join(labtek_dir, feature_files[0]), 'r')
                feat = pickle.load(fp)
                fp.close()
                
                values = feat[0][:,feature_index]
                colors = [cm.getColorFromMap(x, cr, FEATURE_RANGE[feature][0], FEATURE_RANGE[feature][1])
                          for x in values.tolist()]
                #tvec = feat[1][0]
                #xvec = feat[1][1]
                #yvec = feat[1][2]
                markers = {}
                for i, track_coord in enumerate(feat[1]):                    
                    tvec, xvec, yvec = track_coord
                    
                    for t, x, y in zip(tvec, xvec, yvec):
                        if not t in markers:
                            markers[t] = []
                        markers[t].append((x, y, colors[i]))
                        
                #disp_features = dict(zip(tvec, [x,y,v for x,y,v in zip(xvec, yvec, values)]))
                
        
        # convert images
        print 'converting images ... '
        modDir = in_path.replace('(', '\(')
        modDir = modDir.replace(')', '\)')
        
        img_list = filter(lambda x: os.path.splitext(x)[-1].lower() in ['.png', '.tif', '.tiff', '.jpg'], 
                          os.listdir(in_path))
        img_list.sort()
        
        for i, img_name in enumerate(img_list):
            img = vigra.readImage(os.path.join(in_path, img_name))
            new_filename = os.path.join(tempDir, 
                                        '%s.jpg' % os.path.splitext(img_name)[0][2:])
            vigra.impex.writeImage(img.astype(np.dtype('uint8')), 
                                   new_filename)
            if not feature_movie_dir is None:
                f_temp_dir = os.path.join(feature_movie_dir, 'temp')
                if not os.path.isdir(f_temp_dir):
                    os.makedirs(f_temp_dir)
                new_filename = os.path.join(f_temp_dir, 
                                            'feature_projection_%s.png' % os.path.splitext(img_name)[0][2:])                
                
                # generate the color image                
                width = img.shape[0]
                height = img.shape[1]
                img_rgb = vigra.RGBImage((width, height))
                
                for channel in img_rgb.channelIter():
                    channel.copyValues(img)
                
                if i in markers:
                    for m in markers[i]:
                        x = m[0]
                        y = m[1]
                        color = [255 * ccc for ccc in m[2]]
                        
                        xvec = co[0] + x
                        yvec = co[1] + y
                        for i in range(3):
                            img_rgb[xvec, yvec, i] = color[i]
                            #print x, y, xvec, yvec, color[i]
                        #for c in color:
                        #    for xo, yo in co:
                        #        img_rgb[x, y, c] = color[c]
                vigra.impex.writeImage(img_rgb.astype(np.dtype('uint8')), 
                                       new_filename)
                
        #shell_command = "mito_convert -s png -o %s %s/*.tif"
        #shell_command %= (tempDir, modDir)
        #os.system(shell_command)

        # movie filename
        movieName = id
        if not gene is None:
            movieName += ('--%s' % gene)
        if not sirna is None:
            movieName += ('--%s' % sirna)
        movieName += '.avi'

        # make output directory
        if outDir is None:
            outDir = os.path.join(out_path, ltId)
        if not os.path.isdir(outDir):
            os.makedirs(outDir)

        # encode command
        try:
            print 'generating movie ... '
            encode_command = 'mencoder "mf://%s/*.jpg" -mf fps=3 -o %s -ovc xvid -oac copy -xvidencopts fixed_quant=2.5'
            encode_command %= (tempDir, os.path.join(outDir, movieName))
            print encode_command
            print 'movie generated: %s' % os.path.join(outDir, movieName)
            os.system(encode_command)
            
            if not feature_movie_dir is None and not feature is None:
                encode_command = 'mencoder "mf://%s/*.png" -mf fps=3 -o %s -ovc xvid -oac copy -xvidencopts fixed_quant=2.5'
                target_dir = feature_movie_dir
                if not os.path.exists(target_dir):
                    os.makedirs(target_dir)
                encode_command %= (f_temp_dir, os.path.join(target_dir, '%s_%s' % (feature.replace(" ", ""), movieName)))
                print encode_command
                print 'movie generated: %s' % os.path.join(target_dir, '%s_%s' % (feature.replace(" ", ""), movieName))
                os.system(encode_command)
                
        except:
            print 'MOVIE GENERATION FAILED'
            pass    
        
        # cleaning up temporary directory
        try:
            print 'cleaning %s' % tempDir
            shell_command = 'rm %s/*.jpg' % tempDir
            print shell_command
            os.system(shell_command)

            #shell_command = 'rmdir %s' % tempDir
            #os.system(shell_command)
            if not feature_movie_dir is None and not feature is None:
                shell_command = 'rm %s/*.png' % f_temp_dir
                print shell_command
                os.system(shell_command)
                #shell_command = 'rmdir %s' % f_temp_dir
                #os.system(shell_command)

        except:
            print 'temp directory not removed'
            pass

        
        print 'Done: %s' % id
        print
        
        return

   
if __name__ ==  "__main__":
    '''
    Pour le lancer, mail de Thomas:
    python make_movies_mito_cbio.py --out_path ./mini_out --pickle_file mini.pkl --feature_target_dir ./mini_overlay
    '''
    parser = OptionParser()
    parser.add_option("-i", "--in_path", dest="in_path",
                      help="base input directory")
    parser.add_option("-o", "--out_path", dest="out_path",
                      help="base output directory")
    parser.add_option("-p", "--pickle_file", dest="pickle_file",
                      help="file containing the names of ids to process")
    parser.add_option("-f", "--feature_target_dir", dest="feature_target_dir",
                      help="directory for the feature projection movies")
                      
    (options, args) = parser.parse_args()
    
    if (options.pickle_file is None or options.out_path is None):
        parser.error("incorrect number of arguments!")

    out_path = options.out_path
    if options.in_path is None:
        in_path = BASE_DIR
    else:
        in_path = options.in_path

    fp = open(options.pickle_file, 'r')
    movies = pickle.load(fp)
    fp.close()

    mm = MovieMaker(in_path)
    l=['signed turning angle',
 'movement type',
 'corrected straightness index',
 'ball number2',
 'norm convex hull area',
 'effective space length',
 'effective speed',
 'diffusion coefficient',
 'largest move',
 'ball number1', 'entropy1', 'entropy2']
    for category in l:
        if category !='mean squared displacement':
            ids = movies[category]
    
            cat_out_path = os.path.join(out_path, category.replace(' ', '_'))
            if not os.path.exists(cat_out_path):
                os.makedirs(cat_out_path)
    
            for i,id in enumerate(ids[:20]):
                print 'making movie number %i %s' %(i, id)
                mm.make_movie(id, cat_out_path, 
                              feature_movie_dir=options.feature_target_dir,
                              feature=category)
        