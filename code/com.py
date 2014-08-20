import pyfits as pf
import numpy as np
import pysao
import matplotlib.pyplot as plt
from sys import argv

def com(fname, maxrad = 10):
    """
        find centroid of a star, given a small thumbnail with the star close to the center
        
        fname: name of file to use including .fits extension
        maxrad: maximum number of pixels around the center to integrate
        """
    hdulist = pf.open(fname)
    im = hdulist[0].data
    #print np.mean(im), np.min(im), np.max(im)
    #print im[95:105, 95:105]
    x0, y0 = 100, 100
    
    xi, yi = np.indices(im.shape)
    R = np.sqrt( (yi - y0)**2. + (xi - x0)**2. )
    
    im = im * (R < maxrad)
    M = np.sum(im)
    x_ctr = np.sum(im * xi) / M
    y_ctr = np.sum(im * yi) / M
    return x_ctr, y_ctr

def com_dir(dirname, maxrad = 15):
    import glob
    """
        wrapper around com for reading in all files in a given directory and computing all COMs at once
        
        dirname: string of directory location of all the .fits files to be used
        maxrad: maximum number of pixels around the center to integrate, passes to com_dir
        """
    flist = glob.glob(dirname + '/*.fits')
    flist_short = [item[len(dirname)+1:-6] for item in flist]
    #print flist_short
    
    x_ctr_agg, y_ctr_agg = 0., 0.
    
    for i, filtername in enumerate(flist):
        x_ctr, y_ctr = com(filtername, maxrad = maxrad)
        #print 'COM of %s in filter %s is at %.2f, %.2f' % (dirname, flist_short[i], x_ctr, y_ctr)
        x_ctr_agg += x_ctr
        y_ctr_agg += y_ctr
    print 'Computed average COM of %s is %.2f, %.2f \n' % (dirname, x_ctr_agg/len(flist), y_ctr_agg/len(flist)) #remember to reverse these if interfacing with pyfits!