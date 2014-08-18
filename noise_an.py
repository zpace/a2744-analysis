import pyfits as pf
import numpy as np
import pysao
import matplotlib.pyplot as plt
from sys import argv
from scipy.stats import anderson, skewtest, normaltest, kurtosistest

def noise(fname, x0 = 100, y0 = 100, maxrad = 30):
    from astroML.plotting import hist
    hdulist = pf.open(fname)
    im = hdulist[0].data
    #print np.mean(im), np.min(im), np.max(im)
    #print im[95:105, 95:105]
    # x0, y0 = 100, 100
    xi, yi = np.indices(im.shape)
    R = np.sqrt( (yi - int(y0))**2. + (xi - int(x0))**2. )
    phot_a = np.zeros(maxrad + 1)
    phot_a[0] = 0
    
    bmasked = im * ((R > maxrad) * (R < maxrad + 20.))
    bdata = bmasked.flatten()
    #print bdata[bdata != 0.]
    #print len(bdata[bdata != 0.])
    #print len(bdata)
    
    plt.subplot(3, 1, 1)
    hist(bdata[bdata != 0.], bins = 'blocks')
    plt.xlabel('Flux')
    plt.ylabel('(Bayesian Blocks)')
    plt.title('Noise')
    #plt.show()
    
    plt.subplot(3, 1, 2)
    hist(bdata[bdata != 0.], bins = 50)
    plt.xlabel('Flux')
    plt.ylabel('(50 bins)')
    #plt.title('Noise (50 bins)')
    #plt.show()
    
    plt.subplot(3, 1, 3)
    hist(bdata[bdata != 0.], bins = 'knuth')
    plt.xlabel('Flux')
    plt.ylabel('(Knuth\'s Rule)')
    #plt.title('Noise (Knuth\'s Rule)')
    plt.show()
    
    A2, crit, sig = anderson(bdata[bdata != 0.], dist = 'norm')
    print 'A-D Statistic:', A2
    print ' CVs \t  Sig.'
    print np.vstack((crit, sig)).T

    normality = normaltest(bdata[bdata != 0.])
    print 'Normality:', normality

    skewness = skewtest(bdata[bdata != 0.])
    print 'Skewness:', skewness

    kurtosis = kurtosistest(bdata[bdata != 0.])
    print 'Kurtosis:', kurtosis

    print 'Mean:', np.mean(bdata[bdata != 0.])
    print 'Median:', np.median(bdata[bdata != 0.])