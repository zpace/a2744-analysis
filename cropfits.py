import pyfits as pf
import numpy as np
import pysao
from sys import argv

def cropfits(dirname, xctr, yctr, xbdr, ybdr):
    files = ['105_aligned.fits', '125_aligned.fits', '140_aligned.fits', '160_aligned.fits', '435_aligned.fits', '606_aligned.fits', '814_aligned.fits']
    shortfiles = ['435c.fits', '606c.fits', '814c.fits', '105c.fits', '125c.fits', '140c.fits', '160c.fits']
    for i, fname in enumerate(files):
        hdulist = pf.open(fname)
        im = hdulist[0].data
        im = im[yctr - ybdr:yctr + ybdr + 1, xctr - xbdr:xctr + xbdr + 1]
        hdulist[0].data = im
        #print type(hdulist[0].data)
        outname = dirname + '/' + shortfiles[i]
        hdulist.writeto(outname, clobber = True)
        print outname, 'successfully cropped!'