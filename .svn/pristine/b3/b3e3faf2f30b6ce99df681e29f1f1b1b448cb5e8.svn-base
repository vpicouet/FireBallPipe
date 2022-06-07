#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 15:29:27 2018

@author: Vincent
"""


import numpy as np
from astropy import wcs
from astropy.io import fits
import sys

def load_wcs_from_file(fitsname,wcsname):
    hdulist = fits.open(fitsname)
    w = wcs.WCS(wcsname)#w = wcs.WCS(hdulist[0].header)
    w.wcs.print_contents()
    pixcrd = np.array([[0, 0], [24, 38], [45, 98]], np.float_)

    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)
    print(world)

    # Convert the same coordinates back to pixel coordinates.
    pixcrd2 = w.wcs_world2pix(world, 1)
    print(pixcrd2)
    
    header = w.to_header()
    #print (hdulist[0].header)
    #print (header)
    hdulist[0].header.extend(header)#fits.PrimaryHDU(header=header)
    print (hdulist[0].header)

    hdulist.writeto(fitsname[:-5]+'_wcs'+'.fits',overwrite=True)
    return hdulist


if __name__ == '__main__':
    print('''\n\n\n\n      START IMAGE ANALYSIS \n\n\n\n''')

    #wcsname = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018/wcs.fits'
    #fitsname = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018/stack8103716_pa+078_2018-06-11T06-21-24.fits'
    wcsname = sys.argv[-1]
    fitsname = sys.argv[-2]
    print ('fits=   ',fitsname,'\n wcs=     ',wcsname)
    a = load_wcs_from_file(fitsname,wcsname)
