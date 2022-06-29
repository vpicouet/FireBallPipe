# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy.table import Table
from astropy import units as u
import matplotlib.pyplot as plt
from scipy      import optimize

from MaskAstrometry import LocalGuiderMaskProjector, LocalScienceMaskProjector


f3_astrometry_filename = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170901_perf3.0/sky_00_astrometry/stack475788.new"
hdu = fits.open(f3_astrometry_filename)
w  = wcs.WCS(hdu[0].header)

#F3_converter = LocalGuiderMaskProjector(w)
#
#m13tab = Table.read('stack475788.starpos.fits')
#
#
#local = F3_converter.pix2local(np.array([m13tab['xcentroid'],m13tab['ycentroid']]).T)
#plt.figure()
#plt.plot(local.lon, local.lat, 'o')
#
#pix = F3_converter.local2pix(local)
#plt.figure()
#plt.plot(m13tab['xcentroid'], m13tab['ycentroid'],'+r')
#plt.plot(pix[0], pix[1], 'x')

path = "/data/FireBall/alignment"

#tag6
thispath = path + '/tag6/'
files = ['stack2279546.fits_table.fits',
         'stack2278216.fits_table.fits',
            'stack2277818.fits_table.fits',
            'stack2276284.fits_table.fits',
            'stack2273428.fits_table.fits',]
            # 'stack2271620.fits_table.fits']

#tag7
thispath = path + '/tag7/'
files = [#'stack2283451.fits_table.fits',
         'stack2282438.fits_table.fits',
         'stack2281281.fits_table.fits',
         ]

#tag8 
thispath = path + '/tag8/'
files = ['stack2292466.fits_table.fits',
        'stack2292015.fits_table.fits',
        'stack2291415.fits_table.fits',
        'stack2289503.fits_table.fits',
        # 'stack2288725.fits_table.fits'
         ]

#tag 9
thispath = path + '/tag9/'
files = ['stack2294939.fits_table.fits',
         'stack2294197.fits_table.fits',
         'stack2293593.fits_table.fits'
         ]

#tag10
thispath = path + '/tag10/'
files = ['stack2303021.fits_table.fits',
         'stack2302286.fits_table.fits',
         'stack2301250.fits_table.fits',
         ]

#tag11
thispath = path + '/tag11/'
files = ['stack2306383.fits_table.fits',
         'stack2305293.fits_table.fits',
         'stack2304556.fits_table.fits',]


rotstar = []
for f in files:
    tab = Table.read(thispath + f)[ 'xcentroid','ycentroid']
    tab.sort('xcentroid')
    rotstar.append( tab[0].data )

rotstar = np.array(rotstar)
plt.figure()
plt.plot(rotstar['xcentroid'], rotstar['ycentroid'],'o')

# find center rot


x = rotstar['xcentroid'] 
y = rotstar[ 'ycentroid'] 

def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f_2(c):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

center_estimate = 1280, 540
center_2, ier = optimize.leastsq(f_2, center_estimate)

xc_2, yc_2 = center_2
print(xc_2, yc_2)

plt.plot(xc_2, yc_2, '+r')
Ri_2       = calc_R(*center_2)
R_2        = Ri_2.mean()
residu_2   = sum((Ri_2 - R_2)**2)

print(R_2)

theta = np.linspace(0,2*np.pi,100)
xx = xc_2 + R_2*np.cos(theta)
yy = yc_2 + R_2*np.sin(theta)

plt.plot(xx,yy,'xg')
plt.xlim(0,1400)
plt.ylim(0,1100)