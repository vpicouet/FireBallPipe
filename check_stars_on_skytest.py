#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 08:46:34 2018

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
import os, glob
from astropy.table import Table
from matplotlib import pyplot as plt



cloudpath = '/home/dvibert/ownCloud/FIREBALL/'


target_stars = {'F2':np.array([[308.5,	137.8],
                               [740.7,	277.0],
                               [628.0,	615.5],
                               [506.4,	188.6]]),

                'F3':np.array([[341.7,	903.7],
                              [946.6,	743.4],
                               [50.6,	237.2]]),
	
                'F1': np.array([[989.5,	435.2],
                               [840.4,	298.6],
                                [942.5,	911.1]]),
	
        	        'F4': np.array([[370.5,	229.1],
                                [326.5,	76.9],
                                [727.5,	724.6],
                                [889.2,	685.1]]),
	
                'QSO MgII': np.array([[1024.7,	819.9],
                                      [687.3,	596.6],
                                      [847.1,	427.8]])
                }

##############################
## Field 2 180908
##############################

#path = cloudpath  + 'TestsFTS2018-Flight/E2E-AIT-Flight/SkyTest080918/'
path = cloudpath  + 'TestsFTS2018-Flight/E2E-AIT-Flight/SkyTest110918/TF_-161/'
#start = 22469325
#end   = 22530861

pa = -161.0

tab = Table.read(path + 'MergedCatalog.csv', format='csv')

pa_flag = (tab['ROTENC'] == pa)
tab = tab[pa_flag]

range_flag = (tab['FRAMESTA'] >= start ) & (tab['FRAMESTA'] <= end )
tab = tab[range_flag]

# sort by image number
framesta = np.unique(tab['FRAMESTA'])

plt.figure()
for j,n in enumerate(framesta):
    framesta_flag = (tab['FRAMESTA'] ==  n)
    line_idx = np.nonzero(framesta_flag)[0][0]
    line = tab[line_idx]
    # asume all lines build by Vincent contains same info from headers
    valid = np.zeros(8, dtype='bool')
    cx = np.zeros(8)
    cy = np.zeros(8)
    tx = np.zeros(8)
    ty = np.zeros(8)
    for i in range(8):
        istr = "{:1d}".format(i)
        valid[i] = (line['USE' + istr] & line['VALID' + istr] & 1) == 1
        cx[i] = line['CX' + istr]
        cy[i] = line['CY' + istr]
        tx[i] = line['TX' + istr]
        ty[i] = line['TY' + istr]
    # valid are actually 0->3 ($ stars)
    # checked targets are ok
    valid = np.arange(4)
    cx = cx[valid]
    cy = cy[valid]
    tx = tx[valid]
    ty = ty[valid]
    deltax = cx - tx
    deltay = cy - ty
    deltax_mean = deltax.mean()
    deltay_mean = deltay.mean()
    deltax_rms = np.sqrt(np.mean(np.square(deltax - deltax_mean)))
    deltay_rms = np.sqrt(np.mean(np.square(deltay - deltay_mean)))

    #plt.figure()
    plt.subplot(3,4,j+1)
    plt.scatter(tx, ty)
    plt.xlim([0,1280])
    plt.ylim([1080,0])
    q = plt.quiver(tx, ty, deltax, deltay)
    plt.quiverkey(q, .8,.8, 1, '1 pixel', color='k')
    legend =  "error mean in x,y {:.1f}, {:.1f} pixels\n".format(deltax_mean, deltay_mean)
    legend += "error rms in EL,CE {:.1f}, {:.1f} pixels".format(deltax_rms, deltay_rms)
    plt.text(50,800, legend)
    plt.xlabel('x pixels')
    plt.ylabel('y pixels')
    #plt.title('F2 stars throughfocus')

##############################
## Field 1 180911
##############################

#path = cloudpath  + 'TestsFTS2018-Flight/E2E-AIT-Flight/SkyTest080918/'
path = cloudpath  + 'TestsFTS2018-Flight/E2E-AIT-Flight/SkyTest110918/TF_119/'
pa = 119.0

tab = Table.read(path + 'MergedCatalog.csv', format='csv')

pa_flag = (tab['ROTENC'] == pa)
tab = tab[pa_flag]

#range_flag = (tab['FRAMESTA'] >= start ) & (tab['FRAMESTA'] <= end )
#tab = tab[range_flag]

# sort by image number
framesta = np.unique(tab['FRAMESTA'])

plt.figure()
for j,n in enumerate(framesta):
    framesta_flag = (tab['FRAMESTA'] ==  n)
    line_idx = np.nonzero(framesta_flag)[0][0]
    line = tab[line_idx]
    # asume all lines build by Vincent contains same info from headers
    valid = np.zeros(8, dtype='bool')
    cx = np.zeros(8)
    cy = np.zeros(8)
    tx = np.zeros(8)
    ty = np.zeros(8)
    for i in range(8):
        istr = "{:1d}".format(i)
        valid[i] = (line['USE' + istr] & line['VALID' + istr] & 1) == 1
        cx[i] = line['CX' + istr]
        cy[i] = line['CY' + istr]
        tx[i] = line['TX' + istr]
        ty[i] = line['TY' + istr]
    # valid are actually 0->2 (3 stars)
    # checked targets are ok
    valid = np.arange(3)
    cx = cx[valid]
    cy = cy[valid]
    tx = tx[valid]
    ty = ty[valid]
    deltax = cx - tx
    deltay = cy - ty
    deltax_mean = deltax.mean()
    deltay_mean = deltay.mean()
    deltax_rms = np.sqrt(np.mean(np.square(deltax - deltax_mean)))
    deltay_rms = np.sqrt(np.mean(np.square(deltay - deltay_mean)))

    #plt.figure()
    plt.subplot(3,4,j+1)
    plt.scatter(tx, ty)
    plt.xlim([0,1280])
    plt.ylim([1080,0])
    q = plt.quiver(tx, ty, deltax, deltay)
    plt.quiverkey(q, .8,.8, 3, '3 pixels', color='k')
    legend =  "error mean in x,y {:.1f}, {:.1f} pixels\n".format(deltax_mean, deltay_mean)
    legend += "error rms in EL,CE {:.1f}, {:.1f} pixels".format(deltax_rms, deltay_rms)
    plt.text(50,800, legend)
    plt.xlabel('x pixels')
    plt.ylabel('y pixels')
    #plt.title('F2 stars throughfocus')
    

#files = glob.glob(os.path.join(path, 'stack*.fits')).sort()
#filenumbers = np.array([int(os.path.basename(f).replace('stack','').replace('.fits','')) for f in files])
#select = (filenumbers >= start) & (filenumbers <= end)
   
#for f in files[select]:
#    with fits.open(f) as hdu:
#        h = hdu[0].header
#        if h['ROTENC'] != pa:
#            continue
#        for i in range(8):
#            valid = h['VALID{:1d}'.format(i)]
#            use = h['USE{:1d}'.format(i)]
#            tx = h['TX{:1d}'.format(i)]
#            ty = h['TY{:1d}'.format(i)]
#            flux = h['FLUX{:1d}'.format(i)]
#            sigmax = h['SIGMAX{:1d}'.format(i)]
#            sigmay = h['SIGMAY{:1d}'.format(i)]
#            cx = h['CX{:1d}'.format(i)]
#            cy = h['CY{:1d}'.format(i)]
#            





