#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 15:14:30 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy import units as u
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from astropy.table import Table
from astroML.crossmatch import crossmatch
from guider2UV import Guider2UV

path_SCGUI01 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/"

#G2UV= Guider2UV(filename=path_SCGUI01  + 'Guider2UV_F1.pkl')
G2UV= Guider2UV(filename=path_SCGUI01  + 'Guider2UV_F4.pkl')


corners = np.array([[1, 1], [1280,1080]])
corner_angles = G2UV.GuiderP.w.all_pix2world(corners, 1)

field = corner_angles[1,:]-corner_angles[0,:]
cdelt = field/np.array([1208,1080])


# distortion map

x = np.linspace(1,1280,21)
y = np.linspace(1,1080,21)
xx,yy = np.meshgrid(x, y)

xxd, yyd = G2UV.GuiderP.w.sip_pix2foc(xx, yy, 1)
xxd -=  - G2UV.GuiderP.w.wcs.crpix[0]
yyd -=  - G2UV.GuiderP.w.wcs.crpix[1]

plt.figure()
qv = plt.quiver(xx, yy, xxd-xx, yyd-yy, scale=500 )
plt.quiverkey(qv, .5, -.1, 20, "20 pix")
plt.xlim([0,1280])
plt.ylim([1080,0])



plt.figure()
for i in range(x.size):
    plt.plot(xx[i,:], yy[i,:],':g' )
    plt.plot(xx[:,i], yy[:,i], ':g' )
    plt.plot(xxd[i,:], yyd[i,:],'b' )
    plt.plot(xxd[:,i], yyd[:,i], 'b' )
    