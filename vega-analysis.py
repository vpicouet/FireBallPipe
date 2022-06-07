#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 12:03:45 2018

@author: dvibert
"""

from Calibration.focustest import *
from astropy.io import fits
from astropy.visualization import MinMaxInterval, SqrtStretch, AsinhStretch, AsymmetricPercentileInterval, ImageNormalize     

from maptplotlib import pyplot as plt
import numpy as np

path = '/data/ownCloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180611/'

# stack dark img12-21
#dark, darkname = stackImages(path, all=False, DS=0, function = 'mean', numbers=np.arange(12,22).astype(int), save=True, name="dark")
#
##stack an sub dark
#stackImages(path, all=False, DS=dark, function = 'mean', numbers=np.arange(2,12).astype(int), save=True, name="subdark-vega")

dark_name = 'image-000012-000021-dark-000012-000021-dark-stack.fits'
with fits.open(path + dark_name) as f:
    dark = f[0].data
    dark_header = f[0].header

plt.figure()
norm = ImageNormalize(dark, interval=AsymmetricPercentileInterval(10, 90), stretch=AsinhStretch())   
plt.imshow(dark, cmap='cubehelix', norm=norm, origin='lower')
plt.colorbar()

vega_name = 'image-000002-000011-subdark-vega-stack.fits'
with fits.open(path + vega_name) as f:
    vega = f[0].data
    vega_header = f[0].header

increase = vega/dark*100.

plt.figure()
norm = ImageNormalize(vega, interval=AsymmetricPercentileInterval(10, 90), stretch=AsinhStretch())   
plt.imshow(vega, cmap='cubehelix', norm=norm, origin='lower')

# stack lambda

#vega_lambda_sum = vega[:,1069:2125].sum(axis=1)
vega_lambda_sum = vega[:,1069:2001].mean(axis=1)
plt.figure()
plt.plot(vega_lambda_sum)

increase = vega[:,1069:2001]/dark[:,1069:2001]*100.
increase_sum = increase.mean(axis=1)


plt.figure()
plt.plot(increase_sum)