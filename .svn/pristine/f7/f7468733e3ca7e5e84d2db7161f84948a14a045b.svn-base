#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 16:17:45 2018

@author: dvibert
"""

from Calibration.focustest import *
from astropy.io import fits
from astropy.visualization import MinMaxInterval, SqrtStretch, AsinhStretch, AsymmetricPercentileInterval, ImageNormalize     

from matplotlib import pyplot as plt
import numpy as np
import glob



path = '/data/ownCloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018/'

cloudpath = '/data/ownCloud/FIREBALL/'

filenames = glob.glob(path + '*.fits') 

#name = filenames[0]
for name in filenames:
    test = Focus(filename = name, threshold = [10], fwhm = [5,7,9,12.5,15, 17], 
                 quick=False, reversex=False, plot=True, figsize=12, cloudpath=cloudpath, 
                 shape='gaussian', min_fwhm=2.5, cut=[50,99] )#, HumanSupervision=True)


#name = path  + 'stack8103079_pa+159_2018-06-11T06-18-16.fits'
#test = Focus(filename = name, threshold = [10], fwhm = [5,7,9,12.5,15, 17], 
#                 quick=False, reversex=False, plot=True, figsize=12, cloudpath=cloudpath, 
#                 shape='gaussian', min_fwhm=2.5, cut=[50,99] )#, HumanSupervision=True)
#
#test.plot_FWHM_guider('test', figsize=12, cut=[50,99])
