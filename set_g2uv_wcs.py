#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 22:04:13 2018

@author: dvibert
"""
from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from guider2UV.guider2UV import Guider2UV
from guider2UV.MaskAstrometry import LocalGuiderMaskProjector

cloudpath = '/home/dvibert/ownCloud/FIREBALL/'
pklpath = cloudpath + 'TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/'
wcspath = cloudpath + 'TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018_NEW/OpenCluster/'

wcs_files  = {'F1': 'stack8102222_pa+119_2018-06-11T06-13-59_wcs.fits',
              'F2': 'stack8099885_pa-161_2018-06-11T06-02-18_wcs.fits',
              'F2': 'stack8099885_pa-161_2018-06-11T06-02-18_wcs.fits',
              'F3': 'stack8100949_pa-121_2018-06-11T06-07-37_wcs.fits' ,
              'F4': 'stack8103027_pa+159_2018-06-11T06-18-01_wcs.fits' }

g2uv_files = {'F1': 'F1_180907.pkl',
              'F2': 'F2_180907.pkl',
              'F3': 'F3_180904.pkl',
              'F4': 'F4_180904.pkl'}

for f in ['F1','F2','F3','F4']:
    G2UV = Guider2UV(filename = pklpath + g2uv_files[f]) 
    print(G2UV.GuiderP.w)
    with fits.open(wcspath + f + '_wcs/' +  wcs_files[f]) as hdu:
        h = hdu[0].header
    w = WCS(h)
    scale = proj_plane_pixel_scales(w)
    del h['CD1_1']
    del h['CD1_2']
    del h['CD2_2']
    del h['CD2_1']
    h['CDELT1'] = scale[0]
    h['CDELT2'] = scale[1]
    guider_wcs = WCS(h)        
    #G2UV.GuiderP = LocalGuiderMaskProjector(guider_wcs)
    print(guider_wcs)