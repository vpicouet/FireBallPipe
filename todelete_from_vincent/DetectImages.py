#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 22:18:58 2018

@author: Vincent
"""
import sys
#sys.path.append('/Users/Vincent/Documents/FireBallPipe')
from Calibration.focustest import Focus

for path in sys.argv[1:]:
    print('path = ', path)
#    F = Focus(filename = path, HumanSupervision=False, source='Zn', shape='holes', windowing=False, peak_threshold=50,quick=False, threshold = [7], fwhm = [10], plot=True)
    F = Focus(filename = path, HumanSupervision=False, source='Zn', shape='holes', windowing=False, peak_threshold=50,quick=True, plot=True)

#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/XYCalib180904/ThermalDriftGuiderCenter/stack20149520.fits'
#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/XYCalib180904/F4/F4Guider/GS2_s14/stack20719325.fits'#'/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/XYCalib180904/F4/F4Guider/GS2_s14/stack20730132.fits'
#F = Focus(filename = path, HumanSupervision=False,quick=False, threshold = [7], fwhm = [12], source='Zn', shape='holes', windowing=False, peak_threshold=50,plot=True)
##
# 
