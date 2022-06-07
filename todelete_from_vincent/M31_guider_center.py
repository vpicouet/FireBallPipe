#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 21:54:02 2018

@author: dvibert
"""

import numpy as np
from matplotlib import pyplot as plt

from astropy import coordinates
from astropy import units as u
from astropy.table import Table

from  guider2UV.guider2UV import Guider2UV
from guider2UV.MaskAstrometry import LocalScienceMaskProjector

cloudpath = '/home/dvibert/ownCloud/FIREBALL/'

G2UV = Guider2UV(filename=cloudpath + 'TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/F2_180907.pkl')

target_filename = 'Calibration/Targets/targets_F2.txt'
targets = Table.read(target_filename, format='ascii')


tx = np.array(targets['xmm']) 
ty = np.array(targets['ymm'])

g=G2UV.FieldP.pix2world(np.array([tx, ty]))

s31 = targets['Internal-count']=='31'
g[s31]

qso1 = targets['Internal-count']=='QSO1'
g[qso1]
