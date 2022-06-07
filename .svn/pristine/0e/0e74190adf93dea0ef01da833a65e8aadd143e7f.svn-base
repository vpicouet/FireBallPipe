#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 17:12:11 2018

@author: dvibert
"""

from os import path

from Calibration.mapping import Mapping

cloudpath = '/data/ownCloud/FIREBALL/'
pklpath = 'TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/Detector_Mask_mappings/'

pkls = {'F1': 'mapping-mask-det-180612-F1.pkl', 
        'F2': 'mapping-mask-det-180612-F2.pkl',
        'F3': 'mapping-mask-det-180612-F3.pkl',
        'F4': 'mapping-mask-det-180612-F4.pkl'}

for val in pkls.values():
    pkl_file = path.join(cloudpath, pklpath, val)
    m = Mapping(pkl_file)
    
    name, ext = path.splitext(pkl_file)
    outname = name + '_p2' + ext
    m.save(path.join(outname))

#check 
test = Mapping(outname)