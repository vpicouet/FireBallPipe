#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 11:05:52 2018

@author: dvibert
"""

import os
import numpy as np
from astropy.table import Table

from FireBallIMO.PSFInterpoler.PSFImageHyperCube import PSFImageHyperCube
from FireBallIMO.PSFInterpoler.PSFInterpoler import PSFInterpoler
from FireBallIMO.PSFInterpoler.SkySlitMapping        import SkySlitMapping

from mapping import Mapping, load_pickle
from polyfit import PolyFit

try:
   import cPickle as pickle # for python 2 only CPickle is faster
except:
   import pickle

cloudpath = '/data/ownCloud/FIREBALL/'

sky2maskCube = os.path.join(cloudpath, 'InstrumentModel/PSFCubes/Baseline-April-2015/PSFCube-Baseline-April-2015-IMAGE-MASK.pkl')
sky2detCube = os.path.join(cloudpath, 'InstrumentModel/PSFCubes/Baseline-April-2015/PSFCube-Baseline-April-2015-IMAGE-DETECTOR.pkl')

mask_psf_cube = PSFImageHyperCube(filename = sky2maskCube)
mapping = SkySlitMapping(cubein = mask_psf_cube , sym_centroid_deg=[2,3], symmetric=True) 

det_psf_cube = PSFImageHyperCube(filename = sky2detCube)
psfi = PSFInterpoler(det_psf_cube, centroid_deg=[2,5,2]) 

mask='F1'
if mask == 'F1':
    slit_location = os.path.join(cloudpath, 'Target_selection/targets_F1.txt')
if mask == 'grid':
    slit_location = os.path.join(cloudpath, 'Target_selection/grid_mask.txt')
if mask == 'F2':
    slit_location = os.path.join(cloudpath, 'Target_selection/targets_F2.txt')
if mask == 'F3':
    slit_location = os.path.join(cloudpath, 'Target_selection/targets_F3.txt')
if mask == 'F4':
    slit_location = os.path.join(cloudpath, 'Target_selection/targets_F4.txt')

slitfile = slit_location
if 'targets_F1' in slit_location:
    print('Computing predicted location for F1')
    slits = Table.read(slitfile, format='ascii')
    # remove multi object slits
    idok = (slits['slit_length_right'] != 0) &  (slits['slit_length_left'] !=0)
    slits = slits[idok]
    xmask = slits['xmask'] + (slits['slit_length_right'] - slits['slit_length_left'])/2.
    ymask = slits['ymask'] + slits['offset']
    z = slits['z'] 
    internalCount = slits['Internal-count']

elif 'grid_mask' in slit_location:
    print('Computing predicted location for grid mask')
    slits =  Table.read(slitfile, format='ascii')
    xmask = slits['x']
    ymask = slits['y']
    #x = xmask
    #y = -ymask

elif ('F2' in slit_location) or ('F3' in slit_location) or ('F4' in slit_location) or ('Tilted' in slit_location):
    print('Computing predicted location for science mask')
    slits =  Table.read(slitfile, format='ascii')
    xmask = slits['xmm']
    ymask = slits['ymm']
    internalCount = slits['Internal-count']
    #mainID = slits['#main_id']
    if not ('Tilted' in slit_location):
        z = slits['Z'] 

## eventually rotate mask along pa
#if pa != pa0:
#    rotate_pa_mask(pa)
    
# convert mask axis convention to zemax/IMO convention
xzmx = xmask
yzmx = - ymask  

zinc_lines = np.array([0.20255, 0.20619, 0.21382])
w = zinc_lines[:,np.newaxis]            

xzmx = xzmx[np.newaxis,:]
yzmx = yzmx[np.newaxis,:]
   
    
# compute detector position
xsky, ysky = mapping.SlitToSkyPosition(w, xzmx, yzmx)
xdet_zmx, ydet_zmx = psfi.Centroid(w, xsky, ysky)

# convert zemax/IMO to fits detector
xdet, ydet = -ydet_zmx, xdet_zmx
            

# set mapping object for sky2mask
##################################
pfit =  mapping._cubei.symmetric_centroid_polyfit
m = Mapping(symmetric=True)
m.mapping = PolyFit(deg=pfit.deg, type='2d_woct', coeffs=pfit.coeffs, nbset=pfit.nbset)
ipfit =  mapping._cubei.inverse_symmetric_centroid_polyfit
m.inv_mapping = PolyFit(deg=ipfit.deg, type='2d_woct', coeffs=ipfit.coeffs, nbset=ipfit.nbset)

#######3 CHECK mapping object
xysky = m.inv_map(w, xzmx, yzmx)
xsky_ = xysky[0]
ysky_ = xysky[1]

print(np.all(xsky - xsky_ == 0.))
print(np.all(ysky - ysky_ == 0.))

# save mapping 
m.save('mapping-sky2mask-Baseline-April-2015.pkl')

test = Mapping('mapping-sky2mask-Baseline-April-2015.pkl')
xysky = test.inv_map(w, xzmx, yzmx)
xsky_ = xysky[0]
ysky_ = xysky[1]

print(np.all(xsky - xsky_ == 0.))
print(np.all(ysky - ysky_ == 0.))

# save xsky, ysky for check in python3
with  open('xsky.pkl','wb') as file_handle:            
    pickle.dump(xsky, file_handle)
with  open('ysky.pkl','wb') as file_handle:            
    pickle.dump(ysky, file_handle)



# set mapping object for sky2det
##################################
pfit =  psfi.centroid_polyfit
m = Mapping()
m.mapping = PolyFit(deg=pfit.deg, type='3d', coeffs=pfit.coeffs, nbset=pfit.nbset)
m.inv_mapping = None

#######3 CHECK mapping object
xydet = m.map(w, xsky, ysky)
xdet_ = xydet[0]
ydet_ = xydet[1]

print(np.all(xdet_zmx - xdet_ == 0.))
print(np.all(ydet_zmx - ydet_ == 0.))

# save mapping 
m.save('mapping-sky2det-Baseline-April-2015.pkl')

test = Mapping('mapping-sky2det-Baseline-April-2015.pkl')
xydet = test.map(w, xsky, ysky)
xdet_ = xydet[0]
ydet_ = xydet[1]

print(np.all(xdet_zmx - xdet_ == 0.))
print(np.all(ydet_zmx - ydet_ == 0.))

# save xdet_zmx, ydet_zmx for check in python3
with  open('xdet_zmx.pkl','wb') as file_handle:            
    pickle.dump(xdet_zmx, file_handle)
with  open('ydet_zmx.pkl','wb') as file_handle:            
    pickle.dump(ydet_zmx, file_handle)










# test in python 3 
#xsky = load_pickle('xsky.pkl')
#ysky = load_pickle('ysky.pkl')
#
#test = Mapping('mapping-sky2mask-Baseline-April-2015.pkl')
#xysky_ = test.inv_map(w, xzmx, yzmx)
#xsky_ = xysky_[0]
#ysky_ = xysky_[1]
#
#print(np.all(xsky - xsky_ == 0.))
#print(np.all(ysky - ysky_ == 0.))
#
#xdet_zmx = load_pickle('xdet_zmx.pkl')
#ydet_zmx = load_pickle('ydet_zmx.pkl')
#test = Mapping('mapping-sky2det-Baseline-April-2015.pkl')
#xydet_ = test.map(w, xsky, ysky)
#xdet_ = xydet_[0]
#ydet_ = xydet_[1]
#
#print(np.all(xdet_zmx - xdet_ == 0.))
#print(np.all(ydet_zmx - ydet_ == 0.))
