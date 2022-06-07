#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Twilight Bright Star

Created on Fri Sep 14 18:29:14 2018

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy import coordinates
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.table import Table
import matplotlib.patches as patches

from guider2UV.guider2UV import Guider2UV
from guider2UV.MaskAstrometry import LocalScienceMaskProjector
from read_ds9_ctr import read_ds9_ctr

cloudpath = '/home/dvibert/ownCloud/FIREBALL/'
path = cloudpath + 'TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/'

def plot_stars(guider_star_pos, select, ctrs, title):
    gc = np.array([640, 540]) # guider center

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    
    plt.scatter(guider_star_pos[0][select[0]], guider_star_pos[1][select[0]], 200, 
                marker="*", color='r', label='master star')
    plt.scatter(guider_star_pos[0][select[1:]], guider_star_pos[1][select[1:]], 200, 
                marker="*", color='g', label='guide stars')
    
    plt.xlim([0,1280])
    plt.ylim([1080,0])
    # plot guider center
    plt.plot(gc[0], gc[1], 'xg')
    plt.text(gc[0], gc[1], 'Gcenter')
        
    for i,s in  enumerate(select):
        plt.text(guider_star_pos[0][s]+30, guider_star_pos[1][s]-30, str(i+1))
    
    plt.xlabel('x guider pix (Dec axis)')
    plt.ylabel('y guider pix (-Ra axis)')
    plt.title(title)
    
    for c in ctrs:
        if c.shape[0] > 10:
            plt.plot(c[:,0], c[:,1],':k')        
    
    plt.show()

G2UV = Guider2UV(filename=path + 'F2_180907.pkl')

# centre guider: 
Ra = 245.674
Dec = 33.7333

Gcenter = coordinates.SkyCoord(Ra*u.deg, Dec*u.deg)
Gcenter_frame = Gcenter.skyoffset_frame(rotation=70*u.deg - 90*u.deg)

FOV_center = coordinates.SkyCoord(G2UV.FOV_center_guider_coord.lon, G2UV.FOV_center_guider_coord.lat, frame=Gcenter_frame)
FOV_center = FOV_center.transform_to('icrs')
print(FOV_center)
#FOV_center_guider_coord = coordinates.SkyCoord(G2UV.FOV_center_guider_coord.lon, 
#                                               G2UV.FOV_center_guider_coord.lat,frame=G2UV_TBS.GuiderP.localframe)

#G2UV_TBS = G2UV.copy()
G2UV_TBS = Guider2UV(Field_center = FOV_center,
                     Field_rotation = G2UV.FieldP.rotation.copy(), 
                     Field_gamma = G2UV.FieldP.gamma,
                     mask_rotation = G2UV.mask_rotation.copy(),
                     FOVcenter_guider_coord = G2UV.FOV_center_guider_coord.copy(),
                     guider_wcs = G2UV.GuiderP.w.copy() )
  

#G2UV_TBS.FieldP = LocalScienceMaskProjector(FOV_center, G2UV.FieldP.rotation, G2UV.FieldP.gamma)
#G2UV_TBS.FOV_center_guider_coord = coordinates.SkyCoord(G2UV.FOV_center_guider_coord.lon, 
#                                                        G2UV.FOV_center_guider_coord.lat,frame=G2UV_TBS.GuiderP.localframe)

print(G2UV_TBS)

#v2 CrB 245.62177, 33.7038 — mag 5.4
#v1 CrB 245.5893, 33.7989 — mag 5.2
#HIP80279 245.81093, 33.6976 --- mag 7.65  (UVmag 12.5)

TBS_stars_radec = np.array([[245.62177, 33.7038],
                            [245.5893, 33.7989],
                            [245.81093, 33.6976]])

TBS_stars = Table(TBS_stars_radec, names=[ 'RA', 'DEC'])
TBS_stars['Internal count'] = [1,2,3]

stars = [1,2,3]

#star_target_path = cloudpath + 'Target_selection/GuidingStars/'
#F1_stars = Table.read(star_target_path + "F1_guidingstars.fits", format='fits')
#mag = np.fmin(np.fmin(F1_stars['GAIA gband'].filled(99), F1_stars['SDSS gband'].filled(99)),
#        F1_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(TBS_stars['RA']*u.deg, TBS_stars['DEC']*u.deg)
#coords = [c for c in coords]

#guider_star_pos = G2UV_TBS.SienceMask2guider(coords, world=True, angle=False)
guider_star_pos = np.array([ np.array(G2UV_TBS.SienceMask2guider(c, world=True, angle=False)) for c in coords]).squeeze()

ctrs = read_ds9_ctr('Calibration/Slits/F2.ctr')

select = np.arange(3)

# plot in guider pix
plot_stars(guider_star_pos.T, select, ctrs, 'Field mask 2 - Twilight Bright Stars')
