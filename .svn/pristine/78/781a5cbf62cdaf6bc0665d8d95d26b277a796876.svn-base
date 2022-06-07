#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 17:53:07 2018

@author: dvibert
"""
import numpy as np
from matplotlib import pyplot as plt

from astropy import coordinates
from astropy import units as u
from astropy.table import Table

from  guider2UV.guider2UV import Guider2UV
from guider2UV.MaskAstrometry import LocalScienceMaskProjector

from read_ds9_ctr import read_ds9_ctr

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


cloudpath = '/home/dvibert/ownCloud/FIREBALL/'

G2UV = Guider2UV(filename=cloudpath + 'TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/F3_180904.pkl')

G2UV_M31center = G2UV.copy()

Field_center = coordinates.SkyCoord(10.68*u.deg, 41.27*u.deg)
Field_rotation = -75*u.deg
Field_gamma = G2UV.FieldP.gamma

FieldP_M31center = LocalScienceMaskProjector(Field_center, Field_rotation, Field_gamma)

G2UV_M31center.FieldP = FieldP_M31center
print(G2UV_M31center)

# guider center in Ra/Dec
gc_local_coord = G2UV_M31center.guider_to_FieldLocal(np.array([[640.,540.]]), angle=False)
print(gc_local_coord)
G2UV_M31center.FieldP.local2world(gc_local_coord)



M31center_stars_radec = np.array([[11.125, 41.35]])

#M31center_stars_table = Table(M31center_stars_radec, names=[ 'RA', 'DEC'])
#M31center_stars_table['Internal count'] = [1]
#stars = [1]

stars_radec = coordinates.SkyCoord(M31center_stars_radec[:,0]*u.deg, M31center_stars_radec[:,1]*u.deg)
stars_pos_guider_pix = G2UV_M31center.SienceMask2guider(stars_radec, world=True, angle=False)
print(stars_pos_guider_pix)

ctrs = read_ds9_ctr('Calibration/Slits/F3.ctr')

select = np.arange(1)

# plot in guider pix
plot_stars(stars_pos_guider_pix, select, ctrs, 'Field mask 3 - M31 center star')

###################### 
G2UV_M31N206 = G2UV.copy()

Field_center = coordinates.SkyCoord(10.13*u.deg, 40.74*u.deg)
Field_rotation = -75*u.deg
Field_gamma = G2UV.FieldP.gamma

FieldP_M31N206 = LocalScienceMaskProjector(Field_center, Field_rotation, Field_gamma)

G2UV_M31N206.FieldP = FieldP_M31N206
print(G2UV_M31N206)

# guider center in Ra/Dec
gc_local_coord = G2UV_M31N206.guider_to_FieldLocal(np.array([[640.,540.]]), angle=False)
print(gc_local_coord)
G2UV_M31N206.FieldP.local2world(gc_local_coord)
