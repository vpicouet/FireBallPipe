#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 16:43:20 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy import units as u
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib import pyplot as plt
from astropy.table import Table
import matplotlib.patches as patches

from guider2UV import Guider2UV


path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/'
F2_astrometry_filename = path + 'Distortion/F2_stack486607.new'
Field_center = coordinates.SkyCoord(253.0624*u.deg, 34.9699*u.deg)
G2UV = Guider2UV(F2_astrometry_filename, Field_center, Field_rotation=(90-20)*u.deg )
#G2UVang = Guider2UV(F2_astrometry_filename, Field_center , Field_rotation=20*u.deg )


#################
################# registration  with SC-GUI01
### slit pos in mm

#path_SCGUI01 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/Pattern/done/"
# 
#F2_star_pattern_filename = path  + 'F2_stack7899665.fits_table.fits'
#
#F2_star_pattern_tab = Table.read(F2_star_pattern_filename)  
#F2_star_pattern_tab.sort('xcentroid')
#F2_star_pattern_pos = np.array([F2_star_pattern_tab['xcentroid'], F2_star_pattern_tab['ycentroid']]).T
#
#starid = 3
#slit_pos = np.array([0.5525186, -4.7601449]) #30


### register with new data

path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"

F2_star_pattern_filename = path_SCGUI02  + 'Pattern/done/F2_stack3060946.fits_table.fits'

F2_star_pattern_tab = Table.read(F2_star_pattern_filename)  
F2_star_pattern_tab.sort('xcentroid')
F2_star_pattern_pos = np.array([F2_star_pattern_tab['xcentroid'], F2_star_pattern_tab['ycentroid']]).T

slit_pos = np.array([0.32258, -4.0767]) #30

star_UV_filename = path_SCGUI02 + 'F2_stack3077872.fits_table.fits' #F2
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T


G2UV.set_FOV_center(F2_star_pattern_pos, star_UV_pos, slit_pos, UVStar_id=3, Stars_id=[0,1,2], world=False)
#FOV angular position in guider <SkyCoord (SkyOffsetICRS: rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
#    ( 250.39837843,  36.42897653)>): (lon, lat) in deg
#    ( 0.1831965, -0.01420325)>
#FOV pixel position in guider [array(1358.7625953703882), array(471.2734739991174)]

#G2UV.save(path + 'Guider2UV_F2.pkl')
G2UV.save(path + 'Guider2UV_F2.new.pkl')

G2UV= Guider2UV(filename=path + 'Guider2UV_F2.new.pkl')


# predictions3
starid = 3
pred = G2UV.pattern_in_slit(starid, slit_pos, world=False)
print(pred)
#[array([ 264.05437916,  605.78836908,  655.09592817,  986.94901275]), 
#array([ 510.27128209,  118.65678102,  837.3107609 ,  449.17732682])]

slit_pos1 = np.array([1.7145, -1.06426]) # 32
starid = 3
pred1 = G2UV.pattern_in_slit(starid, slit_pos1, world=False)
print(pred1)
#[array([  561.67622515,   898.61676342,   937.77163519,  1265.76327428]),
#array([ 378.91810739,  -14.26427139,  701.06725506,  311.6064834 ])]

slit_pos2 = np.array([]) 
pred2 = G2UV.pattern_in_slit(starid, slit_pos2, world=False)
print(pred2)


#slit_pos3 = np.array([]) # 
#pred3 = G2UV.pattern_in_slit(starid, slit_pos3, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
#print(pred3)
#
#
#slit_pos4 = np.array([]) #
#pred4 = G2UV.pattern_in_slit(starid, slit_pos4, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
#print(pred4)
#
#slit_pos5 = np.array([]) #
#pred5 = G2UV.pattern_in_slit(starid, slit_pos5, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
#print(pred5)







this_pred = pred1
plt.figure()
plt.plot(F2_star_pattern_pos[:,0], -F2_star_pattern_pos[:,1], '+')
#plt.gca().invert_xaxis()
plt.axis('equal')
plt.xlim([0,1400])
plt.ylim([-1100, 0])
for i in  range(4):
    plt.text(F2_star_pattern_pos[i,0], -F2_star_pattern_pos[i,1], str(i))
plt.plot(G2UV.FOV_center_guider_pix[0], -G2UV.FOV_center_guider_pix[1],'o')
plt.grid()
plt.plot(this_pred[0], -this_pred[1],'*r')
for i in  range(4):
    plt.text(this_pred[0][i], -this_pred[1][i], str(i))


plt.figure()
plt.plot(F1_star_pattern_pos[:,0], -F1_star_pattern_pos[:,1], '+')
plt.xlim([0,1400])
plt.ylim([0,-1100])
for i in  range(4):
    plt.text(F1_star_pattern_pos[i,0], -F1_star_pattern_pos[i,1], str(i))

plt.plot(star_UV_pos[:,0], -star_UV_pos[:,1], 'x')
for i in  range(3):
    plt.text(star_UV_pos[i,0], -star_UV_pos[i,1], str(i))


plt.plot(G2UV.FOV_center_guider_pix[0], -G2UV.FOV_center_guider_pix[1],'o')
plt.grid()




