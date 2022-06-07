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

from guider2UV import Guider2UV


path = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/perf3.0-09-01/sky_00_astrometry/"
F1_astrometry_filename = path + 'stack500513.new'
Field_center = coordinates.SkyCoord(32.19*u.deg, -5.688*u.deg)
G2UV = Guider2UV(F1_astrometry_filename, Field_center )
#G2UVang = Guider2UV(F1_astrometry_filename, Field_center )


#################
################# registration  with SC-GUI01
### slit pos in mm

path_SCGUI01 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/"
 
F1_star_pattern_filename = path_SCGUI01  + 'F1_stack7971979.fits_table.fits'
F1_star_pattern_tab = Table.read(F1_star_pattern_filename)  
F1_star_pattern_tab.sort('xcentroid')
F1_star_pattern_pos = np.array([F1_star_pattern_tab['xcentroid'], F1_star_pattern_tab['ycentroid']]).T

starid = 3
slit_pos = np.array([0.5525186, -4.7601449]) #30
slit_pos_ang = coordinates.SkyCoord(32.202194*u.deg, -5.800675*u.deg) #30


### register with new data
#star_UV_filename = path + 'starsinguiderwhenuvsource/stack1399407_table.fits' #F?
star_UV_filename = path_SCGUI01 + '/Guider/img/stack8074052.fits_table.fits' #F1
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T


G2UV.set_FOV(F1_star_pattern_pos, star_UV_pos, slit_pos, UVStar_id=3, world=False)
G2UV.set_FOV(F1_star_pattern_pos, star_UV_pos, slit_pos_ang, UVStar_id=3, world=True)

#G2UV.save(path_SCGUI01 + "test.pkl")
#
#G2UV.restore(path_SCGUI01 + "test.pkl")
#

plt.figure()
plt.plot(F1_star_pattern_pos[:,0], F1_star_pattern_pos[:,1], '+')
ax = plt.gca()
ax.invert_yaxis()
plt.xlim([0,1400])
plt.ylim([1100, 0])
for i in  range(4):
    plt.text(F1_star_pattern_pos[i,0], F1_star_pattern_pos[i,1], str(i))

plt.plot(star_UV_pos[:,0], star_UV_pos[:,1], 'x')
for i in  range(3):
    plt.text(star_UV_pos[i,0], star_UV_pos[i,1], str(i))


plt.plot(G2UV.FOV_guider_pix[0], G2UV.FOV_guider_pix[1],'o')
plt.grid()


    
    
# predict all slits with  another pattern
##########################################
path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI02/"

new_F1_star_pattern_filename = path_SCGUI02  + 'Pattern/done/F1_stack2800289.fits_table.fits'

new_F1_star_pattern_tab = Table.read(new_F1_star_pattern_filename)  
new_F1_star_pattern_tab.sort('xcentroid')
new_F1_star_pattern_pos = np.array([new_F1_star_pattern_tab['xcentroid'], new_F1_star_pattern_tab['ycentroid']]).T

slit_pos = np.array([0.5525186, -4.7601449]) # slit 30
pred = G2UV.pattern_in_slit(starid, slit_pos, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred)

slit_pos2 = np.array([ 3.7831, -2.1463])
pred2 = G2UV.pattern_in_slit(starid, slit_pos2, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred2)
#array([  427.21333233,   765.99792323,   812.50340864,  1140.19497765]), 
#array([ 293.84745642,  -86.46704138,  624.13171706,  248.28738213])]

slit_pos3 = np.array([4.376, 2.644]) # 
pred3 = G2UV.pattern_in_slit(starid, slit_pos3, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred3)
#[array([  888.72746867,  1215.36282835,  1259.6109785 ,  1576.52166497]), 
#array([ 244.32217839, -125.01150783,  567.65426811,  202.43595081])]


slit_pos4 = np.array([-3.655,-.894]) #
pred4 = G2UV.pattern_in_slit(starid, slit_pos4, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred4)

slit_pos5 = np.array([-3.655,-.894]) #
pred5 = G2UV.pattern_in_slit(starid, slit_pos5, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred5)


this_pred = pred4
plt.figure()
plt.plot(new_F1_star_pattern_pos[:,0], -new_F1_star_pattern_pos[:,1], '+')
#plt.gca().invert_xaxis()
plt.axis('equal')
plt.xlim([0,1400])
plt.ylim([-1100, 0])
for i in  range(4):
    plt.text(new_F1_star_pattern_pos[i,0], -new_F1_star_pattern_pos[i,1], str(i))
plt.plot(G2UV.FOV_guider_pix[0], -G2UV.FOV_guider_pix[1],'o')
plt.grid()
plt.plot(this_pred[0], -this_pred[1],'*r')
for i in  range(4):
    plt.text(this_pred[0], -this_pred[1], str(i))

# predict all slits same pattern
################################
target_file = "/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/targets_F1.txt"
target_tab = Table.read(target_file, format='ascii')

all_pred = np.zeros((4,2,len(target_tab)))
for target in range(len(target_tab)):
    this_slit_pos = np.array([target_tab['col4'][target], target_tab['col5'][target]])
    pred = G2UV.pattern_in_slit(starid, this_slit_pos, world=False)
    all_pred[:, 0, target] = pred[0]
    all_pred[:, 1, target] = pred[1]


all_pred = np.zeros((4,2,len(target_tab)))
for target in range(len(target_tab)):
    this_slit_pos = np.array([target_tab['col4'][target], target_tab['col5'][target]])
    pred = G2UV.pattern_in_slit(starid, this_slit_pos, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
    all_pred[:, 0, target] = pred[0]
    all_pred[:, 1, target] = pred[1]
 

#################################
#
# predict guiding stars
#
##################################
star_target_path = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/GuidingStars/'
F1_stars = Table.read(star_target_path + "F1_guidingstars.txt", format='ascii')

coords = coordinates.SkyCoord(F1_stars['RA']*u.deg, F1_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)


plt.figure()
plt.plot(guider_star_pos[0], guider_star_pos[1], '*')
#plt.xlim([0,1400])
#plt.ylim([0,-1100])

#write table
F1_stars['Xguider'] = guider_star_pos[0]
F1_stars['Yguider'] = guider_star_pos[1]

F1_stars.write(star_target_path + "F1_guidingstars.txt", format='ascii', overwrite=True)
