#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 16:43:20 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy import coordinates
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.table import Table
from guider2UV import Guider2UV


path = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170901_perf3.0/sky_00_astrometry/"
F1_astrometry_filename = path + 'stack500513.new'
Field_center = coordinates.SkyCoord(32.19*u.deg, -5.688*u.deg)
G2UV = Guider2UV(F1_astrometry_filename, Field_center )
#G2UVang = Guider2UV(F1_astrometry_filename, Field_center )

################
################# registration  with SC-GUI01
### slit pos in mm

path_SCGUI01 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/"
 
F1_star_pattern_filename = path_SCGUI01  + 'F1_stack7971979.fits_table.fits'
F1_star_pattern_tab = Table.read(F1_star_pattern_filename)  
F1_star_pattern_tab.sort('xcentroid')
F1_star_pattern_pos = np.array([F1_star_pattern_tab['xcentroid'], F1_star_pattern_tab['ycentroid']]).T

starid = 3
slit_pos = np.array([0.5525186, -4.7601449]) #30


### register with new data
#star_UV_filename = path + 'starsinguiderwhenuvsource/stack1399407_table.fits' #F?
star_UV_filename = path_SCGUI01 + '/Guider/img/stack8074052.fits_table.fits' #F1
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T


G2UV.set_FOV_center(F1_star_pattern_pos, star_UV_pos, slit_pos, UVStar_id=3, world=False)
#FOV angular position in guider <SkyCoord (SkyOffsetICRS: rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
#    ( 250.39272815,  36.41856381)>): (lon, lat) in deg
#    ( 0.17572098,  0.0154012)>
#FOV pixel position in guider [array(1335.9526770494988), array(600.7896965721665)]

#G2UV.save(path + 'Guider2UV_F1.pkl')
#G2UV.save(path_SCGUI01 + 'Guider2UV_F1.new.pkl')

G2UV= Guider2UV(filename=path_SCGUI01  + 'Guider2UV_F1.new.pkl')

plt.figure()
plt.plot(F1_star_pattern_pos[:,0], F1_star_pattern_pos[:,1], '+')
plt.xlim([0,1400])
plt.ylim([1100, 0])
for i in  range(4):
    plt.text(F1_star_pattern_pos[i,0], F1_star_pattern_pos[i,1], str(i))

plt.plot(star_UV_pos[:,0], star_UV_pos[:,1], 'x')
for i in  range(3):
    plt.text(star_UV_pos[i,0], star_UV_pos[i,1], str(i))


plt.plot(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1],'o')
plt.grid()


    
    
# predict all slits with  another pattern
##########################################
path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"

new_F1_star_pattern_filename = path_SCGUI02  + 'Pattern/done/F1_stack2800289.fits_table.fits'

new_F1_star_pattern_tab = Table.read(new_F1_star_pattern_filename)  
new_F1_star_pattern_tab.sort('xcentroid')
new_F1_star_pattern_pos = np.array([new_F1_star_pattern_tab['xcentroid'], new_F1_star_pattern_tab['ycentroid']]).T

slit_pos = np.array([0.5525186, -4.7601449]) # slit 30
pred = G2UV.pattern_in_slit(starid, slit_pos, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred)
#[array([166.3051032 , 507.3500569 , 564.94180604, 894.13619392]), 
#array([609.33069497, 216.83893945, 936.5268878 , 549.3900237 ])]

slit_pos2 = np.array([ 3.7831, -2.1463]) #slit 38
pred2 = G2UV.pattern_in_slit(starid, slit_pos2, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred2)
#array([  427.21333233,   765.99792323,   812.50340864,  1140.19497765]), 
#array([ 293.84745642,  -86.46704138,  624.13171706,  248.28738213])]

slit_pos3 = np.array([4.376, 2.644]) # 39
pred3 = G2UV.pattern_in_slit(starid, slit_pos3, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred3)
#[array([  888.72746867,  1215.36282835,  1259.6109785 ,  1576.52166497]), 
#array([ 244.32217839, -125.01150783,  567.65426811,  202.43595081])]


slit_pos4 = np.array([-3.655,-.894]) #18
pred4 = G2UV.pattern_in_slit(starid, slit_pos4, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred4)
#[array([ 553.16188681,  876.162271  ,  945.51900347, 1258.12583329]), 
#array([1012.89897261,  622.23598526, 1323.87584448,  939.0695285 ])]

slit_pos5 = np.array([-3.655,-.894]) #18
pred5 = G2UV.pattern_in_slit(starid, slit_pos5, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred5)
#[array([ 553.16188681,  876.162271  ,  945.51900347, 1258.12583329]), 
#array([1012.89897261,  622.23598526, 1323.87584448,  939.0695285 ])]

this_pred = pred4
plt.figure()
plt.plot(new_F1_star_pattern_pos[:,0], -new_F1_star_pattern_pos[:,1], '+')
#plt.gca().invert_xaxis()
plt.axis('equal')
plt.xlim([0,1400])
plt.ylim([-1100, 0])
for i in  range(4):
    plt.text(new_F1_star_pattern_pos[i,0], -new_F1_star_pattern_pos[i,1], str(i))
plt.plot(G2UV.FOV_center_guider_pix[0], -G2UV.FOV_center_guider_pix[1],'o')
plt.grid()
plt.plot(this_pred[0], -this_pred[1],'*r')
for i in  range(4):
    plt.text(this_pred[0][i], -this_pred[1][i], str(i))

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
 

