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
import matplotlib.pyplot as plt
from astropy.table import Table

from guider2UV import Guider2UV


path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/'
F3_astrometry_filename = path + 'Distortion/F3_stack475788.new'

Field_center = coordinates.SkyCoord(352.3424*u.deg, 0.21245*u.deg)
G2UV = Guider2UV(F3_astrometry_filename, Field_center )
#G2UVang = Guider2UV(F2_astrometry_filename, Field_center  )


#################
################# registration  with SC-GUI01
### slit pos in mm

path_SCGUI01 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/"
 
F3_star_pattern_filename = path  + '/Pattern/done/F3_stack7903393.fits_table.fits'

F3_star_pattern_tab = Table.read(F3_star_pattern_filename)  
F3_star_pattern_tab.sort('xcentroid')
F3_star_pattern_pos = np.array([F3_star_pattern_tab['xcentroid'], F3_star_pattern_tab['ycentroid']]).T

starid = 3
slit_pos = np.array([2.7178, -6.0325]) #30


### register with new data
#star_UV_filename = path + 'starsinguiderwhenuvsource/stack1399407_table.fits' #F?
star_UV_filename = path_SCGUI01 + 'Guider/guidingStarsWhenInSlit/stack1340907.fits_table.fits' #F3
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T


G2UV.set_FOV_center(F3_star_pattern_pos, star_UV_pos, slit_pos, UVStar_id=3, Stars_id=[0,1], world=False)

# => wrong star_UV_image ??? bad FOV



plt.figure()
plt.plot(F3_star_pattern_pos[:,0], -F3_star_pattern_pos[:,1], '+')
plt.xlim([0,1400])
plt.ylim([0,-1100])
for i in  range(4):
    plt.text(F3_star_pattern_pos[i,0], -F3_star_pattern_pos[i,1], str(i))

plt.plot(star_UV_pos[:,0], -star_UV_pos[:,1], 'x')
for i in  range(3):
    plt.text(star_UV_pos[i,0], -star_UV_pos[i,1], str(i))


plt.plot(G2UV.FOV_center_guider_pix[0], -G2UV.FOV_center_guider_pix[1],'o')
plt.grid()

 
    
# read new  pattern
##########################################
path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"

new_F3_star_pattern_filename = path_SCGUI02  + 'Pattern/done/F3_stack2963676.fits_table.fits'

new_F3_star_pattern_tab = Table.read(new_F3_star_pattern_filename)  
new_F3_star_pattern_tab.sort('xcentroid')
new_F3_star_pattern_pos = np.array([new_F3_star_pattern_tab['xcentroid'], new_F3_star_pattern_tab['ycentroid']]).T

slit_pos = np.array([2.7178, -6.0325]) #38

### register with new data
#star_UV_filename = path + 'starsinguiderwhenuvsource/stack1399407_table.fits' #F?
star_UV_filename = path_SCGUI02 + '/F3_stack2998035.fits_table.fits' #F3
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T


G2UV.set_FOV_center(new_F3_star_pattern_pos, star_UV_pos, slit_pos, UVStar_id=3, Stars_id=[0,2], world=False)
#FOV angular position in guider <SkyCoord (SkyOffsetICRS: rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
#    ( 250.42402654,  36.43499699)>): (lon, lat) in deg
#    ( 0.18348537, -0.00658703)>
#FOV pixel position in guider [array(1361.697758698235), array(515.5679367651493)]


#G2UV.save(path + 'Guider2UV_F3.new.pkl')

G2UV= Guider2UV(filename=path + 'Guider2UV_F3.new.pkl')



# predictions
pred = G2UV.pattern_in_slit(starid, slit_pos, world=False)
print(pred)

slit_pos1 = np.array([-2.7599, -2.93116]) #18
starid = 3
pred1 = G2UV.pattern_in_slit(starid, slit_pos1, world=False)
print(pred1)
#[array([  384.09060501,   714.68203387,   778.23701319,  1095.78279597]), 
#array([  840.978344  ,   448.74864922,  1157.86611749,   772.21272886])]

slit_pos2 = np.array([ 4.3915166, -1.8161]) #32
pred2 = G2UV.pattern_in_slit(starid, slit_pos2, world=False)
print(pred2)
#[array([  495.32468487,   826.86997206,   877.5643129 ,  1196.62037055]), 
#array([ 146.46010697, -228.44870725,  477.32957977,  107.47521061])]

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







this_pred = pred2
plt.figure()
plt.plot(new_F3_star_pattern_pos[:,0], -new_F3_star_pattern_pos[:,1], '+')
#plt.gca().invert_xaxis()
plt.axis('equal')
plt.xlim([0,1400])
plt.ylim([-1100, 0])
for i in  range(4):
    plt.text(new_F3_star_pattern_pos[i,0], -new_F3_star_pattern_pos[i,1], str(i))
plt.plot(G2UV.FOV_center_guider_pix[0], -G2UV.FOV_center_guider_pix[1],'o')
plt.grid()
plt.plot(this_pred[0], -this_pred[1],'*r')
for i in  range(4):
    plt.text(this_pred[0][i], -this_pred[1][i], str(i))



all_pred = np.zeros((4,2,len(target_tab)))
for target in range(len(target_tab)):
    this_slit_pos = np.array([target_tab['col4'][target], target_tab['col5'][target]])
    pred = G2UV.pattern_in_slit(starid, this_slit_pos, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
    all_pred[:, 0, target] = pred[0]
    all_pred[:, 1, target] = pred[1]
 





