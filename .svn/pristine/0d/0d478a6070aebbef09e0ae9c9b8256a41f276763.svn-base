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


path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/'
F4_astrometry_filename = path + 'Distortion/F4_stack502282.new'
Field_center = coordinates.SkyCoord(36.9049*u.deg, 0.65245*u.deg)
G2UV = Guider2UV(F4_astrometry_filename, Field_center )
#G2UVang = Guider2UV(F4_astrometry_filename, Field_center  )



### register with new data

path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"

F4_star_pattern_filename = path_SCGUI02  + 'Pattern/done/F4_stack3111879.fits_table.fits'

F4_star_pattern_tab = Table.read(F4_star_pattern_filename)  
F4_star_pattern_tab.sort('xcentroid')
F4_star_pattern_pos = np.array([F4_star_pattern_tab['xcentroid'], F4_star_pattern_tab['ycentroid']]).T

slit_pos = np.array([-3.556, -6.04266]) #15

star_UV_filename = path_SCGUI02 + 'F4_stack3131287.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T


G2UV.set_FOV_center(F4_star_pattern_pos, star_UV_pos, slit_pos, UVStar_id=2, Stars_id=[0,1,3], world=False)
#FOV angular position in guider <SkyCoord (SkyOffsetICRS: rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
#    ( 250.39242605,  36.41307253)>): (lon, lat) in deg
#    ( 0.18491328,  0.01856568)>
#FOV pixel position in guider [array(1373.624249970326), array(605.2129609883267)]

#
G2UV.save(path + 'Guider2UV_F4.new.pkl')

G2UV= Guider2UV(filename=path + 'Guider2UV_F4.new.pkl')


# predictions
starid = 2
pred = G2UV.pattern_in_slit(starid, slit_pos, world=False)
print(pred)
#pred was: 
#[array([  422.81368264,   756.34545346,   814.58498277,  1134.38432638]), 
#array([ 632.81407971,  242.74162261,  954.08454469,  570.29015194])]  

slit_pos1 = np.array([4.85394, -3.81254]) #31
starid = 3
pred1 = G2UV.pattern_in_slit(starid, slit_pos1, world=False)
print(pred1)
#[array([  306.13009931,   652.69360359,   692.675947  ,  1024.63859803]), 
#array([ 195.00421859, -197.84825757,  531.28746886,  144.45233252])]

slit_pos2 = np.array([1.36652, 2.832]) #0vi9
starid = 3
pred2 = G2UV.pattern_in_slit(starid, slit_pos2, world=False)
print(pred2)
#[array([  944.41426349,  1266.9315126 ,  1319.73178159,  1629.86032581]), 
#array([ 540.17822684,  150.62622214,  852.32649324,  469.26267328])]

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
plt.plot(F4_star_pattern_pos[:,0], -F4_star_pattern_pos[:,1], '+')
#plt.gca().invert_xaxis()
plt.axis('equal')
plt.xlim([0,1400])
plt.ylim([-1100, 0])
for i in  range(4):
    plt.text(F4_star_pattern_pos[i,0], -F4_star_pattern_pos[i,1], str(i))
plt.plot(G2UV.FOV_center_guider_pix[0], -G2UV.FOV_center_guider_pix[1],'o')
plt.grid()
plt.plot(this_pred[0], -this_pred[1],'*r')
for i in  range(4):
    plt.text(this_pred[0][i], -this_pred[1][i], str(i))






