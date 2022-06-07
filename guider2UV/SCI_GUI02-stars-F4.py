#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 15:08:17 2017

@author: dvibert
"""


from __future__ import division, print_function
from scipy.optimize import curve_fit

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy import units as u
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib import pyplot as plt
from astropy.table import Table
import matplotlib.patches as patches

from guider2UV import Guider2UV


path = '/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/'
G2UV= Guider2UV(filename=path + 'Guider2UV_F4.pkl')


#################################
#
# predict guiding stars
#
##################################
star_target_path = '/Users/vincent/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/GuidingStars/'
#F4_stars = Table.read(star_target_path + "F4_guidingstars.txt", format='ascii')
F4_stars = Table.read(star_target_path + "F4_guidingstars.fits", format='fits')

mag = np.fmin(np.fmin(F4_stars['GAIA gband'].filled(99), F4_stars['SDSS gband'].filled(99)),
        F4_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(F4_stars['RA']*u.deg, F4_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)

#new positions
path_SCGUI02 = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI02/'
G2UVnew= Guider2UV(filename=path_SCGUI02 + 'Guider2UV_F4_nogamma.pkl')

guider_star_pos_new = G2UVnew.SienceMask2guider(coords, world=True, angle=False)

plt.figure()
plt.plot(guider_star_pos[0], guider_star_pos[1], '*')
#plt.xlim([0,1400])
#plt.ylim([0,-1100])

#write table
F4_stars['Xguider'] = guider_star_pos[0]
F4_stars['Yguider'] = guider_star_pos[1]


skytest_path = '/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/FlightConfigurationSKyTest/img/'
seen_stars = Table.read(skytest_path + "stack0537972_table.csv", format='csv')

F4_stars.write(star_target_path + "F4_guidingstars.fits", format='fits', overwrite=True)
F4_stars.write(star_target_path + "F4_guidingstars.txt", format='ascii.fixed_width', overwrite=True)

skytest_path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/FlightConfigurationSKyTest/img/'
seen_stars = Table.read(skytest_path + "stack0537972_table.csv")

###
select = mag <=13

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(F4_stars['Xguider'][select], F4_stars['Yguider'][select], '+r', label='pred stars')
plt.plot(F4_stars['Xguider2'][select], F4_stars['Yguider2'][select], 'xb', label='pred stars, rotation')

plt.plot(seen_stars['xcentroid'], seen_stars['ycentroid'],'*', label='observed stars')

plt.ylim([1200,-200]) #reverse y axis
plt.xlabel('Xguider in pixels')
plt.ylabel('Yguider in pixels')
plt.plot(G2UV.FOV_guider_pix[0], G2UV.FOV_guider_pix[1], 'og')
plt.text(G2UV.FOV_guider_pix[0], G2UV.FOV_guider_pix[1], 'FOV')

ax.add_patch(
        patches.Rectangle(
        (0, 0),
        1280,
        1080,
        fill=False      # remove background
    )
   )
for s in  np.nonzero(select)[0]:
    plt.text(F4_stars['Xguider'][s], F4_stars['Yguider'][s], str(F4_stars['Internal count'][s]))

plt.legend()


######  test 3 ***** 16/09/2017

pattern_pos = np.array([[496.5,323.7,984.6,812.0],[147.4,630.5,314.6,792.2]]).T

plt.figure()
plt.plot(pattern_pos[:,0],pattern_pos[:,1],'+')
plt.ylim([1200,-100])
#for i in range(4):
#    plt.text(pattern_pos[:,0],pattern_pos[:,1], str(i))
    
guid_star_pos_ang = coordinates.SkyCoord(36.987137*u.deg, 0.402799*u.deg) #star 29
star_id = 0
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#[array([ 375.6064398 ,  204.75996257,  869.4419162 ,  698.7445199 ]), 
#array([ 274.0377231 ,  756.39561629,  440.83902922,  917.62890504])]
# guiding on star2:  1 pix error
 


########################### VINCENT'S PART ##############################


select = mag <=12
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')
plt.plot(F4_stars['Xguider'], F4_stars['Yguider'], '+r', label='predicted mag>12')
plt.scatter(F4_stars['Xguider'][select], F4_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F4_stars['Xguider'][s], F4_stars['Yguider'][s], str(F4_stars['Internal count'][s]))

#plt.plot(F1_stars['Xguider'][select], F1_stars['Yguider'][select], '+r')
plt.xlim([-200, 1400])
plt.ylim([1400, -200])
plt.xlabel('Xguider in pixels')
plt.ylabel('Yguider in pixels')
plt.plot(G2UV.FOV_guider_pix[0], G2UV.FOV_guider_pix[1], 'og')
plt.text(G2UV.FOV_guider_pix[0], G2UV.FOV_guider_pix[1], 'FOV')
plt.grid()
ax.add_patch(
        patches.Rectangle(
        (0, 0),
        1280,
        1080,
        fill=False))
plt.legend()


select = mag <=14

plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 400, marker='*', color='m', label='observed')
plt.scatter(F4_stars['Xguider'][select], F4_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F4_stars['Xguider'][s], F4_stars['Yguider'][s], str(F4_stars['Internal count'][s]))
    
    
counts = [14,18,29,34]#F1_stars['Internal count']#
                 
for i in counts:
    select = F4_stars['Internal count'] == i
#    print (F1_stars[select])
    plt.scatter(F4_stars['Xguider'][select], F4_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')  
    print (str(F4_stars['Internal count'][select][0]))
    plt.text(F4_stars['Xguider'][select], F4_stars['Yguider'][select], str(F4_stars['Internal count'][select]))
#plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')
c=3
plt.scatter(seen_stars['xcentroid'][c], seen_stars['ycentroid'][c], 100, marker='*', color='m', label='observed')


def rot_offset((x,y),dx,dy,theta):
    x1 = np.cos(theta) * x - np.sin(theta) * y + dx
    y1 = np.sin(theta) * x + np.cos(theta) * y + dy
    return np.array([x1,y1]).flatten()

#seen_stars.remove_raw(-1)
#seen_stars.remove_raw(3)

counts = [34,29,18,14]#F1_stars['Internal count']#                 
x=[]
y=[]
x1=[]
y1=[]
for i in range(len(counts)):
    select = F4_stars['Internal count'] == counts[i]
    x.append(F4_stars['Xguider'][select][0])
    y.append(F4_stars['Yguider'][select][0])
    x1.append(seen_stars['xcentroid'][i])
    y1.append(seen_stars['ycentroid'][i])

x=np.array(x)
y=np.array(y)
x1=np.array(x1)
y1=np.array(y1)

xy = np.vstack((x,y)).T
x1y1 = np.vstack((x1,y1)).T
xy_angle = G2UV.GuiderP.pix2local(xy)
x1y1_angle = G2UV.GuiderP.pix2local(x1y1)

xy_sec = np.array([xy_angle.lon.arcsec,xy_angle.lat.arcsec]).T
x1y1_sec = np.array([x1y1_angle.lon.arcsec,x1y1_angle.lat.arcsec]).T

                   
                   
a = curve_fit(rot_offset,(xy_sec[:,0],xy_sec[:,1]),np.array([x1y1_sec[:,0],x1y1_sec[:,1]]).flatten())[0]
print('The offset are ({},{}) arcsec, and the rotation is {} arc minutes'.format(a[0],a[1],a[2]*180*60/np.pi ))

xy_pix = coordinates.SkyCoord(xy_sec[:,0]*u.arcsec,xy_sec[:,1]*u.arcsec,frame = G2UV.GuiderP.localframe)
plt.figure()
plt.scatter(x1y1_sec[:,0],x1y1_sec[:,1], 100, marker='*', color='b', label='observed')
plt.scatter(xy_sec[:,0],xy_sec[:,1], 100, marker='*', color='m', label='theoric')
plt.xlim((-60*20,60*+20))
plt.ylim((60*20,-60*20))
plt.scatter(rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,10)[0], rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,10)[1], 100, marker='*', color='r', label='Theoric after offset & rot')
plt.legend()

q1 = plt.quiver(x1y1_sec[:,0],x1y1_sec[:,1],x1y1_sec[:,0]-rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,4)[0],x1y1_sec[:,1]-rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,4)[1])
plt.quiverkey(q1,0.5,0.5,10,color='r', label = '10 arcsec')
for i in range(10):
    plt.text(x1y1_sec[i,0],x1y1_sec[i,1], int(np.sqrt(np.square(x1y1_sec[i,0]-rot_offset((xy_sec[i,0],xy_sec[i,1]),*a)[0]) + np.square(x1y1_sec[i,1]-rot_offset((xy_sec[i,0],xy_sec[i,1]),*a)[1])) ))
plt.xlim((-60*20,60*+20))
plt.ylim((60*20,-60*20))
plt.grid()                   
                  