#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 14:50:02 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy import units as u
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.table import Table

from guider2UV import Guider2UV


path = '/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/'
G2UV= Guider2UV(filename=path + 'Guider2UV_F3.pkl')


#################################
#
# predict guiding stars
#
##################################
star_target_path = '/Users/vincent/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/GuidingStars/'
#F3_stars = Table.read(star_target_path + "F3_guidingstars.txt", format='ascii')
F3_stars = Table.read(star_target_path + "F3_guidingstars.fits", format='fits')
mag = np.fmin(np.fmin(F3_stars['GAIA gband'].filled(99), F3_stars['SDSS gband'].filled(99)),
        F3_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(F3_stars['RA']*u.deg, F3_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)

skytest_path = '/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/FlightConfigurationSKyTest/img/'

seen_stars = Table.read(skytest_path + "stack0363324_table.csv", format='csv')

coords = coordinates.SkyCoord(F3_stars['RA']*u.deg, F3_stars['DEC']*u.deg)
guider_star_pos_new = G2UVnew.SienceMask2guider(coords, world=True, angle=False)


skytest_path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/FlightConfigurationSKyTest/img/'

seen_stars = Table.read(skytest_path + "stack0363324_table.csv")

# as seen on guider image
plt.figure()
plt.plot(F3_stars['DEC'], F3_stars['RA'],'*')
plt.plot(Field_center.dec.deg, Field_center.ra.deg, 'og')
plt.text(Field_center.dec.deg, Field_center.ra.deg, 'FOV')
plt.xlabel('Dec')
plt.ylabel('Ra')

# as seen on guider image
plt.figure()
plt.plot(guider_star_pos[0], guider_star_pos[1], '*')
plt.ylim(plt.ylim()[::-1]) #reverse y axis

#write table
F3_stars['Xguider'] = guider_star_pos[0]
F3_stars['Yguider'] = guider_star_pos[1]

#F3_stars.write(star_target_path + "F3_guidingstars.fits", format='fits', overwrite=True)
#F3_stars.write(star_target_path + "F3_guidingstars.txt", format='ascii.fixed_width', overwrite=True)

#F3_stars.write(star_target_path + "F3_guidingstars.fits", format='fits', overwrite=True)
#F3_stars.write(star_target_path + "F3_guidingstars.txt", format='ascii.fixed_width', overwrite=True)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(F3_stars['Xguider'], F3_stars['Yguider'], '+r')
plt.plot(F3_stars['Xguider2'], F3_stars['Yguider2'], 'xb')

plt.ylim(plt.ylim()[::-1]) #reverse y axis
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


select = mag <=18

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(F3_stars['Xguider'][select], F3_stars['Yguider'][select], '+r', label='pred stars')
plt.plot(F3_stars['Xguider'][select], F3_stars['Yguider'][select], '+r', label='pred stars')

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
    plt.text(F3_stars['Xguider'][s], F3_stars['Yguider'][s], str(F3_stars['Internal count'][s]))

plt.legend()



###### test 3 ***** 16/09/2017
#check guider distortion with these stars.

pattern_pos = np.array([[767.0, 104.2,276.1,600.4],[419.5,740.0,252.4,901.0]]).T

plt.figure()
plt.plot(pattern_pos[:,0],pattern_pos[:,1],'+')
for i in range(4):
    plt.text(pattern_pos[:,0],pattern_pos[:,1], str(i))
    
guid_star_pos_ang = coordinates.SkyCoord(352.228576*u.deg, -0.041124*u.deg) #5
star_id = 3
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#[array([ 523.94532891, -154.60502913,   21.86176671,  353.41231888]), 
#array([ 493.0368871 ,  822.85071116,  326.07031836,  982.97274023])]
# starid @ 356 - 980   @ 15h43


 
guid_star_pos_ang = coordinates.SkyCoord(352.267202 *u.deg,  0.107397*u.deg) # star 9
star_id = 3
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#[array([ 1116.71583691,   475.34516764,   640.74863905,   955.78654743]), 
#array([ 343.02870576,  652.16300192,  176.7783518 ,  813.22707853])]
### guiding on  star 2: starid 957.8 811.4 @16h09

########################### VINCENT'S PART ##############################


select = mag <=12
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')
plt.plot(F3_stars['Xguider'], F3_stars['Yguider'], '+r', label='predicted mag>12')
plt.scatter(F3_stars['Xguider'][select], F3_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F3_stars['Xguider'][s], F3_stars['Yguider'][s], str(F3_stars['Internal count'][s]))

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
plt.scatter(F3_stars['Xguider'][select], F3_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F3_stars['Xguider'][s], F3_stars['Yguider'][s], str(F3_stars['Internal count'][s]))
    
    
counts = [46,42,36,49,38,18,20,27,21,16]#F1_stars['Internal count']#
                 
for i in counts:
    select = F3_stars['Internal count'] == i
#    print (F1_stars[select])
    plt.scatter(F3_stars['Xguider'][select], F3_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')  
    print (str(F3_stars['Internal count'][select][0]))
    plt.text(F3_stars['Xguider'][select], F3_stars['Yguider'][select], str(F3_stars['Internal count'][select][0]))
plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')
c=7
plt.scatter(seen_stars['xcentroid'][c], seen_stars['ycentroid'][c], 100, marker='*', color='m', label='observed')


def rot_offset((x,y),dx,dy,theta):
    x1 = np.cos(theta) * x - np.sin(theta) * y + dx
    y1 = np.sin(theta) * x + np.cos(theta) * y + dy
    return np.array([x1,y1]).flatten()

#seen_stars.remove_raw(-1)
#seen_stars.remove_raw(3)

counts = [49,46,42,38,36,27,21,20,18,16]#F1_stars['Internal count']#                 
x=[]
y=[]
x1=[]
y1=[]
for i in range(len(counts)):
    select = F3_stars['Internal count'] == counts[i]
    x.append(F3_stars['Xguider'][select][0])
    y.append(F3_stars['Yguider'][select][0])
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

q1 = plt.quiver(x1y1_sec[:,0],x1y1_sec[:,1],x1y1_sec[:,0]-rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,10)[0],x1y1_sec[:,1]-rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,10)[1])
plt.quiverkey(q1,0.5,0.5,10,color='r', label = '10 arcsec')
for i in range(10):
    plt.text(x1y1_sec[i,0],x1y1_sec[i,1], int(np.sqrt(np.square(x1y1_sec[i,0]-rot_offset((xy_sec[i,0],xy_sec[i,1]),*a)[0]) + np.square(x1y1_sec[i,1]-rot_offset((xy_sec[i,0],xy_sec[i,1]),*a)[1])) ))
plt.xlim((-60*20,60*+20))
plt.ylim((60*20,-60*20))
plt.grid()                   
                   
                   
                   
                   