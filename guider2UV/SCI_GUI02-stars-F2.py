#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 14:45:47 2017

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

path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/'
path = '/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/'

G2UV= Guider2UV(filename=path + 'Guider2UV_F2.pkl')

#################################
#
# predict guiding stars
#
##################################
star_target_path = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/GuidingStars/'
star_target_path = '/Users/vincent/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/GuidingStars/'
F2_stars = Table.read(star_target_path + "F2_guidingstars.fits", format='fits')
mag = np.fmin(np.fmin(F2_stars['GAIA gband'].filled(99), F2_stars['SDSS gband'].filled(99)),
        F2_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(F2_stars['RA']*u.deg, F2_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)

#F2_stars.write(star_target_path + "F2_guidingstars.fits", format='fits', overwrite=True)
#F2_stars.write(star_target_path + "F2_guidingstars.txt", format='ascii.fixed_width', overwrite=True)

skytest_path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/FlightCinfigurationSKyTest/img/'
skytest_path = '/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/FlightConfigurationSKyTest/img/'
seen_stars = Table.read(skytest_path + "stack0297812.fits_table.fits", format='fits')

skytest_path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/FlightConfigurationSKyTest/img/'
seen_stars = Table.read(skytest_path + "stack0297812_table.csv")

plt.figure()
plt.plot(guider_star_pos[0], guider_star_pos[1], '*')
#plt.xlim([0,1400])
#plt.ylim([0,-1100])

select = mag <=12

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(F2_stars['Xguider'][select], F2_stars['Yguider'][select], '+r', label='pred stars')
plt.plot(seen_stars['xcentroid'], seen_stars['ycentroid'],'*', label='observed stars')

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
for s in  np.nonzero(select)[0]:
    plt.text(F2_stars['Xguider'][s], F2_stars['Yguider'][s], str(F2_stars['Internal count'][s]))

plt.legend()


############# test 3 ***** 16/09/2017
pattern_pos = np.array([[774.5,99.5,272.2,597.8],[326.3,644.5,157.5,808.7]]).T


guid_star_pos = np.array([ 624.1,   655.7]) #star 21
guid_star_pos_ang = coordinates.SkyCoord(252.829951*u.deg, 34.94576*u.deg)
star_id = 3
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#guide étoile 1 --- 14h35 ?
#array([ 805.41583971,  127.63935786,  305.42707496,  624.04280487]), 
#array([ 170.8262446 ,  488.11394605,   -1.28319044,  655.62983843])]
# starid @ ~ 625,654    
    
guid_star_pos_ang = coordinates.SkyCoord(252.895501*u.deg,  34.875333*u.deg) #starv #36
star_id = 3
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#guide étoile 1 ---- 14h42 -> guiding lost (too big move with high gain)
# netttoyage fenêtre ---- 14h56, star @ pos
#[array([ 923.27142733,  242.67667975,  431.34293942,  731.84331909]), 
#array([-178.81360396,  137.49965532, -356.74371233,  310.93065954])]
# star_id @ ~734.5 311.8

guid_star_pos_ang = coordinates.SkyCoord(252.83754*u.deg, 34.836825*u.deg) #star 40
star_id = 2
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#[array([ 986.98366113,  328.23107319,  496.23595395,  814.42774473]), 
# array([ 386.62602009,  710.07894789,  226.44855113,  865.97530951])]
#guide étoile 1 (place étoile 2) - 15h09 guiding
# starid @ ~503 - 228


guid_star_pos_ang = coordinates.SkyCoord(252.787662*u.deg, 34.810634*u.deg ) #star 46
star_id = 2
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#[array([ 799.10256065,  126.39966503,  298.12915269,  623.26969091]), 
#arra([ 345.09711451,  664.05465368,  177.65119515,  826.94050532])]
# starid @ ~ 302 177  guiding @15h20

############################################
pattern_pos_distcor = G2UV.GuiderP.w.sip_pix2foc(pattern_pos,1)
delta = guid_star_pos - G2UV.GuiderP.w.wcs.crpix - pattern_pos_distcor
new_pos = pattern_pos_distcor - delta[star_id]
new_pos = G2UV.GuiderP.w.foc2pix(new_pos,1) + G2UV.GuiderP.w.wcs.crpix



                                
########################### VINCENT'S PART ##############################


select = mag <=12
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')
plt.plot(F2_stars['Xguider'], F2_stars['Yguider'], '+r', label='predicted mag>12')
plt.scatter(F2_stars['Xguider'][select], F2_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F2_stars['Xguider'][s], F2_stars['Yguider'][s], str(F2_stars['Internal count'][s]))
plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'],
            color='g',marker='D', label='observed shifted')#, alpha=.7)
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
plt.scatter(F2_stars['Xguider'][select], F2_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F2_stars['Xguider'][s], F2_stars['Yguider'][s], str(F2_stars['Internal count'][s]))
    
    
counts = [46,42,36,49,38,18,20,27,21,16]#F1_stars['Internal count']#
                 
for i in counts:
    select = F2_stars['Internal count'] == i
#    print (F1_stars[select])
    plt.scatter(F2_stars['Xguider'][select], F2_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')  
    print (str(F2_stars['Internal count'][select][0]))
    plt.text(F2_stars['Xguider'][select], F2_stars['Yguider'][select], str(F2_stars['Internal count'][select][0]))
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
    select = F2_stars['Internal count'] == counts[i]
    x.append(F2_stars['Xguider'][select][0])
    y.append(F2_stars['Yguider'][select][0])
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
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                 
a1 = curve_fit(rot_offset,(x,y),np.array([x1,y1]).flatten())[0]
print('The offset are ({},{}) pixels, and the rotation is {} arc minutes'.format(a[0],a[1],a[2]*180*60/np.pi ))

plt.figure()
plt.scatter(x, y, 100, marker='*', color='m', label='theoric')
plt.scatter(x1, y1, 100, marker='*', color='b', label='observed')
plt.xlim((0,1280))
plt.ylim((1080,0))
plt.scatter(rot_offset((x,y),*a1).reshape(2,5)[0], rot_offset((x,y),*a1).reshape(2,5)[1], 100, marker='*', color='r', label='Theoric after offset & rot')
plt.legend()

plt.figure()
q1 = plt.quiver(x1,y1,x1-rot_offset((x,y),*a1).reshape(2,5)[0],y1-rot_offset((x,y),*a1).reshape(2,5)[1])
plt.quiverkey(q1,0.5,0.5,10,color='r', label = '10 pixel')
plt.xlim((0,1280))
plt.ylim((1080,0))
plt.grid()
                                
                                
                                