#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 14:40:43 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy import units as u
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from astropy.table import Table
from astroML.crossmatch import crossmatch
from guider2UV import Guider2UV
from scipy.optimize import curve_fit


path_SCGUI01 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/"
path_SCGUI01 = "/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/"
G2UV= Guider2UV(filename=path_SCGUI01  + 'Guider2UV_F1.pkl')

#################################
#
# predict guiding stars
#
##################################
star_target_path = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/GuidingStars/'
star_target_path = '/Users/vincent/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/GuidingStars/'
F1_stars = Table.read(star_target_path + "F1_guidingstars.fits", format='fits')
mag = np.fmin(np.fmin(F1_stars['GAIA gband'].filled(99), F1_stars['SDSS gband'].filled(99)),
        F1_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(F1_stars['RA']*u.deg, F1_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)


#new positions
path_SCGUI02 = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI02/'
G2UVnew= Guider2UV(filename=path_SCGUI02 + 'Guider2UV_F1_nogamma.pkl')

guider_star_pos_new = G2UVnew.SienceMask2guider(coords, world=True, angle=False)

# as seen on guider image
#plt.figure()
#plt.plot(F1_stars['DEC'], F1_stars['RA'],'*')
#plt.plot(Field_center.dec.deg, Field_center.ra.deg, 'og')
#plt.text(Field_center.dec.deg, Field_center.ra.deg, 'FOV')
#plt.xlabel('Dec')
#plt.ylabel('Ra')

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(guider_star_pos[0], guider_star_pos[1], '*')
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


#write table
#F1_stars['Xguider'] = guider_star_pos[0]
#F1_stars['Yguider'] = guider_star_pos[1]

#F1_stars.write(star_target_path + "F1_guidingstars.fits", format='fits', overwrite=True)
#F1_stars.write(star_target_path + "F1_guidingstars.txt", format='ascii.fixed_width', overwrite=True)

#F1_stars.write(star_target_path + "F1_guidingstars.fits", format='fits', overwrite=True)
#F1_stars.write(star_target_path + "F1_guidingstars.txt", format='ascii.fixed_width', overwrite=True)

#select brihgt stars
select = mag <=12

# read guide pos from sky test (FightConfiguration)
skytest_path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/FlightConfigurationSKyTest/img/'
skytest_path = '/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/FlightConfigurationSKyTest/img/'
seen_stars = Table.read(skytest_path + "stack0721926.fits_table.fits", format='fits')

# plot 
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(F1_stars['Xguider'][select], F1_stars['Yguider'][select], '+r')
plt.plot(F1_stars['Xguider2'][select], F1_stars['Yguider2'][select], 'xb')

plt.plot(seen_stars['xcentroid'], seen_stars['ycentroid'],'*')
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

# find shift on guiding star
dx = F1_stars['Xguider'][24] - seen_stars['xcentroid'][2]
dy = F1_stars['Yguider'][24] - seen_stars['ycentroid'][2]

select = mag <=12
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')
plt.plot(F1_stars['Xguider'], F1_stars['Yguider'], '+r', label='predicted mag>12')
plt.scatter(F1_stars['Xguider'][select], F1_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F1_stars['Xguider'][s], F1_stars['Yguider'][s], str(F1_stars['Internal count'][s]))

plt.scatter(seen_stars['xcentroid']+dx, seen_stars['ycentroid']+dy,
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
        fill=False      # remove background
    )
   )
plt.legend()


dist, idx = crossmatch(np.array(seen_stars['xcentroid', 'ycentroid']).view((float,2)) + np.array([dx,dy]),
           np.array(F1_stars['Xguider', 'Yguider'][select]).view((float,2)))

#idx = [4,7,24,30,37]
matched = np.nonzero(select)[0][idx]
plt.scatter(F1_stars['Xguider'][matched], F1_stars['Yguider'][matched], 100, 
            marker="D", color='b')

F1_stars[matched]



########  test 3 ***** 16/09/2017
pattern_pos = np.array([[216.1, 390.5, 882.5, 710.0],[635.5, 149.8,318.0, 797.2]]).T

guid_star_pos_ang = coordinates.SkyCoord(32.221491*u.deg, -5.776507*u.deg) # 25
star_id = 2
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#[array([ 332.22978835,  501.07121913,  990.21736915,  822.80549421]), 
#array([ 798.17767841,  314.39163152,  476.91827045,  954.29085556])]
# guiding on star 1: 991 472  



guid_star_pos_ang = coordinates.SkyCoord( 32.329364  *u.deg,  -5.904128*u.deg) # 41
star_id = 1
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)
#[array([ 304.89206244,  480.52324399,  967.9584886 ,  794.25087525]), 
#array([ 517.21886917,   35.23849289,  205.94437842,  681.74042702])]
# guiding on star 2: ok
# but star 41 is to close to top border

guid_star_pos_ang = coordinates.SkyCoord(32.254526*u.deg, -5.815258*u.deg) # 31
star_id = 3
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)


guid_star_pos_ang = coordinates.SkyCoord(32.1014 *u.deg, -5.786504*u.deg) # 8
star_id = 3
pred = G2UV.pattern_in_slit(star_id, guid_star_pos_ang, world=True, FourStarsGuider_pos=pattern_pos)
print(pred)






########################### VINCENT'S PART ##############################



select = mag <=12
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')
plt.plot(F1_stars['Xguider'], F1_stars['Yguider'], '+r', label='predicted mag>12')
plt.scatter(F1_stars['Xguider'][select], F1_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F1_stars['Xguider'][s], F1_stars['Yguider'][s], str(F1_stars['Internal count'][s]))
plt.scatter(seen_stars['xcentroid']+dx, seen_stars['ycentroid']+dy,
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



plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')
plt.scatter(F1_stars['Xguider'][select], F1_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(select)[0]:
    plt.text(F1_stars['Xguider'][s], F1_stars['Yguider'][s], str(F1_stars['Internal count'][s]))
    
    
counts = [38,31,25,8,5]#F1_stars['Internal count']#
                 
for i in counts:
    select = F1_stars['Internal count'] == i
#    print (F1_stars[select])
    plt.scatter(F1_stars['Xguider'][select], F1_stars['Yguider'][select], 100, 
            marker="*", color='r', label='predicted mag<=12')  
    plt.text(F1_stars['Xguider'][select], F1_stars['Yguider'][select], str(F1_stars['Internal count'][select][0]))
plt.scatter(seen_stars['xcentroid'], seen_stars['ycentroid'], 100, marker='*', color='m', label='observed')

c=4
plt.scatter(seen_stars['xcentroid'][c], seen_stars['ycentroid'][c], 100, marker='*', color='m', label='observed')


def rot_offset((x,y),dx,dy,theta):
    x1 = np.cos(theta) * x - np.sin(theta) * y + dx
    y1 = np.sin(theta) * x + np.cos(theta) * y + dy
    return np.array([x1,y1]).flatten()

#rot_offset(1,1,0,0,pi/2)

counts = [38,31,25,8,5]#F1_stars['Internal count']#                 
x=[]
y=[]
x1=[]
y1=[]
for i in range(len(counts)):
    select = F1_stars['Internal count'] == counts[i]
    x.append(F1_stars['Xguider'][select][0])
    y.append(F1_stars['Yguider'][select][0])
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
plt.scatter(rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,5)[0], rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,5)[1], 100, marker='*', color='r', label='Theoric after offset & rot')
plt.legend()

q1 = plt.quiver(x1y1_sec[:,0],x1y1_sec[:,1],x1y1_sec[:,0]-rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,5)[0],x1y1_sec[:,1]-rot_offset((xy_sec[:,0],xy_sec[:,1]),*a).reshape(2,5)[1])
plt.quiverkey(q1,0.5,0.5,10,color='r', label = '10 arcsec')
for i in range(5):
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
