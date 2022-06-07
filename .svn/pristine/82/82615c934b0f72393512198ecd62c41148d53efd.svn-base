#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 18:24:38 2017

@author: vincent
"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy import units as u
from astropy.wcs import WCS
from scipy.optimize import curve_fit
from astropy.coordinates import Angle

path = '/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/perf3.0-09-01/sky_00_astrometry/'
wF3 = WCS(fits.open('stack475788.new')[0]) 
wF2 = WCS(fits.open('stack486607.new')[0])  #-161
wF1 = WCS(fits.open('stack500513.new')[0])
wF4 = WCS(fits.open('stack502282.new')[0])#159 
F3 = fits.open('stack475788.new')[0] 
F2 = fits.open('stack486607.new')[0] 
F1 = fits.open('stack500513.new')[0]
F4 = fits.open('stack502282.new')[0]


wF1b =  WCS('stack500633.new')  #119
wF1c = WCS('stack500731.new') 
wF4b = WCS('stack502179.new') #19
wF4c = WCS('stack502352.new')
wF4d = WCS('stack528075.new')




stars_onF2 = Table.read(path + 'guider/stack1293657_table.fits')['xcentroid','ycentroid']            
stars_onF1 = Table.read(path + 'guider/stack1373757_table.fits')['xcentroid','ycentroid']            


#wF2.sip_pix2foc(stars_onF2,1)
x2,y2 = wF2.sip_pix2foc(stars_onF2['xcentroid'],stars_onF2['ycentroid'],1)
x11,y11 = wF1.sip_pix2foc(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)
x12,y12 = wF1b.sip_pix2foc(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)
x13,y13 = wF1c.sip_pix2foc(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)


plt.figure()
qv1 = plt.quiver(x11, y11, x12-x11, y12-y11, scale_units='xy', angles='xy', color='c', label = 'Distortion 1-2',scale=0.01 )                                         
plt.quiver(x11, y11, x13-x11, y13-y11, scale_units='xy', angles='xy', color='b', label = 'Distortion 1-2',scale=0.01 )                                         
plt.quiverkey(qv1, .2,-0.16, 1, 'scale: 13 microns', coordinates='axes', color='r')                                                                  
plt.legend()
plt.xlim((-640,540))
plt.ylim((-540,540))
#plt.figtext(0.13,0.13,'mean = %.2f - %.2f - %.2f \nmax = %.2f - %.2f - %.2f'%(np.mean(y12-y13),np.mean(x12-x13),np.mean(dist),np.var(y-valy),np.var(x-valx),np.var(dist),np.max(y-valy),np.max(x-valx),np.max(dist)))
plt.xlabel('y')
plt.ylabel('x')
plt.grid(10)
plt.show()


x11,y11 = wF1.all_pix2world(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)
x12,y12 = wF1b.all_pix2world(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)
x13,y13 = wF1c.all_pix2world(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)

plt.figure()
qv1 = plt.quiver(x11, y11, x12-x11, y12-y11, scale_units='xy', angles='xy', color='c', label = 'Distortion 1-2',scale=0.02 )                                         
plt.quiver(x11, y11, x13-x11, y13-y11, scale_units='xy', angles='xy', color='b', label = 'Distortion 1-2',scale=0.02 )                                         
plt.quiverkey(qv1, .2,-0.16, 5/60/60, 'scale: 5 sec', coordinates='axes', color='r')                                                                  
plt.legend()
#plt.figtext(0.13,0.13,'mean = %.2f - %.2f - %.2f \nmax = %.2f - %.2f - %.2f'%(np.mean(y12-y13),np.mean(x12-x13),np.mean(dist),np.var(y-valy),np.var(x-valx),np.var(dist),np.max(y-valy),np.max(x-valx),np.max(dist)))
plt.xlabel('y')
plt.ylabel('x')
plt.grid(10)
plt.show()
       


       
wcsF3 = WCS('stack486607.new')
w = WCS('stack486607.new')
w = WCS('stack500513.new')
w = WCS('stack502282.new')

def sky_Coord((ra,dec),ra_c,dec_c):
    """Given a line sight, the function converts right ascenssion and declination
    into alpha_x alpha_y the small angles to center line of sight (ra_c,dec_c)"""
    numx = -np.cos(dec*np.pi/180)*np.sin((ra - ra_c)*np.pi/180)
    numy = -np.cos((ra - ra_c)*np.pi/180) * np.cos(dec*np.pi/180) * np.sin(dec_c*np.pi/180)  + np.sin(dec*np.pi/180) * np.cos(dec_c*np.pi/180)
    denom = np.cos((ra - ra_c)*np.pi/180) * np.cos(dec*np.pi/180) * np.cos(dec_c*np.pi/180)  + np.sin(dec*np.pi/180) * np.sin(dec_c*np.pi/180)
    alpha_x = - np.arctan(numx/denom)*180/np.pi
    alpha_y = np.arctan(numy/denom)*180/np.pi
    return np.array([alpha_x,alpha_y]).ravel()

def fonction(x,y,ra_center,dec_center):
    coords = SkyCoord(x*u.deg, y*u.deg)
    F1_center = coordinates.SkyCoord(ra_center*u.deg, dec_center*u.deg)
    F1_Frame = coordinates.SkyOffsetFrame(origin=F1_center)    
    return coords.transform_to(F1_Frame)





#################TEST RECUPERATION RA,DEC###########################
rac = -2
decc = 1
x11,y11 = sky_Coord((a11,d11),rac,decc).reshape(2,4)[0,:],sky_Coord((a11,d11),rac,decc).reshape(2,4)[1,:]
ra_center, dec_center = curve_fit(sky_Coord,[a11,d11],np.array([x11,y11]).flat)[0]
#ra_center, dec_center = curve_fit(sky_Coord,[x11,y11],np.array([a11,d11]).flat)[0]
print (ra_center, dec_center)

#################TEST RECUPERATION RA,DEC###########################
#rac = 12.5
#decc = 30
a11 -= np.mean(a11)
d11 -= np.mean(d11)
xx11,yy11 = sky_Coord1((a11,d11),rac,decc).reshape(2,4)[0,:],sky_Coord1((a11,d11),rac,decc).reshape(2,4)[1,:]
scatter(a11,d11);scatter(xx11,yy11)
ra_center, dec_center = curve_fit(sky_Coord1,[a11,d11],np.array([xx11,yy11]).flat)[0]
#ra_center, dec_center = curve_fit(sky_Coord,[x11,y11],np.array([a11,d11]).flat)[0]
print (ra_center, dec_center)



#4 etoiles sur le plan focal
x11,y11 = wF1.sip_pix2foc(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)
#4 etoiles sur le ciel
a11,d11 = wF1.all_pix2world(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)
#4 etoiles sur le plan focal
ra_center, dec_center = curve_fit(sky_Coord,[a11,d11],np.array([x11,y11]).flat)[0]
scatter(x11,y11)

#################TROUVER LA POSITION DE LA 4eme ETOILE DANS LE GUIDER###########################


stars_onF4x= Table.read(path + 'guider/stack1399407_table.fits')['xcentroid']           
stars_onF4y= Table.read(path + 'guider/stack1399407_table.fits')['ycentroid']
#stars_onF4x, stars_onF4y = wF4.sip_pix2foc(stars_onF4x,stars_onF4x,1)
stars_onF4a, stars_onF4d = wF4.all_pix2world(stars_onF4x,stars_onF4y,1)
ra_center, dec_center = curve_fit(sky_Coord,[stars_onF4a,stars_onF4d],np.array([stars_onF4x,stars_onF4y]).flat)[0]
















def sky_Coord1((ra,dec),ra_c,dec_c):
    """Given a line sight, the function converts right ascenssion and declination
    into alpha_x alpha_y the small angles to center line of sight (ra_c,dec_c)"""
    a = SkyCoord(ra, dec, "icrs", unit="deg")
    F1_center = coordinates.SkyCoord(ra_c*u.deg, dec_c*u.deg)
    F1_Frame = coordinates.SkyOffsetFrame(origin=F1_center)
    a2 = a.transform_to(F1_Frame)
    return np.array([a2.lon.deg,a2.lat.deg]).ravel()

def Find_4thStar_on_UV_Visible_Frame(Stars4_on_patern_frame, Stars3_on_UV_Visibe_frame):
    ra_center, dec_center = curve_fit(sky_Coord1,[Stars4_on_patern_frame.lon.deg[:-1],Stars4_on_patern_frame.lat.deg[:-1]],np.array([Stars3_on_UV_Visibe_frame.lon.deg,Stars3_on_UV_Visibe_frame.lat.deg]).flat)[0]
    Stars4_on_UV_Visibe_framex, Stars4_on_UV_Visibe_framey = sky_Coord1((Stars4_on_patern_frame.lon.deg,Stars4_on_patern_frame.lat.deg),ra_center,dec_center).reshape(2,4)[0,:],sky_Coord1((Stars4_on_patern_frame.lon.deg,Stars4_on_patern_frame.lat.deg),ra_center,dec_center).reshape(2,4)[1,:]
    return [Stars4_on_UV_Visibe_framex[-1],Stars4_on_UV_Visibe_framey[-1]]
    
    
stars_onF4x= Table.read(path + 'guider/stack1399407_table.fits')['xcentroid']           
stars_onF4y= Table.read(path + 'guider/stack1399407_table.fits')['ycentroid']
stars_onF4a, stars_onF4d = wF4.all_pix2world(stars_onF4x,stars_onF4y,1)
a4 = SkyCoord(stars_onF4a, stars_onF4d, "icrs", unit="deg")

stars_onF1 = Table.read(path + 'guider/stack1373757_table.fits')['xcentroid','ycentroid']            
a11,d11 = wF1.all_pix2world(stars_onF1['xcentroid'],stars_onF1['ycentroid'],1)
a3 = SkyCoord(a11,d11, "icrs", unit="deg")


patern_Frame = coordinates.SkyOffsetFrame(origin=a4[0])
UV_Visible_Frame = coordinates.SkyOffsetFrame(origin=a4[0])
a4 = a4.transform_to(patern_Frame)
a3 = a3.transform_to(UV_Visible_Frame)


Find_4thStar_on_UV_Visible_Frame(a3,a4)
scatter(a4.lon.deg,a4.lat.deg);scatter(Find_4thStar_on_UV_Visible_Frame(a3,a4)[0],Find_4thStar_on_UV_Visible_Frame(a3,a4)[1])


i=1
scatter(a4[:i].lon.deg,a4[:i].lat.deg);scatter(a3[:i].lon.deg,a3[:i].lat.deg)

sep1 = a4[1].separation(a4[2])
sep2 = a4[0].separation(a4[1])
sep3 = a4[0].separation(a4[2])
print (sep1, sep2, sep3)
sep1 = a3[1].separation(a3[2])
sep2 = a3[0].separation(a3[1])
sep3 = a3[0].separation(a3[2])
print (sep1, sep2, sep3)
