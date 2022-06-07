#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 15:09:11 2017

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
G2UV= Guider2UV(filename=path + 'Guider2UV_F4.new.pkl')

path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"

detpix = 12e-3 #mm
mask_plate_scale = 42.26134 #mm/deg


# not needed, 4star pattern is in G2UV.pattern
#F4_star_pattern_filename = path_SCGUI02  + 'Pattern/done/F4_stack3111879.fits_table.fits'
#F4_star_pattern_tab = Table.read(F4_star_pattern_filename)  
#F4_star_pattern_tab.sort('xcentroid')
#F4_star_pattern_pos = np.array([F4_star_pattern_tab['xcentroid'], F4_star_pattern_tab['ycentroid']]).T

#########  SCI_GUI02 analysis
###########################################
# slit 1
# img155-167 
# best 163 TU: 2017-09-10T06:10:34
# stack3123494.fits
# vertical offset +7 top

slit_pos1 = np.array([-3.556, -6.04266])  #24
slit_coord1 = G2UV.FieldP.pix2local(slit_pos1)

star_UV_filename = path_SCGUI02 + 'Guider/img/stack3123494.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 2
#xoffset = np.array([-7*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1,3], 
                            UVStar_id=starid)

#slit_obs1 = G2UV.guider2ScienceMask(UVstar_coord, angle=True)  + xoffset

slit_obs1 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs1) 

# slit 2
# img171-183
# best 180 TU: 2017-09-10T06:24:46
# stack3146404.fits
# vertical offset +7 pixel (top)

slit_pos2 = np.array([4.85394, -3.81254]) # 40
slit_coord2 = G2UV.FieldP.pix2local(slit_pos2)

star_UV_filename = path_SCGUI02 + 'Guider/img/stack3146404.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

xoffset = np.array([-7*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,2], 
                            UVStar_id=starid)
#slit_obs2 = G2UV.guider2ScienceMask(UVstar_coord, angle=True) + xoffset
slit_obs2 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs2) 


# slit 3
# img191-203
# best 199 TU: 2017-09-10T06:38:07
# stack3167883.fits

# vertical offset -5 pix (bottom)
slit_pos3 = np.array([1.36652, 2.832]) # OVI 9
slit_coord3 = G2UV.FieldP.pix2local(slit_pos3)


star_UV_filename = path_SCGUI02 + 'Guider/img/stack3167883.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
xoffset = np.array([+5*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0], 
                            UVStar_id=starid)
#slit_obs3 = G2UV.guider2ScienceMask(UVstar_coord, angle=True) + xoffset
slit_obs3 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs3)



slit_pos = np.array([slit_pos1, slit_pos2, slit_pos3])
slit_coord = np.array([[slit_coord1.lon.deg, slit_coord1.lat.deg],
                    [slit_coord2.lon.deg, slit_coord2.lat.deg],
                    [slit_coord3.lon.deg, slit_coord3.lat.deg]])
slit_obs = np.array([[slit_obs1.lon.deg, slit_obs1.lat.deg],
                      [slit_obs2.lon.deg, slit_obs2.lat.deg],
                      [slit_obs3.lon.deg, slit_obs3.lat.deg]])
delta = slit_obs - slit_coord
print(delta * 42.26134)

target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection/targets_F4.txt'
targets = Table.read(target_filename, format='ascii')

tx = np.array(targets['xmm']) + np.array([[-.12],[ +.12]])
ty = np.array(targets['ymm'])[np.newaxis,:]
tx, ty = np.broadcast_arrays(tx, ty)

#plt.figure()
#plt.plot(tx, ty,'-b',linewidth=3.0)
#plt.axis('equal')
#plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
#qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], delta[:,0], delta[:,1])
#plt.quiverkey(qv, .9,.9,.0750, "75 mu",color='r')
#plt.xlim([-13, 13])
#plt.ylim([-7,7])
#plt.xlabel('x mm')
#plt.ylabel('y mm')
#plt.title('Slit measured  displacement on Science Mask F4')

plt.figure()
plt.plot(tx, ty,'-b',linewidth=3.0)
plt.axis('equal')
plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], 0., delta[:,0]*mask_plate_scale, scale=1)
plt.quiverkey(qv, .8, .9, 0.1, "100 micron", color='r')
plt.xlim([-13, 13]) 
plt.ylim([-7,7])
plt.xlabel('x mm')
plt.ylabel('y mm')
plt.title('Slit measured vertical displacement on Science Mask F4')



## rotation, ymagnification & yoffset
#####################################

#mat = np.hstack((slit_coord, np.ones((3,1))))
#matinv =  np.linalg.inv(mat)
#sol = matinv.dot(slit_obs[:,0])
#deltax = sol[2]
#gama = sol[0]
#theta_rad = -sol[1]/gama
#theta = theta_rad*180/np.pi*60 #arcmin
#print(gama, theta, deltax*3600)

#xn -x = x dgama - theta y + dx
mat = np.hstack((slit_coord*[1.,-1], np.ones((3,1))))
matinv =  np.linalg.inv(mat)
sol = matinv.dot(slit_obs[:,0] - slit_coord[:,0])
gama = 1 + sol[0]
theta_rad = sol[1]
deltax = sol[2]
theta = theta_rad*180/np.pi*60 #arcmin
print(gama, theta, deltax*3600)
# 1.00122583875 30.728281425 3.67307744168
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[1, 180/np.pi*60, 3600]) #
# [  1.90772112e-03   7.16481624e+00   7.37224670e-01]
#residual
xnew = gama*(np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1]) + deltax
print((slit_obs[:,0] - xnew)*3600)
# [-0.0173179  -0.01748111  0.00838863]
                   
plt.figure()
plt.plot(tx, ty,'-b',linewidth=3.0)
plt.axis('equal')
plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], 0., delta[:,0]*mask_plate_scale)
plt.quiverkey(qv, .8,.9,.1, "100 mu", color='r')
plt.xlim([-13, 13])
plt.ylim([-7,7])
plt.xlabel('x mm')
plt.ylabel('y mm')
plt.title('Slit measured vertical displacement on Science Mask F4')
plt.text(-10,7,"rotation: {:.2f} arcmin\nmagnification {:.4f}\ndeltay: {:.4f} arcsec".format(theta, gama, deltax*3600))


##### check delta alone: 
#deltax = np.mean(slit_obs[:,0] - slit_coord[:,0])
#print((slit_obs[:,0] - slit_coord[:,0] - deltax)*3600)
#theta = 0.
#gama = 1.

## check without gama
# xn - x = - theta y + dx
mat = np.hstack((-slit_coord[:,[1]], np.ones((3,1))))
matinv = np.linalg.pinv(mat)
sol = matinv.dot(slit_obs[:,0] - slit_coord[:,0])
theta_rad = sol[0]
deltax = sol[1]
theta = theta_rad*180/np.pi*60 #arcmin
print(theta, deltax*3600)
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 3600]) #
#residual
xnew = np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1] + deltax
print((slit_obs[:,0] - xnew)*3600)


# Change G2UV object according to these parameters
###################################################
#G2UVnew  = Guider2UV(guider_wcs=G2UV.GuiderP.w, 
#                     Field_center=G2UV.FieldP.center, 
#                     Field_rotation=G2UV.FieldP.rotation,
#                     mask_rotation=theta*u.arcmin,
#                     Field_gamma=gama)

G2UVnew  = Guider2UV(guider_wcs=G2UV.GuiderP.w, 
                     Field_center=G2UV.FieldP.center, 
                     Field_rotation=G2UV.FieldP.rotation,
                     mask_rotation=theta*u.arcmin)


#deltaFOV = - 1/gama*deltax

#G2UVnew  = Guider2UV(guider_wcs=G2UV.GuiderP.w, 
#                     Field_center=G2UV.FieldP.center, 
#                     Field_rotation=G2UV.FieldP.rotation,
#                     Field_gamma=gama)
#                     #mask_rotation=theta*u.arcmin)
#
#deltaFOV = - deltax

newFOVcenter = coordinates.SkyCoord(G2UV.FOV_center_guider_coord.lon + deltax*u.deg, 
                              G2UV.FOV_center_guider_coord.lat,
                              frame = G2UV.GuiderP.localframe)
                     
G2UVnew.FOV_center_guider_coord = newFOVcenter

#G2UVnew.FOV_guider_coord = G2UV.FOV_guider_coord


# check slit 1
G2UVnew.pattern_coord = G2UV.pattern_coord 


star_UV_filename = path_SCGUI02 + 'Guider/img/stack3123494.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 2

UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1,3], 
                            UVStar_id=starid)


print(G2UVnew.SienceMask2guider(slit_pos1, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos1))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))




star_UV_filename = path_SCGUI02 + 'Guider/img/stack3146404.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0,2], 
                            UVStar_id=starid)

print(G2UVnew.SienceMask2guider(slit_pos2, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos2))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))




star_UV_filename = path_SCGUI02 + 'Guider/img/stack3167883.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0], 
                            UVStar_id=starid)

print(G2UVnew.SienceMask2guider(slit_pos3, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos3))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))


######### save G2UVnew
#G2UVnew.save(path_SCGUI02 + 'Guider2UV_F4.new.pkl')

#
G2UVnew.save(path_SCGUI02 + 'Guider2UV_F4_nogamma.new.pkl')

G2UVnew= Guider2UV(filename=path_SCGUI02 + 'Guider2UV_F4.new.pkl')



