#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 15:05:28 2017

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
G2UV= Guider2UV(filename=path + 'Guider2UV_F3.new.pkl')


detpix = 12e-3 #mm
mask_plate_scale = 42.26134 #mm/deg

path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"




#########  SCI_GUI02 analysis
###########################################
# slit 1
# img69-81
# best  76  TU 2017-09-10T04:45:04
# stack2985505.fits (TUguider = TUdet + 1'14'')
# vertical +2 (top)

slit_pos1 = np.array([2.7178, -6.0325]) #29
slit_coord1 = G2UV.FieldP.pix2local(slit_pos1)

star_UV_filename = path_SCGUI02 + 'Guider/img/stack2985505.fits_table.fits' #F3
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
#xoffset = np.array([-2*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,2], 
                            UVStar_id=starid)

slit_obs1 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs1) 


# slit 2
# img86-98
#  best 89 TU:2017-09-10T05:09:47
# stack3025471.fits
# vertical +6 pixel
slit_pos2 = np.array([-2.7599, -2.93116]) # 18
slit_coord2 = G2UV.FieldP.pix2local(slit_pos2)

star_UV_filename = path_SCGUI02 + 'Guider/img/stack3025471.fits_table.fits' #F3
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
#xoffset = np.array([-6*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1], 
                            UVStar_id=starid)

slit_obs2 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs2) 

# slit 3
# img104-116
# best 110   TU: 2017-09-10T05:22:57
# stack3046667.fits
# vertical -1 (bottom)
slit_pos3 = np.array([ 4.3915166, -1.8161])  # 32
slit_coord3 = G2UV.FieldP.pix2local(slit_pos3)

star_UV_filename = path_SCGUI02 + 'Guider/img/stack3046667.fits_table.fits' #F3
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

#xoffset = np.array([+1*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,2], 
                            UVStar_id=starid)

slit_obs3 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs3)

# gather all
slit_pos = np.array([slit_pos1, slit_pos2, slit_pos3])
slit_coord = np.array([[slit_coord1.lon.deg, slit_coord1.lat.deg],
                    [slit_coord2.lon.deg, slit_coord2.lat.deg],
                    [slit_coord3.lon.deg, slit_coord3.lat.deg]])
slit_obs = np.array([[slit_obs1.lon.deg, slit_obs1.lat.deg],
                      [slit_obs2.lon.deg, slit_obs2.lat.deg],
                      [slit_obs3.lon.deg, slit_obs3.lat.deg]])
delta = slit_obs - slit_coord
print(delta * 42.26134)


target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection/targets_F3.txt'
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
#plt.title('Slit measured  displacement on Science Mask F3')

plt.figure()
plt.plot(tx, ty,'-b',linewidth=3.0)
plt.axis('equal')
plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], 0., delta[:,0]*mask_plate_scale, scale=1)
plt.quiverkey(qv, .8, .9, .1, "100 micron", color='r')
plt.xlim([-13, 13])
plt.ylim([-7,7])
plt.xlabel('x mm')
plt.ylabel('y mm')
plt.title('Slit measured vertical displacement on Science Mask F3')


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
# 0.99160754831 43.5895402355 -7.46050683622
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[1, 180/np.pi*60, 3600]) #
# [  3.79925744e-03   7.63255215e+00   1.31780968e+00]
#residual
xnew = gama*(np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1]) + deltax
print((slit_obs[:,0] - xnew)*3600)
# [-0.01625817 -0.04501122  0.02761195]
                   
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
plt.title('Slit measured vertical displacement on Science Mask F3')
plt.text(-10,7,"rotation: {:.2f} arcmin\nmagnification {:.4f}\ndeltay: {:.4f} mm".format(theta, gama, deltax*mask_plate_scale))



## check without gama
# xn - x = - theta y + dx
mat = np.hstack((-slit_coord[:,[1]], np.ones((3,1))))
matinv = np.linalg.pinv(mat)
sol = matinv.dot(slit_obs[:,0] - slit_coord[:,0])
theta_rad = sol[0]
deltax = sol[1]
theta = theta_rad*180/np.pi*60 #arcmin
print(theta, deltax*3600)
# 44.2677914054 -4.91478287366
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 3600]) #
# [ 7.62637373  0.6391298 ]
#residual
xnew = np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1] + deltax
print((slit_obs[:,0] - xnew)*3600)
# [ 1.68153248 -0.42410034 -1.33346867]

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

newFOV_center = coordinates.SkyCoord(G2UV.FOV_center_guider_coord.lon + deltax*u.deg, 
                              G2UV.FOV_center_guider_coord.lat,
                              frame = G2UV.GuiderP.localframe)
                     
G2UVnew.FOV_center_guider_coord = newFOV_center

#G2UVnew.FOV_guider_coord = G2UV.FOV_guider_coord


# check slit 1
G2UVnew.pattern_coord = G2UV.pattern_coord 

star_UV_filename = path_SCGUI02 + 'Guider/img/stack2985505.fits_table.fits' #F3
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0,2], 
                            UVStar_id=starid)


print(G2UVnew.SienceMask2guider(slit_pos1, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos1))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))



star_UV_filename = path_SCGUI02 + 'Guider/img/stack3025471.fits_table.fits' #F3
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1], 
                            UVStar_id=starid)

print(G2UVnew.SienceMask2guider(slit_pos2, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos2))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))



star_UV_filename = path_SCGUI02 + 'Guider/img/stack3046667.fits_table.fits' #F3
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0,2], 
                            UVStar_id=starid)

print(G2UVnew.SienceMask2guider(slit_pos3, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos3))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))


######### save G2UVnew
#G2UVnew.save(path_SCGUI02 + 'Guider2UV_F3.new.pkl')

#G2UVnew.save(path_SCGUI02 + 'Guider2UV_F3_nogamma.new.pkl')

#
G2UVnew= Guider2UV(filename=path_SCGUI02 + 'Guider2UV_F3.new.pkl')











