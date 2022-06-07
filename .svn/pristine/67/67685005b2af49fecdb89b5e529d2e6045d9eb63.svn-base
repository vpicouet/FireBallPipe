#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 14:48:20 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy import coordinates
from matplotlib import pyplot as plt
from astropy.table import Table
from astropy import units as u


from guider2UV import Guider2UV

path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/'
path = '/Users/Vincent/Nextcloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/'

G2UV= Guider2UV(filename=path + 'Guider2UV_F2.new.pkl')

detpix = 12e-3 #mm
mask_plate_scale = 42.26134 #mm/deg

path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"
path_SCGUI02 = "/Users/Vincent/Nextcloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"

#########  SCI_GUI02 analysis
###########################################
# slit 1
# img119-131
#best 125 TU: 2017-09-10T05:35:48
# stack3067356.fits (TUguider = TU det - 1'14'')
# vertical  centered
slit_pos1 = np.array([0.32258, -4.0767])  # 30
slit_coord1 = G2UV.FieldP.pix2local(slit_pos1)

star_UV_filename = path_SCGUI02 + '/Guider/img/stack3067356.fits_table.fits' #F2
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

#xoffset = np.array([0.*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1,2], 
                            UVStar_id=starid)

slit_obs1 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs1)


# slit 2
#img137-149
# best 145 TU:2017-09-10T05:51:49
# stack3093204.fits
# vertical centered

slit_pos2 = np.array([1.7145, -1.06426]) #32
slit_coord2 = G2UV.FieldP.pix2local(slit_pos2)

star_UV_filename = path_SCGUI02 + '/Guider/img/stack3093204.fits_table.fits' #F2
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

#xoffset = np.array([0.*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0], 
                            UVStar_id=starid)

slit_obs2 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs2)



# gather all
slit_pos = np.array([slit_pos1, slit_pos2])
slit_coord = np.array([[slit_coord1.lon.deg, slit_coord1.lat.deg],
                    [slit_coord2.lon.deg, slit_coord2.lat.deg]])
slit_obs = np.array([[slit_obs1.lon.deg, slit_obs1.lat.deg],
                      [slit_obs2.lon.deg, slit_obs2.lat.deg]])
delta = slit_obs - slit_coord
print(delta * mask_plate_scale)

target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection/targets_F2.txt'
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
#plt.title('Slit measured  displacement on Science Mask F2')

plt.figure()
plt.plot(tx, ty,'-b',linewidth=3.0)
plt.axis('equal')
plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], 0., delta[:,0]*mask_plate_scale, scale=1)
plt.quiverkey(qv, .8,.9,.1, "50 micron",color='r')
plt.xlim([-13, 13])
plt.ylim([-7,7])
plt.xlabel('x mm')
plt.ylabel('y mm')
plt.title('Slit measured vertical displacement on Science Mask F2')



### determine rotation & magnification correction ? 
# (ref @ first slit pos -> rot center = 1st slit pos)
# use only vertical shift (dispersion dir)

obs = slit_obs - slit_obs[0]

# magnification alone
gama = (obs[1:,0]*slit_pos[1:,0]).sum()/(slit_pos[1:,0]*slit_pos[1:,0]).sum()
print((1-gama)/gama)
# 22.6389692088

# only 2 slit rotation alone:
# xn - x = - theta y + dx
mat = np.hstack((-slit_coord[:,[1]], np.ones((2,1))))
matinv = np.linalg.inv(mat)
sol = matinv.dot(slit_obs[:,0] - slit_coord[:,0])
theta_rad = sol[0]
deltax = sol[1]
theta = theta_rad*180/np.pi*60 #arcmin
print(theta, deltax*3600)
# 127.789085675 -0.940822987283
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 3600]) #
# [ 41.00217279   1.25342776]
#residual
xnew = np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1] + deltax
print((slit_obs[:,0] - xnew)*3600)
print((slit_obs[:,0] - xnew)*mask_plate_scale)
#[-0.23972451 -0.06138147]
#[-0.00281419 -0.00072057]


plt.figure()
plt.plot(tx, ty,'-b',linewidth=3.0)
plt.axis('equal')
plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], 0., delta[:,0]*mask_plate_scale, scale=1)
plt.quiverkey(qv, .9,.9,.0750, "75 mu", color='r')
plt.xlim([-13, 13])
plt.ylim([-7,7])
plt.xlabel('x mm')
plt.ylabel('y mm')
plt.title('Slit measured vertical displacement on Science Mask F2')
plt.text(-10,8,"rotation ALONE: {:.2f} arcmin\nOR magnification ALONE {:.4f}".format(theta,gama))

# Change G2UV object according to these parameters
###################################################
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

newFOV = coordinates.SkyCoord(G2UV.FOV_center_guider_coord.lon + deltax*u.deg, 
                              G2UV.FOV_center_guider_coord.lat,
                              frame = G2UV.GuiderP.localframe)
                     
G2UVnew.FOV_center_guider_coord = newFOV

#G2UVnew.FOV_guider_coord = G2UV.FOV_guider_coord


# check slit 1
G2UVnew.pattern_coord = G2UV.pattern_coord 

#
#G2UVnew.save(path_SCGUI02 + 'Guider2UV_F2_nogamma.new.pkl')

G2UVnew= Guider2UV(filename=path_SCGUI02 + 'Guider2UV_F2_nogamma.new.pkl')



#


############################################################
# with SCI-GUI03 data
############################################################
path_SCGUI03 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170923_SC_GUI03/"


new_F2_star_pattern_SCGUI03_filename = path_SCGUI03  + 'Pattern/F2.csv'

new_F2_star_pattern_SCGUI03_tab = Table.read(new_F2_star_pattern_SCGUI03_filename)  
new_F2_star_pattern_SCGUI03_tab.sort('xcentroid')
new_F2_star_pattern_SCGUI03_pos = np.array([new_F2_star_pattern_SCGUI03_tab['xcentroid'], new_F2_star_pattern_SCGUI03_tab['ycentroid']]).T

target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection/targets_F2.txt'
F2 = Table.read(target_filename, format='ascii')

# slit #16
slit_pos3 =   np.array([F2[F2['Internal-count']=='16']['xmm'][0], F2[F2['Internal-count']=='16']['ymm'][0]])
slit_coord3 = G2UV.FieldP.pix2local(slit_pos3)

# slit-width scan
star_UV_filename = path_SCGUI03 + 'Guider/img/stack2391542_table.csv'
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 2
UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[1], UVStar_id=starid,
                                        FourStarsGuider_pos=new_F2_star_pattern_SCGUI03_pos)
UVstar_coord_lon = UVstar_coord.lon

# slit-height scan
star_UV_filename = path_SCGUI03 + 'Guider/img/stack2386795_table.csv'
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 2
UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[1], UVStar_id=starid,
                                        FourStarsGuider_pos=new_F2_star_pattern_SCGUI03_pos)
# reset the slit width position the above centered one
UVstar_coord = coordinates.SkyCoord(UVstar_coord_lon, UVstar_coord.lat,
                                    frame=G2UV.GuiderP.localframe)
slit_obs3 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs3) 

# slit #30 again
slit_pos4 = slit_pos1
slit_coord4 = slit_coord1
star_UV_filename = path_SCGUI03 + 'Guider/img/stack2429702_table.csv'
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 2
UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[1], UVStar_id=starid,
                                        FourStarsGuider_pos=new_F2_star_pattern_SCGUI03_pos)

slit_obs4 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs4) 


