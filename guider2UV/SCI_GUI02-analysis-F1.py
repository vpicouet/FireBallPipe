#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 14:43:40 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy import coordinates
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.table import Table
from guider2UV import Guider2UV


path_SCGUI01 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/"
path_SCGUI01 = "/Users/Vincent/Nextcloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/"
G2UV= Guider2UV(filename=path_SCGUI01  + 'Guider2UV_F1.new.pkl')

detpix = 12e-3 #mm
mask_plate_scale = 42.26134 #mm/deg

path_SCGUI02 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"
path_SCGUI02 = "/Users/Vincent/Nextcloud/FIREBALL/Tests-at-FortSumner/170909_SC_GUI02/"


#########  SCI_GUI02 analysis
###########################################
# slit 1
#img dither 12-20
#best img16 2017-09-10T03:09:08
#stack2830609.fits (TUguider = TUdet - 1'14'')
# vertical decentering +7 (top of slit)

slit_pos1 = np.array([0.5525186, -4.7601449]) # slit 30
slit_coord1 = G2UV.FieldP.pix2local(slit_pos1)

star_UV_filename = path_SCGUI02 + 'Guider/img/stack2830609.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
#xoffset = np.array([-7*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1,2], 
                            UVStar_id=starid)

slit_obs1 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs1) 
#<SkyCoord (SkyOffsetICRS: rotation=-90.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
#    (32.19, -5.688)>): (lon, lat) in deg
#    (-0.11232722, -0.01240308)>

#slit 2
#img dither 35 47 
#best 45  TU: 2017-09-10T03:32:15 
# stack2868425.fits
# vertical centering +1  (top)

slit_pos2 = np.array([4.376, 2.644]) # 39
slit_coord2 = G2UV.FieldP.pix2local(slit_pos2)

star_UV_filename = path_SCGUI02 + 'Guider/img/stack2868425.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
#xoffset = np.array([-1*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0], 
                            UVStar_id=starid)

slit_obs2 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs2) 
#<SkyCoord (SkyOffsetICRS: rotation=-90.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
#    (32.19, -5.688)>): (lon, lat) in deg
#    (0.06613708, -0.1027078)>

# slit 3
#img dither img52-64
#best 55 TU:2017-09-10T04:05:00
#stack2920845.fits
# vertical +5 pix
slit_pos3 = np.array([-3.655,-.894]) #18
slit_coord3 = G2UV.FieldP.pix2local(slit_pos3)

star_UV_filename = path_SCGUI02 + 'Guider/img/stack2920845.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3
#xoffset = np.array([-5*detpix, 0])

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1], 
                            UVStar_id=starid)

slit_obs3 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs3)
#<SkyCoord (SkyOffsetICRS: rotation=-90.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
#    (32.19, -5.688)>): (lon, lat) in deg
#    (-0.02228355, 0.08692987)>
    
# gather all
slit_pos = np.array([slit_pos1, slit_pos2, slit_pos3])
slit_coord = np.array([[slit_coord1.lon.deg, slit_coord1.lat.deg],
                    [slit_coord2.lon.deg, slit_coord2.lat.deg],
                    [slit_coord3.lon.deg, slit_coord3.lat.deg]])
slit_obs = np.array([[slit_obs1.lon.deg, slit_obs1.lat.deg],
                      [slit_obs2.lon.deg, slit_obs2.lat.deg],
                      [slit_obs3.lon.deg, slit_obs3.lat.deg]])
delta = slit_obs - slit_coord
print(delta * mask_plate_scale)


target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection/targets_F1.txt'
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
#plt.title('Slit measured  displacement on Science Mask F1')

plt.figure()
plt.plot(tx, ty,'-b',linewidth=3.0)
plt.axis('equal')
plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], 0., delta[:,0]*mask_plate_scale, scale=1)
plt.quiverkey(qv, .8,.9,.1, "100 micron",color='r')
plt.xlim([-13, 13])
plt.ylim([-7,7])
plt.xlabel('x mm')
plt.ylabel('y mm')
plt.title('Slit measured vertical displacement on Science Mask F1')


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
# 1.00720682874 73.6143058008 3.16786600186
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[1, 180/np.pi*60, 3600]) #
# [  2.51377154e-03   7.96717929e+00   6.34744906e-01]
#residual
xnew = gama*(np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1]) + deltax
print((slit_obs[:,0] - xnew)*3600)
# [-0.1008538  -0.00490627  0.0299566 ]
                   
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
plt.title('Slit measured vertical displacement on Science Mask F1')
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
# 83.9681507113 2.44270893215
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 3600]) #
# [ 7.10163893  0.58217143]
#residual
xnew = np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1] + deltax
print((slit_obs[:,0] - xnew)*3600)
# [-2.4608359   1.29394115  1.09065295]

# Change G2UV object according to these parameters
###################################################
G2UVnew  = Guider2UV(guider_wcs=G2UV.GuiderP.w, 
                     Field_center=G2UV.FieldP.center, 
                     Field_rotation=G2UV.FieldP.rotation,
                     mask_rotation=theta*u.arcmin,
                     Field_gamma=gama)

#G2UVnew  = Guider2UV(guider_wcs=G2UV.GuiderP.w, 
#                     Field_center=G2UV.FieldP.center, 
#                     Field_rotation=G2UV.FieldP.rotation,
#                     mask_rotation=theta*u.arcmin)

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

#G2UVnew.FOV_center_guider_coord = G2UV.FOV_center_guider_coord


# check slit 1
G2UVnew.pattern_coord = G2UV.pattern_coord 

star_UV_filename = path_SCGUI02 + 'Guider/img/stack2830609.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1,2], 
                            UVStar_id=starid)


print(G2UVnew.SienceMask2guider(slit_pos1, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos1))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))


# check slit 2
star_UV_filename = path_SCGUI02 + 'Guider/img/stack2868425.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0], 
                            UVStar_id=starid)

print(G2UVnew.SienceMask2guider(slit_pos2, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos2))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))


# check slit 3
star_UV_filename = path_SCGUI02 + 'Guider/img/stack2920845.fits_table.fits' #F4
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 3

UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1], 
                            UVStar_id=starid)

print(G2UVnew.SienceMask2guider(slit_pos3, angle=True).to_string(u'dms')) 
print(UVstar_coord.to_string(u'dms'))

print(G2UVnew.FieldP.pix2local(slit_pos3))
print(G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True))


######### save G2UVnew
#G2UVnew.save(path_SCGUI02 + 'Guider2UV_F1.pkl')
#G2UVnew.save(path_SCGUI02 + 'Guider2UV_F1.new.pkl')

#G2UVnew.save(path_SCGUI02 + 'Guider2UV_F1_nogamma.pkl')
#G2UVnew.save(path_SCGUI02 + 'Guider2UV_F1_nogamma.new.pkl')

G2UVnew= Guider2UV(filename=path_SCGUI02 + 'Guider2UV_F1.new.pkl')

#G2UVnew= Guider2UV(filename=path_SCGUI02 + 'Guider2UV_F1_nogamma.new.pkl')


############################################################
# with SCI-GUI03 data
############################################################
path_SCGUI03 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170923_SC_GUI03/"


new_F1_star_pattern_SCGUI03_filename = path_SCGUI03  + 'Pattern/F1.csv'

new_F1_star_pattern_SCGUI03_tab = Table.read(new_F1_star_pattern_SCGUI03_filename)  
new_F1_star_pattern_SCGUI03_tab.sort('xcentroid')
new_F1_star_pattern_SCGUI03_pos = np.array([new_F1_star_pattern_SCGUI03_tab['xcentroid'], new_F1_star_pattern_SCGUI03_tab['ycentroid']]).T

# slit 1
slit_pos4 = slit_pos1
slit_coord4 = slit_coord1

star_UV_filename = path_SCGUI03 + 'Guider/img/stack2292839_table.csv'
#star_UV_filename = path_SCGUI03 + 'Guider/img/stack2293514_table.csv'
# when arriving (code prediction)
#star_UV_filename = path_SCGUI03 + 'Guider/img/stack2292839_table.csv'

star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 2
#xoffset = np.array([-7*detpix, 0])
#UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[0,3], 
#                            UVStar_id=starid)

UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos[[0]], Stars_id=[0], UVStar_id=starid,
                                        FourStarsGuider_pos=new_F1_star_pattern_SCGUI03_pos)

slit_obs4 = G2UV.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs4) 

print(slit_obs4.lon.arcsec - slit_coord4.lon.arcsec, \
      slit_obs4.lat.arcsec - slit_coord4.lat.arcsec)



#UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos, Stars_id=[0,1,3], UVStar_id=starid,
#                                        FourStarsGuider_pos=new_F1_star_pattern_SCGUI03_pos)

UVstar_coord = G2UVnew.UVstar_guider_coord(star_UV_pos[[0]], Stars_id=[0], UVStar_id=starid,
                                        FourStarsGuider_pos=new_F1_star_pattern_SCGUI03_pos)

slit_obs4 = G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs4) 

print(slit_obs4.lon.arcsec - slit_coord4.lon.arcsec, \
      slit_obs4.lat.arcsec - slit_coord4.lat.arcsec)


# slit 11
target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/targets_F1.txt'
F1 = Table.read(target_filename, format='ascii')
slit_pos5 =   np.array([F1[F1['Internal-count']=='11']['xmm'][0], F1[F1['Internal-count']=='11']['ymm'][0]])
slit_coord5 = G2UV.FieldP.pix2local(slit_pos5)
# slit width scan
star_UV_filename = path_SCGUI03 + 'Guider/img/stack2317935_table.csv'
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 2
UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[1], UVStar_id=starid,
                                        FourStarsGuider_pos=new_F1_star_pattern_SCGUI03_pos)
UVstar_coord_lon = UVstar_coord.lon

# slit height scan
star_UV_filename = path_SCGUI03 + 'Guider/img/stack2335491_table.csv'
star_UV_tab = Table.read(star_UV_filename)  
star_UV_tab.sort('xcentroid')
star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
starid = 2
UVstar_coord = G2UV.UVstar_guider_coord(star_UV_pos, Stars_id=[1], UVStar_id=starid,
                                        FourStarsGuider_pos=new_F1_star_pattern_SCGUI03_pos)
# reset the slit width position the above centered one
UVstar_coord = coordinates.SkyCoord(UVstar_coord_lon, UVstar_coord.lat,
                                    frame=G2UV.GuiderP.localframe)
slit_obs5 = G2UVnew.guider_to_FieldLocal(UVstar_coord, angle=True)
print(slit_obs5) 

print(slit_obs5.lon.arcsec - slit_coord5.lon.arcsec, \
      slit_obs5.lat.arcsec - slit_coord5.lat.arcsec)

    
# reduce all to roatation, magnification, offset
slit_pos = np.array([slit_pos1, slit_pos2, slit_pos3, slit_pos4, slit_pos5])
slit_coord = np.array([[slit_coord1.lon.deg, slit_coord1.lat.deg],
                    [slit_coord2.lon.deg, slit_coord2.lat.deg],
                    [slit_coord3.lon.deg, slit_coord3.lat.deg],
                    [slit_coord4.lon.deg, slit_coord4.lat.deg],
                    [slit_coord5.lon.deg, slit_coord5.lat.deg]
                    ])
slit_obs = np.array([[slit_obs1.lon.deg, slit_obs1.lat.deg],
                      [slit_obs2.lon.deg, slit_obs2.lat.deg],
                      [slit_obs3.lon.deg, slit_obs3.lat.deg],
                      [slit_obs4.lon.deg, slit_obs4.lat.deg],
                      [slit_obs5.lon.deg, slit_obs5.lat.deg],
                      ])
delta = slit_obs - slit_coord
print(delta * mask_plate_scale)
#[[ 0.01470639  0.02854039]
# [ 0.14999359  0.03716569]
# [-0.04753929  0.01798185]
# [-0.0394451   0.12005495]
# [ 0.05992337  0.01464285]]
# 

# slit-width:  xn - x = x dgamma - theta y   + dx  
# slit-height: yn - y = x theta  + y dgamma  + dy
row0   = np.hstack((slit_coord[4,::-1], np.array([0., 1.]))) #slit height
row1_5 = np.hstack((slit_coord*[1.,-1], np.ones((5,1)), np.zeros((5,1)))) # slit width
mat = np.vstack((row0, row1_5))
matinv =  np.linalg.pinv(mat)
data = np.concatenate((slit_obs[[4],1] - slit_coord[[4],1], slit_obs[:,0] - slit_coord[:,0]))
sol = matinv.dot(data)
gama = 1 + sol[0]
theta_rad = sol[1]
deltax = sol[2]
deltay = sol[3]
theta = theta_rad*180/np.pi*60 #arcmin
print(gama, theta, deltax*3600, deltay*3600)
# 1.01869144515 18.4416455612 5.45642003738 -11.331117112
# accuracy, assuming 1 arcsec measurement error
covar =np.linalg.inv(mat.T.dot(mat))
print (np.sqrt(np.diag(covar))/3600*[1, 180/np.pi*60, 3600, 3600]) #
# [  1.84821931e-03   4.27736586e+00   5.27968895e-01   1.60496310e+00]
#residual
xnew = gama*(np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1]) + deltax
print((slit_obs[:,0] - xnew)*3600)
# [ 3.11492133  1.07461917 -6.38159152 -1.49793212  3.73241649]
ynew = gama*(np.sin(theta_rad)*slit_coord[4,0] + np.cos(theta_rad)*slit_coord[4,1]) + deltay
print((slit_obs[4,1] - ynew)*3600)
# 0.0119627939654
print((slit_obs[:,0] - xnew)*mask_plate_scale)
print((slit_obs[4,1] - ynew)*mask_plate_scale)
#[ 0.03656687  0.01261524 -0.07491517 -0.01758462  0.04381581]
#0.000140434361979

plt.figure()
plt.plot(tx, ty,'-b',linewidth=3.0)
plt.axis('equal')
plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], 0., delta[:,0]*mask_plate_scale, scale=.5)
plt.quiverkey(qv, .8,.9,.05, "50 mu", color='r')
plt.xlim([-13, 13])
plt.ylim([-7,7])
plt.xlabel('x mm')
plt.ylabel('y mm')
plt.title('Slit measured vertical displacement on Science Mask F1')
plt.text(-10,7,'''
rotation: {:.2f} arcmin
magnification {:.4f}
deltay: {:.4f} mm
deltax: {:.4f} mm'''.format(theta, gama, deltax*mask_plate_scale, deltay*mask_plate_scale))


# fixed gama=1
# slit-width:  xn - x =  - theta y   + dx  
# slit-height: yn - y = x theta  + dy
row0   = np.hstack((slit_coord[4,0], np.array([0., 1.]))) #slit height
row1_5 = np.hstack((-slit_coord[:,[1]], np.ones((5,1)), np.zeros((5,1)))) # slit width
mat = np.vstack((row0, row1_5))
matinv =  np.linalg.pinv(mat)
data = np.concatenate((slit_obs[[4],1] - slit_coord[[4],1], slit_obs[:,0] - slit_coord[:,0]))
sol = matinv.dot(data)
theta_rad = sol[0]
gama = 1.
deltax = sol[1]
deltay = sol[2]
theta = theta_rad*180/np.pi*60 #arcmin
print(gama, theta, deltax*3600, deltay*3600)
# accuracy, assuming 1 arcsec measurement error
covar =np.linalg.inv(mat.T.dot(mat))
print (np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 3600, 3600]) #
#residual
xnew = gama*(np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1]) + deltax
print((slit_obs[:,0] - xnew)*3600)
ynew = gama*(np.sin(theta_rad)*slit_coord[4,0] + np.cos(theta_rad)*slit_coord[4,1]) + deltay
print((slit_obs[4,1] - ynew)*3600)

print((slit_obs[:,0] - xnew)*mask_plate_scale)
print((slit_obs[4,1] - ynew)*mask_plate_scale)



