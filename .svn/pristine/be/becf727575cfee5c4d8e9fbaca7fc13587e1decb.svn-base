#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 17:42:52 2017

@author: dvibert
"""


from __future__ import division, print_function

import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.table import Table

from guider2UV import Guider2UV
from matplotlib import pyplot as plt

mask_plate_scale = 42.26134 #mm/deg


# encoder gains:
CEg = 1.027 # ? TBC
Elg = 1.001  # ? TBC


#guider center => put source on guider center
gc = np.array([640, 540])

##############################
# Field 4
##############################
path = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/'
G2UV= Guider2UV(filename=path + 'Guider2UV_F4.pkl')
gc_coord = G2UV.GuiderP.pix2local([gc])

#encoder at center: EL 1149361, CE 1050525
# got to slit 1
###############
slit_pos1 = np.array([-3.556, -6.04266])  #24
slit_coord1 = G2UV.FieldP.pix2local(slit_pos1)
print(slit_coord1.to_string(u'dms'))
delta_x_guider = G2UV.FOV_guider_coord.lon - gc_coord.lon + slit_coord1.lon
delta_y_guider = G2UV.FOV_guider_coord.lat - gc_coord.lat + slit_coord1.lat

theta_EL =  delta_x_guider/2./Elg
theta_CE = -delta_y_guider/2./CEg
print(theta_EL.to('arcsec'), theta_CE.to('arcsec'))
# [u'75.2081arcsec'] [u'-180.125arcsec']
#encoder at slit first: CE 1068535, El 1148543
    
# then scan slit X & Y
# note guider frame delta between prediction and best centered position
delta_EL = -10./3600/2/Elg #   -10 arcsec/2
delta_CE = 10./3600/2/CEg  #   +10 arcesc/2

slit_obs1  = coordinates.SkyCoord(slit_coord1.lon + 2*delta_EL*u.deg,
                                        slit_coord1.lat - 2*delta_CE*u.deg,
                                        frame=G2UV.FieldP.localframe)
print(slit_obs1.to_string(u'dms')) 


# go back to guider center and then to slit 2
#############################################
# EL 1149339, CE 1050325
slit_pos2 = np.array([4.85394, -3.81254]) # 40
slit_coord2 = G2UV.FieldP.pix2local(slit_pos2)
print(slit_coord2.to_string(u'dms'))

delta_x_guider = G2UV.FOV_guider_coord.lon - gc_coord.lon + slit_coord2.lon
delta_y_guider = G2UV.FOV_guider_coord.lat - gc_coord.lat + slit_coord2.lat

theta_EL =  delta_x_guider/2./Elg
theta_CE = -delta_y_guider/2./CEg

print(theta_EL.to('arcsec'), theta_CE.to('arcsec'))
#[u'170.196arcsec'] [u'168.881arcsec']
#encoder: El 1147489 CE 1033435


# then scan slit (X dispersion El 
# note guider frame delta between prediction and best centered position
delta_EL = 5./3600/2/Elg # +5 arcsec/2
delta_CE = 0. # TBC deg

slit_obs2  = coordinates.SkyCoord(slit_coord2.lon + 2*delta_EL*u.deg,
                                  slit_coord2.lat - 2*delta_CE*u.deg,
                                  frame=G2UV.FieldP.localframe)
print(slit_obs2.to_string(u'dms')) 


# go to slit 3
###############
# El 1149339, CE 1050525
slit_pos3 = np.array([1.36652, 2.832]) # OVI 9
slit_coord3 = G2UV.FieldP.pix2local(slit_pos3)
print(slit_coord3.to_string(u'dms'))

delta_x_guider = G2UV.FOV_guider_coord.lon - gc_coord.lon + slit_coord3.lon
delta_y_guider = G2UV.FOV_guider_coord.lat - gc_coord.lat + slit_coord3.lat

theta_EL =  delta_x_guider/2./Elg
theta_CE = -delta_y_guider/2./CEg
print(theta_EL.to('arcsec'), theta_CE.to('arcsec'))
# [u'453.03arcsec'] [u'24.1418arcsec']
# EL 1144414, CE 1048115

# then scan slit X (dispersion, EL)
# note guider frame delta between prediction and best centered position
delta_EL = 0. # 0 
delta_CE = 0. # TBC

slit_obs3  = coordinates.SkyCoord(slit_coord3.lon + 2*delta_EL*u.deg,
                                  slit_coord3.lat - 2*delta_CE*u.deg,
                                  frame=G2UV.FieldP.localframe)

print(slit_obs3.to_string(u'dms')) 

# compute orientation/magnification offset of science mask :
############################################################

slit_pos = np.array([slit_pos1, slit_pos2, slit_pos3])
slit_coord = np.array([[slit_coord1.lon.deg, slit_coord1.lat.deg],
                    [slit_coord2.lon.deg, slit_coord2.lat.deg],
                    [slit_coord3.lon.deg, slit_coord3.lat.deg]])
slit_obs = np.array([[slit_obs1.lon.deg, slit_obs1.lat.deg],
                      [slit_obs2.lon.deg, slit_obs2.lat.deg],
                      [slit_obs3.lon.deg, slit_obs3.lat.deg]])
delta = slit_obs - slit_coord
print(delta * 42.26134)

# without using slit2 height scan 
#xn -x = x dgama - theta y + dx 
mat = np.hstack((slit_coord*[1.,-1], np.ones((3,1))))
matinv =  np.linalg.inv(mat)
sol = matinv.dot(slit_obs[:,0] - slit_coord[:,0])
gama = 1 + sol[0]
theta_rad = sol[1]
deltax = sol[2]
theta = theta_rad*180/np.pi*60 #arcmin
print(gama, theta, deltax*3600)
# 1.00189048815 70.1378223203 -2.83144033712
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[1, 180/np.pi*60, 3600]) #
# [  1.90772112e-03   7.16481624e+00   7.37224670e-01]
#residual
xnew = gama*(np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1]) + deltax
print((slit_obs[:,0] - xnew)*3600)
#[-0.09614687 -0.08312773  0.04598416]


xnew_rc6 = np.array([-0.14299218, -0.08832364,  0.06841098])
print((xnew - xnew_rc6)*3600)



target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/targets_F4.txt'
targets = Table.read(target_filename, format='ascii')

tx = np.array(targets['xmm']) + np.array([[-.12],[ +.12]])
ty = np.array(targets['ymm'])[np.newaxis,:]
tx, ty = np.broadcast_arrays(tx, ty)

plt.figure()
plt.plot(tx, ty,'-b',linewidth=3.0)
plt.axis('equal')
plt.plot(slit_pos[:,0], slit_pos[:,1], 'or')
qv = plt.quiver(slit_pos[:,0], slit_pos[:,1], 0., delta[:,0]*mask_plate_scale)
plt.quiverkey(qv, .8,.9,.075, "75 mu", color='r')
plt.xlim([-13, 13])
plt.ylim([-7,7])
plt.xlabel('x mm')
plt.ylabel('y mm')
plt.title('Slit measured vertical displacement on Science Mask F4')
plt.text(-10,7,"rotation: {:.2f} arcmin\nmagnification {:.4f}\ndeltay: {:.4f} arcsec".format(theta, gama, deltax*3600))


#### check with fixed rotation (from RC6 measurement), delta and gama
rot = 30.73/60/180*np.pi
mrot = np.array([[1, -rot], [rot, 1]])
slit_coord_rot = mrot.dot(slit_coord.T).T

#xn -x = x dgama  + dx 
mat = np.hstack((slit_coord_rot[:,[0]], np.ones((3,1))))
matinv =  np.linalg.pinv(mat)
sol = matinv.dot(slit_obs[:,0] - slit_coord_rot[:,0])
gama = 1 + sol[0]
deltax = sol[1]
print(gama, deltax*3600)
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[1, 3600]) #
#residual
xnew = gama*slit_coord_rot[:,0]  + deltax
print((slit_obs[:,0] - xnew)*3600)






# with using slit2 height scan 
# xn -x = x dgama - theta y + dx 
# yn -y = theta x  + y dgama + dy

row1   = np.hstack((slit_coord[0,::-1], np.array([0., 1.])))
row2_4 = np.hstack((slit_coord*[1.,-1], np.ones((3,1)), np.zeros((3,1))))
mat = np.vstack((row1, row2_4))
matinv =  np.linalg.inv(mat)
data = np.concatenate((slit_obs[[0],1] - slit_coord[[0],1], slit_obs[:,0] - slit_coord[:,0]))
sol = matinv.dot(data)
gama = 1 + sol[0]
theta_rad = sol[1]
deltax = sol[2]
deltay = sol[3]
theta = theta_rad*180/np.pi*60 #arcmin
print(gama, theta, deltax*3600, deltay*3600)
covar = matinv.dot(matinv.T)
# accuracy, assuming 1 arcsec measurement error
print(np.sqrt(np.diag(covar))/3600*[1, 180/np.pi*60, 3600]) #
#residual
xnew = gama*(np.cos(theta_rad)*slit_coord[:,0] -  np.sin(theta_rad)*slit_coord[:,1]) + deltax
print((slit_obs[:,0] - xnew)*3600)

xnew_rc6 = np.array([-0.14299218, -0.08832364,  0.06841098])
print((xnew - xnew_rc6)*3600)






# go to guiding stars:
######################
star_target_path = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/GuidingStars/'
#F4_stars = Table.read(star_target_path + "F4_guidingstars.txt", format='ascii')
F4_stars = Table.read(star_target_path + "F4_guidingstars.fits", format='fits')

mag = np.fmin(np.fmin(F4_stars['GAIA gband'].filled(99), F4_stars['SDSS gband'].filled(99)),
        F4_stars['SDSS rband'].filled(99))


# star #29
F4_stars[28]
guid_star_pos_ang = coordinates.SkyCoord(36.987137*u.deg, 0.402799*u.deg) #star 29

star_coord = G2UV.FieldP.world2local(guid_star_pos_ang)

delta_x_guider = G2UV.FOV_guider_coord.lon - gc_coord.lon + star_coord.lon
delta_y_guider = G2UV.FOV_guider_coord.lat - gc_coord.lat + star_coord.lat

theta_EL =  delta_x_guider/2./Elg
theta_CE = -delta_y_guider/2./CEg
print(theta_EL.to('arcsec'), theta_CE.to('arcsec'))

# star pixel seen on guider: x=376.5 y=281.5

# star #34
F4_stars[33]

guid_star_pos_ang = coordinates.SkyCoord(37.023638*u.deg, 0.390844*u.deg) #star 34

star_coord = G2UV.FieldP.world2local(guid_star_pos_ang)

delta_x_guider = G2UV.FOV_guider_coord.lon - gc_coord.lon + star_coord.lon
delta_y_guider = G2UV.FOV_guider_coord.lat - gc_coord.lat + star_coord.lat

theta_EL =  delta_x_guider/2./Elg
theta_CE = -delta_y_guider/2./CEg
print(theta_EL.to('arcsec'), theta_CE.to('arcsec'))

# angle cadre El 1149341 CE 1050325    at center

# star pixel seen on guider: x=326 y=131.5
#    x=327  y=132
    

###################3
# test with RC6
CEg = 1.0232 * 0.8898 # ? TBC
Elg = 1.00537  # ? TBC

F4_stars[28]
guid_star_pos_ang = coordinates.SkyCoord(36.987137*u.deg, 0.402799*u.deg) #star 29

star_coord = G2UV.FieldP.world2local(guid_star_pos_ang)

delta_x_guider = G2UV.FOV_guider_coord.lon - gc_coord.lon + star_coord.lon
delta_y_guider = G2UV.FOV_guider_coord.lat - gc_coord.lat + star_coord.lat

theta_EL =  delta_x_guider/2./Elg
theta_CE = -delta_y_guider/2./CEg
print(theta_EL.to('arcsec'), theta_CE.to('arcsec'))
#(<Angle [-116.41074448] arcsec>, <Angle [ 111.591738] arcsec>)

# star pixel seen on guider: x=375 y=312
# predicted
# 375.6064398
# 274.0377231 

