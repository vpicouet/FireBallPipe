#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 21:28:20 2017

@author: dvibert
"""


from __future__ import division, print_function

import numpy as np
from astropy import coordinates
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.table import Table
import matplotlib.patches as patches

from guider2UV import Guider2UV

cloudpath = '/home/dvibert/ownCloud/FIREBALL/'
cloudpath = '/Users/Vincent/Nextcloud/FIREBALL/'
#path_SCGUI01 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170908_SC_GUI01/"
#G2UV= Guider2UV(filename=path_SCGUI01  + 'Guider2UV_F1.new.pkl')

#path_SC_GUI02 = cloudpath + 'Tests-at-FortSumner/170909_SC_GUI02/'
path = cloudpath + 'TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/'
#detpix = 12e-3 #mm

gc = np.array([640, 540]) # guider center

#path_SCGUI03 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/170923_SC_GUI03/"

#######################################################
# F1
##########################################################
#G2UV = Guider2UV(filename=path_SC_GUI02 + 'Guider2UV_F1_nogamma.new.pkl')

G2UV = Guider2UV(filename=path + 'F1_180907.pkl')

target_filename = cloudpath  + 'Target_selection/targets_F1.txt'
targets = Table.read(target_filename, format='ascii')

tx = np.array(targets['xmm']) + np.array([[-.12],[ +.12]])
ty = np.array(targets['ymm'])[np.newaxis,:]
tx, ty = np.broadcast_arrays(tx, ty)

#get stars
star_target_path = cloudpath + 'Target_selection/GuidingStars/'
F1_stars = Table.read(star_target_path + "F1_guidingstars.fits", format='fits')
mag = np.fmin(np.fmin(F1_stars['GAIA gband'].filled(99), F1_stars['SDSS gband'].filled(99)),
        F1_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(F1_stars['RA']*u.deg, F1_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)


# plot in mask coord  
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(ty, tx,'-b',linewidth=3.0)
plt.xlim(plt.xlim()[::-1])
plt.ylim(plt.ylim()[::-1])
plt.xlabel('y mm (Dec axis)')
plt.ylabel('x mm (Ra axis)')
plt.title('Field mask 1')        
for s in range(len(targets)):
    plt.text(targets['ymm'][s], targets['xmm'][s], targets['Internal-count'][s], color='k')
plt.show()

# plot in guider pix

gx = np.zeros(tx.shape)
gy = np.zeros(ty.shape)
for i,(xi,yi) in  enumerate(np.nditer([tx, ty])):    
    g = G2UV.SienceMask2guider(np.array([xi, yi]), angle=False)
    gx.ravel()[i] = g[0][0]
    gy.ravel()[i] = g[1][0]

 
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(gx, gy,'-b',linewidth=3.0)
#plt.xlim(plt.xlim()[::-1])
#plt.ylim(plt.ylim()[::-1])
plt.xlabel('x guider pix (Dec axis)')
plt.ylabel('y guider pix (-Ra axis)')
plt.title('Field mask 1')
for s in range(len(targets)):
    plt.text(gx[1,s], gy[1,s], targets['Internal-count'][s], color='k')
ax.add_patch(patches.Rectangle((0, 0),1280,1080,
        fill=False, color='r'))
plt.xlim([0,plt.xlim()[1]])
plt.ylim(plt.ylim()[::-1])
# plot guider center
plt.plot(gc[0], gc[1], 'xg')
plt.text(gc[0], gc[1], 'Gcenter')
# plot Field center
plt.plot(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'og')
plt.text(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'Fcenter')
    
# plot stars    
select = mag <=12    
bright = mag <=12
select = np.in1d(F1_stars['Internal count'], [8,25,31])
plt.scatter(guider_star_pos[0][bright], guider_star_pos[1][bright], 100, 
            marker="*", color='k', label='predicted mag<=12')
plt.scatter(guider_star_pos[0][select], guider_star_pos[1][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(bright)[0]:
    plt.text(guider_star_pos[0][s], guider_star_pos[1][s], str(F1_stars['Internal count'][s]))
    plt.text(guider_star_pos[0][s], guider_star_pos[1][s]-50, str(mag[s])+'m')
  
plt.show()


#######################################################
# F1 QSOMgII
##########################################################
#G2UV = Guider2UV(filename=path_SC_GUI02 + 'Guider2UV_F1_nogamma.new.pkl')

G2UV = Guider2UV(filename=path + 'F1_QSOMgII_180907.pkl')

target_filename = cloudpath  + 'Target_selection/targets_F1.txt'
targets = Table.read(target_filename, format='ascii')

tx = np.array(targets['xmm']) + np.array([[-.12],[ +.12]])
ty = np.array(targets['ymm'])[np.newaxis,:]
tx, ty = np.broadcast_arrays(tx, ty)

#get stars

QSO_stars_radec = np.array([[19.014794, 32.334118],
                        [19.13107, 32.31136],
                        [19.134911, 32.369127]])

QSO_stars = Table(QSO_stars_radec, names=[ 'RA', 'DEC'])
QSO_stars['Internal count'] = [1,2,3]

stars = [1,2,3]

#star_target_path = cloudpath + 'Target_selection/GuidingStars/'
#F1_stars = Table.read(star_target_path + "F1_guidingstars.fits", format='fits')
#mag = np.fmin(np.fmin(F1_stars['GAIA gband'].filled(99), F1_stars['SDSS gband'].filled(99)),
#        F1_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(QSO_stars['RA']*u.deg, QSO_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)


# plot in mask coord  
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(ty, tx,'-b',linewidth=3.0)
plt.xlim(plt.xlim()[::-1])
plt.ylim(plt.ylim()[::-1])
plt.xlabel('y mm (Dec axis)')
plt.ylabel('x mm (Ra axis)')
plt.title('Field mask 1')        
for s in range(len(targets)):
    plt.text(targets['ymm'][s], targets['xmm'][s], targets['Internal-count'][s], color='k')
plt.show()

# plot in guider pix

gx = np.zeros(tx.shape)
gy = np.zeros(ty.shape)
for i,(xi,yi) in  enumerate(np.nditer([tx, ty])):    
    g = G2UV.SienceMask2guider(np.array([xi, yi]), angle=False)
    gx.ravel()[i] = g[0][0]
    gy.ravel()[i] = g[1][0]

 
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(gx, gy,'-b',linewidth=3.0)
#plt.xlim(plt.xlim()[::-1])
#plt.ylim(plt.ylim()[::-1])
plt.xlabel('x guider pix (Dec axis)')
plt.ylabel('y guider pix (-Ra axis)')
plt.title('Field mask 1 - QSO MgII stars')
for s in range(len(targets)):
    plt.text(gx[1,s], gy[1,s], targets['Internal-count'][s], color='k')
ax.add_patch(patches.Rectangle((0, 0),1280,1080,
        fill=False, color='r'))
plt.xlim([0,plt.xlim()[1]])
plt.ylim(plt.ylim()[::-1])
# plot guider center
plt.plot(gc[0], gc[1], 'xg')
plt.text(gc[0], gc[1], 'Gcenter')
# plot Field center
plt.plot(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'og')
plt.text(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'Fcenter')
    
# plot stars    
plt.scatter(guider_star_pos[0], guider_star_pos[1], 100, 
            marker="*", color='r', label='')
for s in  range(3):
    plt.text(guider_star_pos[0][s], guider_star_pos[1][s], str(QSO_stars['Internal count'][s]))
    #plt.text(guider_star_pos[0][s], guider_star_pos[1][s]-50, str(mag[s])+'m')
  
plt.show()


# slit scann;
# gc - 46 - 51 -FHC- 20 - 10 - FHC - s8 - s31 - FHC - s25
    
#######################################################
# F2
##########################################################

#G2UV = Guider2UV(filename=path_SC_GUI02 + 'Guider2UV_F2_nogamma.new.pkl')

G2UV = Guider2UV(filename=path + 'F2_180907.pkl')

target_filename = cloudpath  + 'Target_selection/targets_F2.txt'
targets = Table.read(target_filename, format='ascii')

tx = np.array(targets['xmm']) + np.array([[-.12],[ +.12]])
ty = np.array(targets['ymm'])[np.newaxis,:]
tx, ty = np.broadcast_arrays(tx, ty)

#get stars
star_target_path = cloudpath + 'Target_selection/GuidingStars/'
F2_stars = Table.read(star_target_path + "F2_guidingstars.fits", format='fits')
mag = np.fmin(np.fmin(F2_stars['GAIA gband'].filled(99), F2_stars['SDSS gband'].filled(99)),
        F2_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(F2_stars['RA']*u.deg, F2_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)


# plot in mask coord  
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(ty, tx,'-b',linewidth=3.0)
plt.xlim(plt.xlim()[::-1])
plt.ylim(plt.ylim()[::-1])
plt.xlabel('y mm (Dec axis)')
plt.ylabel('x mm (Ra axis)')
plt.title('Field mask 2')        
for s in range(len(targets)):
    plt.text(targets['ymm'][s], targets['xmm'][s], targets['Internal-count'][s], color='k')
plt.show()


# plot in guider pix

gx = np.zeros(tx.shape)
gy = np.zeros(ty.shape)
for i,(xi,yi) in  enumerate(np.nditer([tx, ty])):    
    g = G2UV.SienceMask2guider(np.array([xi, yi]), angle=False)
    gx.ravel()[i] = g[0][0]
    gy.ravel()[i] = g[1][0]

 
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(gx, gy,'-b',linewidth=3.0)
#plt.xlim(plt.xlim()[::-1])
#plt.ylim(plt.ylim()[::-1])
plt.xlabel('x guider pix (Dec axis)')
plt.ylabel('y guider pix (-Ra axis)')
plt.title('Field mask 2')
for s in range(len(targets)):
    plt.text(gx[1,s], gy[1,s], targets['Internal-count'][s], color='k')
ax.add_patch(patches.Rectangle((0, 0),1280,1080,
        fill=False, color='r'))
plt.xlim([0,plt.xlim()[1]])
plt.ylim(plt.ylim()[::-1])
# plot guider center
plt.plot(gc[0], gc[1], 'xg')
plt.text(gc[0], gc[1], 'Gcenter')
# plot Field center
plt.plot(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'og')
plt.text(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'Fcenter')
    
# plot stars    
select = mag <=12    
bright = mag <=12
select = np.in1d(F2_stars['Internal count'], [21,36,40,46])
plt.scatter(guider_star_pos[0][bright], guider_star_pos[1][bright], 100, 
            marker="*", color='k', label='predicted mag<=12')
plt.scatter(guider_star_pos[0][select], guider_star_pos[1][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(bright)[0]:
    plt.text(guider_star_pos[0][s], guider_star_pos[1][s], str(F2_stars['Internal count'][s]))
   # plt.text(guider_star_pos[0][s], guider_star_pos[1][s]-50, str(mag[s])+'m')
plt.show()
    
# slit scann;
# gc ???- s36 - s40 - s46 - FHC - s21
    


#######################################################
# F3
##########################################################

#G2UV = Guider2UV(filename=path_SC_GUI02 + 'Guider2UV_F3_nogamma.new.pkl')
G2UV = Guider2UV(filename=path + 'F3_180904.pkl')

target_filename = cloudpath  + 'Target_selection/targets_F3.txt'
targets = Table.read(target_filename, format='ascii')

tx = np.array(targets['xmm']) + np.array([[-.12],[ +.12]])
ty = np.array(targets['ymm'])[np.newaxis,:]
tx, ty = np.broadcast_arrays(tx, ty)

#get stars
star_target_path = cloudpath + 'Target_selection/GuidingStars/'
F3_stars = Table.read(star_target_path + "F3_guidingstars.fits", format='fits')
mag = np.fmin(np.fmin(F3_stars['GAIA gband'].filled(99), F3_stars['SDSS gband'].filled(99)),
        F3_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(F3_stars['RA']*u.deg, F3_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)


# plot in mask coord  
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(ty, tx,'-b',linewidth=3.0)
plt.xlim(plt.xlim()[::-1])
plt.ylim(plt.ylim()[::-1])
plt.xlabel('y mm (Dec axis)')
plt.ylabel('x mm (Ra axis)')
plt.title('Field mask 3')        
for s in range(len(targets)):
    plt.text(targets['ymm'][s], targets['xmm'][s], targets['Internal-count'][s], color='k')
plt.show()


# plot in guider pix

gx = np.zeros(tx.shape)
gy = np.zeros(ty.shape)
for i,(xi,yi) in  enumerate(np.nditer([tx, ty])):    
    g = G2UV.SienceMask2guider(np.array([xi, yi]), angle=False)
    gx.ravel()[i] = g[0][0]
    gy.ravel()[i] = g[1][0]

 
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(gx, gy,'-b',linewidth=3.0)
#plt.xlim(plt.xlim()[::-1])
#plt.ylim(plt.ylim()[::-1])
plt.xlabel('x guider pix (Dec axis)')
plt.ylabel('y guider pix (-Ra axis)')
plt.title('Field mask 3')
for s in range(len(targets)):
    plt.text(gx[1,s], gy[1,s], targets['Internal-count'][s], color='k')
ax.add_patch(patches.Rectangle((0, 0),1280,1080,
        fill=False, color='r'))
plt.xlim([0,plt.xlim()[1]])
plt.ylim(plt.ylim()[::-1])
# plot guider center
plt.plot(gc[0], gc[1], 'xg')
plt.text(gc[0], gc[1], 'Gcenter')
# plot Field center
plt.plot(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'og')
plt.text(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'Fcenter')
    
# plot stars    
select = mag <=12    
bright = mag <=12
select = np.in1d(F3_stars['Internal count'], [5,9,20])
plt.scatter(guider_star_pos[0][bright], guider_star_pos[1][bright], 100, 
            marker="*", color='k', label='predicted mag<=12')
plt.scatter(guider_star_pos[0][select], guider_star_pos[1][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(bright)[0]:
    plt.text(guider_star_pos[0][s], guider_star_pos[1][s], str(F3_stars['Internal count'][s]))
#    plt.text(guider_star_pos[0][s], guider_star_pos[1][s]-50, str(mag[s])+'m')
plt.show()
   
# slit scann;
# gc - ????????? - s9 - s20 - FHC - s5 
    
#######################################################
# F4
##########################################################

#G2UV = Guider2UV(filename=path_SC_GUI02 + 'Guider2UV_F4_nogamma.new.pkl')

#G2UV = Guider2UV(filename=path_SC_GUI02 + 'Guider2UV_F4_nogamma.new.pkl')

G2UV = Guider2UV(filename=path + 'F4_180904.pkl')

target_filename = cloudpath  + 'Target_selection/targets_F4.txt'
targets = Table.read(target_filename, format='ascii')

tx = np.array(targets['xmm']) + np.array([[-.12],[ +.12]])
ty = np.array(targets['ymm'])[np.newaxis,:]
tx, ty = np.broadcast_arrays(tx, ty)

#get stars
star_target_path = cloudpath + 'Target_selection/GuidingStars/'
F4_stars = Table.read(star_target_path + "F4_guidingstars.fits", format='fits')
mag = np.fmin(np.fmin(F4_stars['GAIA gband'].filled(99), F4_stars['SDSS gband'].filled(99)),
        F4_stars['SDSS rband'].filled(99))

coords = coordinates.SkyCoord(F4_stars['RA']*u.deg, F4_stars['DEC']*u.deg)
guider_star_pos = G2UV.SienceMask2guider(coords, world=True, angle=False)


# plot in mask coord  
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(ty, tx,'-b',linewidth=3.0)
plt.xlim(plt.xlim()[::-1])
plt.ylim(plt.ylim()[::-1])
plt.xlabel('y mm (Dec axis)')
plt.ylabel('x mm (Ra axis)')
plt.title('Field mask 4')        
for s in range(len(targets)):
    plt.text(targets['ymm'][s], targets['xmm'][s], targets['Internal-count'][s], color='k')

plt.show()


# plot in guider pix

gx = np.zeros(tx.shape)
gy = np.zeros(ty.shape)
for i,(xi,yi) in  enumerate(np.nditer([tx, ty])):    
    g = G2UV.SienceMask2guider(np.array([xi, yi]), angle=False)
    gx.ravel()[i] = g[0][0]
    gy.ravel()[i] = g[1][0]

 
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(gx, gy,'-b',linewidth=3.0)
#plt.xlim(plt.xlim()[::-1])
#plt.ylim(plt.ylim()[::-1])
plt.xlabel('x guider pix (Dec axis)')
plt.ylabel('y guider pix (-Ra axis)')
plt.title('Field mask 4')
for s in range(len(targets)):
    plt.text(gx[1,s], gy[1,s], targets['Internal-count'][s], color='k')
ax.add_patch(patches.Rectangle((0, 0),1280,1080,
        fill=False, color='r'))
plt.xlim([0,plt.xlim()[1]])
plt.ylim(plt.ylim()[::-1])
# plot guider center
plt.plot(gc[0], gc[1], 'xg')
plt.text(gc[0], gc[1], 'Gcenter')
# plot Field center
plt.plot(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'og')
plt.text(G2UV.FOV_center_guider_pix[0], G2UV.FOV_center_guider_pix[1], 'Fcenter')
    
# plot stars    
bright = mag <=12
select = np.in1d(F4_stars['Internal count'], [29,14,18,34])
plt.scatter(guider_star_pos[0][bright], guider_star_pos[1][bright], 100, 
            marker="*", color='k', label='predicted mag<=12')
plt.scatter(guider_star_pos[0][select], guider_star_pos[1][select], 100, 
            marker="*", color='r', label='predicted mag<=12')
for s in  np.nonzero(bright)[0]:
    plt.text(guider_star_pos[0][s], guider_star_pos[1][s], str(F4_stars['Internal count'][s]))
#    plt.text(guider_star_pos[0][s], guider_star_pos[1][s]-50, str(mag[s])+'m')
plt.show()
    
# slit scann;
# gc -?? - s18 - s29 - s34 - FHC - s14    

