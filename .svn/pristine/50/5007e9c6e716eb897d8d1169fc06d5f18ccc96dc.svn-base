#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 14:03:25 2018

@author: dvibert
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from astropy.table import Table
from astropy.io import fits

path = '/data/ownCloud/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/XYCalib180904/'
#path = '/data/ownCloud/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/XYCalib180907/'
 
##############################
## Field 4
##############################

path_F4 = os.path.join(path, 'F4', 'F4Guider')
pa = 159.0

star_folders =  [ 'GC', 'GS1_s29', 'GS2_s14', 'GS3_s18', 'GS4_s34'] 
filename = 'MergedCatalog.csv'

#for sf in star_folders:
#    star_tab = Table.read(os.path.join(path_F4, sf, filename), format='csv')
#    scmask = (star_tab['ROTENC'] == pa)
#    plt.figure()
#    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
#    plt.plot(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'],'-*' )
    
   
stars = []    
for sf in star_folders:
    files = glob.glob(os.path.join(path_F4, sf, 'stack*.fits'))
    x = []
    y = []
    for f in files:
        with fits.open(f) as hdu:
            h = hdu[0].header
            if h['ROTENC'] == pa:
                cx0 = h['CX0']
                cy0 = h['CY0']
                x.append(cx0)
                y.append(cy0)
    x = np.array(x)
    y = np.array(y)
    xy = np.vstack([x,y])
    stars.append(xy)

plt.figure()
for i,s in enumerate(stars):
    plt.subplot(2,3,i+1)
    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
    plt.plot(s[0], s[1],'*' )
    m = s.mean(axis=1)
    plt.plot(m[0], m[1], 'o')

for s in stars:
    print("mean position: {}".format(np.array2string(s.mean(axis=1), formatter={'float_kind':lambda x: "%.1f" % x})))
    
##############################
## Field 3
##############################

pa = -121.0
path_F3 = os.path.join(path, 'F3', 'F3Guider')

star_folders =  [ 'GC', 'GS1', 'GS2', 'GS3', 'GS4'] 
filename = 'MergedCatalog.csv'

#for sf in star_folders:
#    star_tab = Table.read(os.path.join(path_F4, sf, filename), format='csv')
#    scmask = (star_tab['ROTENC'] == pa)
#    plt.figure()
#    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
#    plt.plot(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'],'-*' )
    
   
stars = []    
for sf in star_folders:
    files = glob.glob(os.path.join(path_F3, sf, 'stack*.fits'))
    x = []
    y = []
    for f in files:
        with fits.open(f) as hdu:
            h = hdu[0].header
            if h['ROTENC'] == pa:
                cx0 = h['CX0']
                cy0 = h['CY0']
                x.append(cx0)
                y.append(cy0)
    x = np.array(x)
    y = np.array(y)
    xy = np.vstack([x,y])
    stars.append(xy)

plt.figure()
for i,s in enumerate(stars):
    plt.subplot(2,3,i+1)
    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
    plt.plot(s[0], s[1],'*' )
    m = s.mean(axis=1)
    plt.plot(m[0], m[1], 'o')

for s in stars:
    print("mean position: {}".format(np.array2string(s.mean(axis=1), formatter={'float_kind':lambda x: "%.1f" % x})))
 
    
##############################
## Field 3 QSO5
##############################

pa = -121.0
path_QSO5 = os.path.join(path, 'QS0', 'QSOGuider')

star_folders =  [ 'GC', 'GS1', 'GS2'] 
filename = 'MergedCatalog.csv'

#for sf in star_folders:
#    star_tab = Table.read(os.path.join(path_F4, sf, filename), format='csv')
#    scmask = (star_tab['ROTENC'] == pa)
#    plt.figure()
#    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
#    plt.plot(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'],'-*' )
    
   
stars = []    
for sf in star_folders:
    files = glob.glob(os.path.join(path_QSO5, sf, 'stack*.fits'))
    x = []
    y = []
    for f in files:
        with fits.open(f) as hdu:
            h = hdu[0].header
            if h['ROTENC'] == pa:
                cx0 = h['CX0']
                cy0 = h['CY0']
                x.append(cx0)
                y.append(cy0)
    x = np.array(x)
    y = np.array(y)
    xy = np.vstack([x,y])
    stars.append(xy)

plt.figure()
for i,s in enumerate(stars):
    plt.subplot(2,3,i+1)
    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
    plt.plot(s[0], s[1],'*' )
    m = s.mean(axis=1)
    plt.plot(m[0], m[1], 'o')

for s in stars:
    print("mean position: {}".format(np.array2string(s.mean(axis=1), formatter={'float_kind':lambda x: "%.1f" % x})))
     
    
##############################
## Field 2
##############################

pa = -161.0
#path_F2 = os.path.join(path, 'F2', 'F2Guider')
path_F2 = os.path.join(path, 'F2_-161', 'F2Guider')

star_folders =  [ 'GC', 'GS1', 'GS2', 'GS3', 'GS4'] 
#filename = 'MergedCatalog.csv'

#for sf in star_folders:
#    star_tab = Table.read(os.path.join(path_F4, sf, filename), format='csv')
#    scmask = (star_tab['ROTENC'] == pa)
#    plt.figure()
#    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
#    plt.plot(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'],'-*' )
    
   
stars = []    
for sf in star_folders:
    files = glob.glob(os.path.join(path_F2, sf, 'stack*.fits'))
    x = []
    y = []
    for f in files:
        with fits.open(f) as hdu:
            h = hdu[0].header
            if h['ROTENC'] == pa:
                cx0 = h['CX0']
                cy0 = h['CY0']
                x.append(cx0)
                y.append(cy0)
    x = np.array(x)
    y = np.array(y)
    xy = np.vstack([x,y])
    stars.append(xy)

plt.figure()
for i,s in enumerate(stars):
    plt.subplot(2,3,i+1)
    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
    plt.plot(s[0], s[1],'*' )
    m = s.mean(axis=1)
    plt.plot(m[0], m[1], 'o')

for s in stars:
    print("mean position: {}".format(np.array2string(s.mean(axis=1), formatter={'float_kind':lambda x: "%.1f" % x})))
     
    
##############################
## Field 1
##############################

pa = 119.0
path_F1 = os.path.join(path, 'F1_119', 'F1Guider')

star_folders =  [ 'GC', 'GS1', 'GS2', 'GS3'] 
filename = 'MergedCatalog.csv'

#for sf in star_folders:
#    star_tab = Table.read(os.path.join(path_F4, sf, filename), format='csv')
#    scmask = (star_tab['ROTENC'] == pa)
#    plt.figure()
#    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
#    plt.plot(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'],'-*' )
    
   
stars = []    
for sf in star_folders:
    files = glob.glob(os.path.join(path_F1, sf, 'stack*.fits'))
    x = []
    y = []
    for f in files:
        with fits.open(f) as hdu:
            h = hdu[0].header
            if h['ROTENC'] == pa:
                cx0 = h['CX0']
                cy0 = h['CY0']
                x.append(cx0)
                y.append(cy0)
    x = np.array(x)
    y = np.array(y)
    xy = np.vstack([x,y])
    stars.append(xy)

plt.figure()
for i,s in enumerate(stars):
    plt.subplot(2,3,i+1)
    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
    plt.plot(s[0], s[1],'*' )
    m = s.mean(axis=1)
    plt.plot(m[0], m[1], 'o')

for s in stars:
    print("mean position: {}".format(np.array2string(s.mean(axis=1), formatter={'float_kind':lambda x: "%.1f" % x})))
         
    
##############################
## Field 1 QSO
##############################

pa = 119.0
path_QSO = os.path.join(path, 'QS0', 'QSOGuider')

star_folders =  [ 'GC', 'GS1', 'GS2', 'GS3'] 
filename = 'MergedCatalog.csv'

#for sf in star_folders:
#    star_tab = Table.read(os.path.join(path_F4, sf, filename), format='csv')
#    scmask = (star_tab['ROTENC'] == pa)
#    plt.figure()
#    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
#    plt.plot(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'],'-*' )
    
   
stars = []    
for sf in star_folders:
    files = glob.glob(os.path.join(path_QSO, sf, 'stack*.fits'))
    x = []
    y = []
    for f in files:
        with fits.open(f) as hdu:
            h = hdu[0].header
            if h['ROTENC'] == pa:
                cx0 = h['CX0']
                cy0 = h['CY0']
                x.append(cx0)
                y.append(cy0)
    x = np.array(x)
    y = np.array(y)
    xy = np.vstack([x,y])
    stars.append(xy)

plt.figure()
for i,s in enumerate(stars):
    plt.subplot(2,3,i+1)
    #plt.scatter(star_tab[scmask]['CX0'], star_tab[scmask]['CY0'])
    plt.plot(s[0], s[1],'*' )
    m = s.mean(axis=1)
    plt.plot(m[0], m[1], 'o')

for s in stars:
    print("mean position: {}".format(np.array2string(s.mean(axis=1), formatter={'float_kind':lambda x: "%.1f" % x})))
