#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 21:16:07 2018

@author: dvibert
"""

from astropy.table import Table
from Calibration.mapping_mask_det import map_mask_detector, recompute_mask_pos, create_DS9regions
from Calibration.focustest import Focus

lambda_lya = 0.121567

################ F1 +119
######################################################

path = '/data/FireBall/FTS-06-2018/180612/'
#filename = 'image-000275-000284-Zinc-with_dark119-stack_table.csv'
imagename = 'image-000275-000284-Zinc-with_dark119-stack.fits'

F1 = Focus(filename = path + imagename, 
           threshold = [7], fwhm = [9,12.5], HumanSupervision=False,
           quick=False, reversex=False, plot=False, source='Zn',
           shape='slits',windowing=True, mask='F1', dist_min=30, date=180612, peak_threshold=50, min_fwhm=1.5)

tabname = imagename.replace('.fits','_table.csv')
#t=Table.read(path + tabname, format='csv')
t = F1.table

xmask, ymask, internalCount, z = recompute_mask_pos(t, 'F1')

mappings, centers = map_mask_detector(t)
print(centers)


# remove 213nm  do linear interp in lambda
only_202_206 = (t['wavelength'] == 0.20255) | (t['wavelength'] == 0.20619)
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,5,2])
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,5,5])
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,2,2])
mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,3,3])

#mappings.plot()


mappings_linw.save('mapping-mask-det-w-1806012-F1.pkl')

# create mapped region file
wcolors = {0.20255:'blue', 0.20619:'green', 0.21382:'red'}
colors = [wcolors[i] for i in mappings.w]
#colors = ['cyan']*3
ID =[internalCount]*3
xdetpix = []
ydetpix = []
for w in mappings.w:
    #xydetpix = mappings.map(w, xmask, ymask)
    xydetpix = mappings_linw.map(w, xmask, ymask)
    xdetpix.append(xydetpix[0])
    ydetpix.append(xydetpix[1])

#ds9reg_filename = tabname.replace('table.csv', 'mapped')
ds9reg_filename = tabname.replace('table.csv', 'mapped_w')

create_DS9regions(xdetpix, ydetpix, form=['bpanda']*3, radius=10, save=True, 
                      savename=path + ds9reg_filename, color = colors, ID=ID)

# create Lya reg
# filter redshifts
w = lambda_lya * (1+z)
xydetpix_lya = mappings_linw.map(w, xmask, ymask)
ds9reg_filename = tabname.replace('table.csv', 'mapped_lya')
create_DS9regions([xydetpix_lya[0]], [xydetpix_lya[1]], form=['bpanda'], radius=10, save=True, 
                      savename=path + ds9reg_filename, color = ['white'], ID=[internalCount])



################ F2 -161
######################################################

#filename = 'image-000025-000034-Zinc-with_dark-161-stack_table.csv'
imagename = 'image-000025-000034-Zinc-with_dark-161-stack.fits'

F2 = Focus(filename = path + imagename, 
           threshold = [7], fwhm = [9,12.5], HumanSupervision=False,
           quick=False, reversex=False, plot=False, source='Zn',
           shape='slits',windowing=True, mask='F2', dist_min=30, date=180612, peak_threshold=50, min_fwhm=1.5)

tabname = imagename.replace('.fits','_table.csv')
#t=Table.read(path + tabname, format='csv')
t = F2.table

xmask, ymask, internalCount, z = recompute_mask_pos(t, 'F2')

mappings, centers = map_mask_detector(t)
print(centers)


# remove 213nm  do linear interp in lambda
only_202_206 = (t['wavelength'] == 0.20255) | (t['wavelength'] == 0.20619)
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,5,2])
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,5,5])
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,2,2])
mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,3,3])

#mappings.plot()

#mappings.save('mapping-mask-det-1806012-F2.pkl')
mappings_linw.save('mapping-mask-det-w-1806012-F2.pkl')

# create mapped region file
wcolors = {0.20255:'blue', 0.20619:'green', 0.21382:'red'}
colors = [wcolors[i] for i in mappings.w]
#colors = ['cyan']*3
ID =[internalCount]*3
xdetpix = []
ydetpix = []
for w in mappings.w:
    #xydetpix = mappings.map(w, xmask, ymask)
    xydetpix = mappings_linw.map(w, xmask, ymask)
    xdetpix.append(xydetpix[0])
    ydetpix.append(xydetpix[1])

#ds9reg_filename = tabname.replace('table.csv', 'mapped')
ds9reg_filename = tabname.replace('table.csv', 'mapped_w')
create_DS9regions(xdetpix, ydetpix, form=['bpanda']*3, radius=10, save=True, 
                      savename=path + ds9reg_filename, color = colors, ID=ID)

# create Lya reg
# filter redshifts
w = lambda_lya * (1+z)
xydetpix_lya = mappings_linw.map(w, xmask, ymask)
ds9reg_filename = tabname.replace('table.csv', 'mapped_lya')
create_DS9regions([xydetpix_lya[0]], [xydetpix_lya[1]], form=['bpanda'], radius=10, save=True, 
                      savename=path + ds9reg_filename, color = ['white'], ID=[internalCount])

################ F3 -121
######################################################

#filename = 'image-000075-000084-Zinc-with_dark-121-stack_table.csv'
imagename = 'image-000075-000084-Zinc-with_dark-121-stack.fits'

F3 = Focus(filename = path + imagename, 
           threshold = [7], fwhm = [9,12.5], HumanSupervision=False,
           quick=False, reversex=False, plot=False, source='Zn',
           shape='slits',windowing=True, mask='F3', dist_min=30, date=180612, peak_threshold=50, min_fwhm=1.5)

tabname = imagename.replace('.fits','_table.csv')
#t=Table.read(path + tabname, format='csv')
t = F3.table

xmask, ymask, internalCount, z = recompute_mask_pos(t, 'F3')

mappings, centers = map_mask_detector(t)
print(centers)


# remove 213nm  do linear interp in lambda
only_202_206 = (t['wavelength'] == 0.20255) | (t['wavelength'] == 0.20619)
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,5,2])
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,5,5])
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,2,2])
mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,3,3])

#mappings.plot()


mappings_linw.save('mapping-mask-det-w-1806012-F3.pkl')

# create mapped region file
wcolors = {0.20255:'blue', 0.20619:'green', 0.21382:'red'}
colors = [wcolors[i] for i in mappings.w]
#colors = ['cyan']*3
ID =[internalCount]*3
xdetpix = []
ydetpix = []
for w in mappings.w:
    #xydetpix = mappings.map(w, xmask, ymask)
    xydetpix = mappings_linw.map(w, xmask, ymask)
    xdetpix.append(xydetpix[0])
    ydetpix.append(xydetpix[1])

#ds9reg_filename = tabname.replace('table.csv', 'mapped')
ds9reg_filename = tabname.replace('table.csv', 'mapped_w')
create_DS9regions(xdetpix, ydetpix, form=['bpanda']*3, radius=10, save=True, 
                      savename=path + ds9reg_filename, color = colors, ID=ID)

# create Lya reg
# filter redshifts
w = lambda_lya * (1+z)
xydetpix_lya = mappings_linw.map(w, xmask, ymask)
ds9reg_filename = tabname.replace('table.csv', 'mapped_lya')
create_DS9regions([xydetpix_lya[0]], [xydetpix_lya[1]], form=['bpanda'], radius=10, save=True, 
                      savename=path + ds9reg_filename, color = ['white'], ID=[internalCount])

################ F4 159
######################################################

#filename = 'image-000325-000334-Zinc-with_dark159-stack_table.csv'
imagename = 'image-000325-000334-Zinc-with_dark159-stack.fits'

F4 = Focus(filename = path + imagename, 
           threshold = [7], fwhm = [9,12.5], HumanSupervision=False,
           quick=False, reversex=False, plot=False, source='Zn',
           shape='slits',windowing=True, mask='F4', dist_min=30, date=180612, peak_threshold=50, min_fwhm=1.5)

tabname = imagename.replace('.fits','_table.csv')
#t=Table.read(path + tabname, format='csv')
t = F4.table

xmask, ymask, internalCount, z = recompute_mask_pos(t, 'F4')

mappings, centers = map_mask_detector(t)
print(centers)


# remove 213nm  do linear interp in lambda
only_202_206 = (t['wavelength'] == 0.20255) | (t['wavelength'] == 0.20619)
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,5,2])
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,5,5])
#mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,2,2])
mappings_linw, centers_linw =  map_mask_detector(t[only_202_206], bywave=False, deg=[1,3,3])

#mappings.plot()


mappings_linw.save('mapping-mask-det-w-1806012-F4.pkl')

# create mapped region file
wcolors = {0.20255:'blue', 0.20619:'green', 0.21382:'red'}
colors = [wcolors[i] for i in mappings.w]
#colors = ['cyan']*3
ID =[internalCount]*3
xdetpix = []
ydetpix = []
for w in mappings.w:
    #xydetpix = mappings.map(w, xmask, ymask)
    xydetpix = mappings_linw.map(w, xmask, ymask)
    xdetpix.append(xydetpix[0])
    ydetpix.append(xydetpix[1])

#ds9reg_filename = tabname.replace('table.csv', 'mapped')
ds9reg_filename = tabname.replace('table.csv', 'mapped_w')
create_DS9regions(xdetpix, ydetpix, form=['bpanda']*3, radius=10, save=True, 
                      savename=path + ds9reg_filename, color = colors, ID=ID)


w = lambda_lya * (1+z)
xydetpix_lya = mappings_linw.map(w, xmask, ymask)
ds9reg_filename = tabname.replace('table.csv', 'mapped_lya')
create_DS9regions([xydetpix_lya[0]], [xydetpix_lya[1]], form=['bpanda'], radius=10, save=True, 
                      savename=path + ds9reg_filename, color = ['white'], ID=[internalCount])
