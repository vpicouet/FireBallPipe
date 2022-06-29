#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 21:16:07 2018

@author: dvibert
"""


F4 = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Slits/F4_2022.csv")

plt.plot(F4["xmm_L"],F4["ymm_L"],".")


for file in glob.glob("/Users/Vincent/Github/notebooks/clean_fields/*_220411.csv"):
    cat=Table.read(file)
    plt.figure(figsize(12,6))
    plt.plot(-cat["x_mm"],-cat["y_mm"],".")
    # plt.text(cat["x_mm"],cat["y_mm"],")
    plt.title(os.path.basename(file))
    plt.savefig(file.replace(".csv",".png"))
    plt.show()


#%%
from astropy.table import Table
from Calibration.mapping_mask_det import map_mask_detector, recompute_mask_pos, create_DS9regions
from Calibration.focustest import Focus
import numpy as np
import sys 
lambda_lya = 0.121567

################ F1 +119
######################################################

# path = '/data/FireBall/FTS-06-2018/180612/'
#filename = 'image-000275-000284-Zinc-with_dark119-stack_table.csv'
# imagename = 'image-000275-000284-Zinc-with_dark119-stack.fits'
t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/Detector_Mask_mappings/image-000275-000284-Zinc-with_dark119-stack_table.csv")
# fn = "/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F1/F1_2022_5_-105.fits"
t =Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F1/F1_2022_6_-106.csv")
t = t[t["amp_x"]>5]

mask_table = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Slits/F1_new.csv")
mask_table1 = mask_table.copy()
mask_table2 = mask_table.copy()
mask_table["wavelength"] = 0.20255
mask_table1["wavelength"] = 0.20619
mask_table2["wavelength"] = 0.21382
from astropy.table import vstack
mask_table=vstack([mask_table,mask_table1,mask_table2])

cols = [ 'line', 'x', 'y', 'amp_x', 'lx', 'x0', 'fwhm_x', 'off_x', 'amp_y', 'ly', 'y0', 'fwhm_y', 'off_y', 'smearing', 'fwhm_x_unsmear', 'lx_unsmear', 'x0_unsmear', 'amp_x_unsmear']#, 'l203', 'l214', 'l206']
for col in cols:
    mask_table[col] = np.nan

t["wavelength"]=np.nan
t["wavelength"][t["line"]==203] = 0.20255
t["wavelength"][t["line"]==206] = 0.20619
t["wavelength"][t["line"]==214] = 0.21382


# mask_table["X_IMAGE"] = np.nan
# mask_table["Y_IMAGE"] = np.nan
# mask = (mask_table==0.20255) | (mask_table==0.20619) | (mask_table==0.20619)
for i,line in enumerate(mask_table):
    mask = (t["name"]==line["Internal-count"]) & (t["wavelength"]==line["wavelength"])
    n=len(t[mask])
    print(n)
    if n>1:
        print(t[(t["name"]==line["Internal-count"])])
        # sys.exit()
    elif n==1:
        for col in cols:
            mask_table[col][i] = t[mask][0][col]
            
mask_table["X_IMAGE"] = mask_table["x"]      
mask_table["Y_IMAGE"] = mask_table["y"]      
mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) #& (mask_table["wavelength"]==0.20619)

mappings_linw, centers_linw =  map_mask_detector(mask_table[mask_mapping], bywave=False, deg=[1,3,3])
mappings_linw.save('/Users/Vincent/Github/FireBallPipe/Calibration/Mappings/mapping-mask-det-w-2022-5-F1.pkl')

sys.exit()

#%%
mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) & (mask_table["wavelength"]==0.20619)

fig, (ax0,ax1) = plt.subplots(1,2)
ax0.scatter(mask_table['ymask'][mask_mapping], mask_table['xmask'][mask_mapping],c=mask_table[mask_mapping]['wavelength'])
ax0.set_xlabel("x_mm")
ax0.set_ylabel("y_mm")
ax1.set_xlabel("x PIX")
ax1.set_ylabel("y PIX")

ax1.scatter(mask_table['X_IMAGE'][mask_mapping],-mask_table['Y_IMAGE'][mask_mapping],c=mask_table[mask_mapping]['wavelength'],label="slit center")
ax1.plot(mappings_linw.map(0.20619,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[0],-mappings_linw.map(0.20619,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[1],"+",c="orange",label="mask to det mapping")
ax1.legend()
plt.show()

#%%
d=DS9n()
regs = getregion(d, selected=True)
region = regs[0]
x, y = int(region.xc), int(region.yc)
w, h = int(region.w), int(region.h)
x_inf, x_sup, y_inf, y_sup = lims_from_region(region=region, coords=None)
image = d.get_pyfits()[0].data
n=10
subim3 = image[y_inf - n : y_sup + n, x_inf - n : x_sup + n]
ly,lx = subim3.shape
# xs = x - 1 + np.arange(-lx/ 2, lx / 2)
xs = x_inf -n + np.arange(lx+1)
ys = y_inf -n + np.arange(ly+1)
plt.imshow(subim3,extent=[xs.min(),xs.max(),ys.min(),ys.max()]);plt.colorbar()
plt.grid()
#%%
cat = Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F1/F1_2022_5_-105.csv")
create_ds9_regions(
    [cat["X_IMAGE_unsmear"]],
    # [cat["X_IMAGE"]],
    [cat["Y_IMAGE"]],
    # radius=[table_to_array(cat["h", "w"]).T],
    # radius=[np.array(cat["lx"]),np.array(cat["ly"])],
    radius=[np.array(cat["lx_unsmear"]),np.array(cat["ly"])],
    save=True,
    savename=tmp_region,
    form=["box"],
    color=cat["color"],
    ID=None,  # cat["name"],
)
d=DS9n()
d.set("regions %s" % (tmp_region))

#%%
t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/Detector_Mask_mappings/image-000275-000284-Zinc-with_dark119-stack_table.csv")


# F1 = Focus(filename = fn, 
#            threshold = [7], fwhm = [9,12.5], HumanSupervision=False,
#            quick=False, reversex=False, plot=False, source='Zn',
#            shape='slits',windowing=True, mask='F1', dist_min=30, date=180612, peak_threshold=50, min_fwhm=1.5)

# tabname = imagename.replace('.fits','_table.csv')
# #t=Table.read(path + tabname, format='csv')
# t["id_slit"] =t["name"] 
# t["id_slit"] = np.arange(len(t))#t["name"] 
# 
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




# #ds9reg_filename = tabname.replace('table.csv', 'mapped')
# ds9reg_filename = fn.replace('table.csv', 'mapped_w')

# create_DS9regions(xdetpix, ydetpix, form=['bpanda']*3, radius=10, save=True, 
#                       savename=path + ds9reg_filename, color = colors, ID=ID)

# # create Lya reg
# # filter redshifts
# w = lambda_lya * (1+z)
# xydetpix_lya = mappings_linw.map(w, xmask, ymask)
# ds9reg_filename = tabname.replace('table.csv', 'mapped_lya')
# create_DS9regions([xydetpix_lya[0]], [xydetpix_lya[1]], form=['bpanda'], radius=10, save=True, 
#                       savename=path + ds9reg_filename, color = ['white'], ID=[internalCount])




#%%

################ F2 -161
######################################################

#filename = 'image-000025-000034-Zinc-with_dark-161-stack_table.csv'
imagename = 'image-000025-000034-Zinc-with_dark-161-stack.fits'

F2 = Focus(filename = path + imagename, 
           threshold = [7], fwhm = [9,12.5], HumanSupervision=False,
           quick=False, reversex=False, plot=False, source='Zn',
           shape='slits',windowing=True, mask='F2', dist_min=30, date=180612, peak_threshold=50, min_fwhm=1.5)

tabname = fn.replace('.fits','_table.csv')
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
