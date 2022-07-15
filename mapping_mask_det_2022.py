#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 21:16:07 2018

@author: vpicouet
"""


from astropy.table import Table
from Calibration.mapping_mask_det import map_mask_detector, recompute_mask_pos, create_DS9regions
from Calibration.focustest import Focus
import numpy as np
import sys 




def plt_mapping(mask_table, mappings_linw,field):
    wl = 0.21382
    wl = 0.20619
    
    mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) & (mask_table["wavelength"]==wl)
    
    fig, (ax0,ax1) = plt.subplots(1,2)
    ax0.scatter(mask_table['ymask'][mask_mapping], mask_table['xmask'][mask_mapping],c=mask_table[mask_mapping]['wavelength'])
    ax0.set_xlabel("x_mm")
    ax0.set_ylabel("y_mm")
    ax1.set_xlabel("x PIX")
    ax1.set_ylabel("y PIX")
    ax1.scatter(mask_table['X_IMAGE'][mask_mapping],-mask_table['Y_IMAGE'][mask_mapping],c=mask_table[mask_mapping]['wavelength'],label="slit center")
    ax1.plot(mappings_linw.map(wl,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[0],-mappings_linw.map(0.20619,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[1],"+",c="orange",label="mask to det mapping")
    ax1.legend()
    plt.savefig("/Users/Vincent/Github/FireBallPipe/Calibration/Mappings/2022/%s.png"%(field))
    plt.show()
    return
def save_region_file(t,field,mappings_linw):
    wl = 0.21382
    # x, y  = mappings_linw.map(wl,x=t['x_mask_corrected'],y=t['y_mask_corrected'])
    x, y  = mappings_linw.map(wl,x=t['x_mm'],y=t['y_mm'])
    wl = 0.20619
    # x2, y2  = mappings_linw.map(wl,x=t['x_mask_corrected'],y=t['y_mask_corrected'])
    x2, y2  = mappings_linw.map(wl,x=t['x_mm'],y=t['y_mm'])
    wl = 0.20255
    # x3, y3  = mappings_linw.map(wl,x=t['x_mask_corrected'],y=t['y_mask_corrected'])
    x3, y3  = mappings_linw.map(wl,x=t['x_mm'],y=t['y_mm'])
    
    
    create_ds9_regions(
        [list(x),list(x2),list(x3)],
        [list(y),list(y2),list(y3)],
        radius=[15, 40] ,
        save=True,
        savename="/Users/Vincent/Github/FireBallPipe/Calibration/Mappings/2022/%s.reg"%(field),
        form=["box","box","box"],
        color=["red","yellow","blue"],
        ID=[list(t["Internal-count"]),list(t["Internal-count"]),list(t["Internal-count"])],
    )
    return


def fit_magnfication(mask_table, field):
    print(mask_table, field)
    xmm = mask_table["ymm"]
    ymm = -mask_table["xmm"]
    xpix = mask_table["X_IMAGE"]*0.013
    ypix = mask_table["Y_IMAGE"]*0.013
    def guider_to_detector(xy, x0, y0, theta, mx, my):
        x,y = xy
        # print(x)
        xrot = x * np.cos(theta*np.pi/180) + y * np.sin(theta*np.pi/180)
        yrot = y * np.cos(theta*np.pi/180) - x * np.sin(theta*np.pi/180) 
        xpix = (xrot+x0)*mx 
        ypix = (yrot+y0)*my  
        return (xpix,ypix)
    
    
    def get_mean_distance(x1, y1, x2, y2):
        return np.sqrt((x1 - x2)**2 + (y1 - y2)**2).mean()
    
    def err_func(params, x, y, xpix, ypix):
        x0, y0, theta, mx, my = params
        xpred,ypred = guider_to_detector((x,y), x0, y0, theta, mx, my)
        return get_mean_distance(xpix,ypix, xpred,ypred)
    
    from scipy.optimize import leastsq, minimize
    a = minimize(err_func, x0=(10,10,0,0.9,1), args=(xmm,ymm,xpix,ypix))
    print(a)
    
    ##%%
    fig, (ax0,ax1, ax2) = plt.subplots(1,3,figsize=(9,5))
    ax0.scatter(xmm,ymm,label="x0,x0=%0.1f, %0.1f \nθ=%0.2f \nmx,my=%0.2f,%0.2f"%(*a["x"],))
    for i, l in enumerate(mask_table):
        ax0.text(xmm[i],ymm[i],l["Internal-count"])
    ax0.set_xlabel("x_mm")#,c=mask_table['wavelength']
    ax0.set_ylabel("y_mm")
    ax1.set_xlabel("x PIX * 0.013")
    ax1.set_ylabel("y PIX * 0.013")
    ax2.set_xlabel("x PIX ")
    ax2.set_ylabel("y PIX ")
    ax1.scatter(xpix,ypix,label="slit center")
    x_pred,y_pred = guider_to_detector((xmm,ymm), *a["x"])
    ax1.plot(x_pred,y_pred,"+",c="orange",label="mask to det mapping")
    # ax1.plot(mappings_linw.map(wl,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[0],-mappings_linw.map(0.20619,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[1],"+",c="orange",label="mask to det mapping")
    ax1.legend(loc="lower right")
    ax0.legend(loc="lower right")
    print(np.max(abs(x_pred-xpix)))
    # scale=0.1
    q = ax2.quiver( mask_table["X_IMAGE"], mask_table["Y_IMAGE"],(x_pred-xpix)/0.013,(y_pred-ypix)/0.013,scale=0.02, scale_units='xy', angles='xy' )        
    ax2.quiverkey(q, .5, 0.1, 5, 'scale: %0.1f pix'%(5), coordinates='axes', color='r')                                                                  
    ax0.set_title(field)
    ax2.set_title("Residuals")
    fig.tight_layout()
    plt.show()
    return 

def create_mapping(field):
    lambda_lya = 0.121567
    # field="grid"
    # field="F4"
    ################ F1 +119
    ######################################################
    
    # path = '/data/FireBall/FTS-06-2018/180612/'
    #filename = 'image-000275-000284-Zinc-with_dark119-stack_table.csv'
    # imagename = 'image-000275-000284-Zinc-with_dark119-stack.fits'
    # t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/Detector_Mask_mappings/image-000275-000284-Zinc-with_dark119-stack_table.csv")
    # fn = "/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F1/F1_2022_5_-105.fits"
    try:
        t =Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/%s/%s_2022_6_-106.csv"%(field,field))
        t =Table.read("/Volumes/ExtremePro/LAM/FIREBALL/2022/DetectorData/220708/xycalib/diffuse/85_%s_-60.csv"%(field))
        t=Table.read("/Volumes/GoogleDrive-105248178238021002216/.shortcut-targets-by-id/1ZgB7kY-wf7meXrq8v-1vIzor75aRdLDn/FIREBall-2/FB2_2022/Detector_Data/220712/xycalib/stack_nanmean_1-5.csv")
        t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F3_60.csv")
        t = t[t["amp_x"]>5]
    except FileNotFoundError:
        t =Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/%s/%s_2022_6_-106_cat.fits"%(field,field))
        t= t[t["FLUX_MAX"]>200]
         
    # mask_table = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Slits/%s_new.csv"%(field))
    mask_table = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Targets/2022/targets_%s.txt"%(field),format="ascii")
    try:
        mask_table["xmask"],mask_table["ymask"] = mask_table["xmm"],mask_table["ymm"]
    except KeyError:
        mask_table["xmm"],mask_table["ymm"] = mask_table["x_mm"],mask_table["y_mm"]
        mask_table["xmask"],mask_table["ymask"] = mask_table["xmm"],mask_table["ymm"]
        
    mask_table1 = mask_table.copy()
    mask_table2 = mask_table.copy()
    mask_table["wavelength"] = 0.20255
    mask_table1["wavelength"] = 0.20619
    mask_table2["wavelength"] = 0.21382
    from astropy.table import vstack
    mask_table=vstack([mask_table,mask_table1,mask_table2])
    
    cols = [ 'line', 'x', 'y', 'amp_x', 'lx', 'x0', 'fwhm_x', 'off_x', 'amp_y', 'ly', 'y0', 'fwhm_y', 'off_y', 'smearing', 'fwhm_x_unsmear', 'lx_unsmear', 'x0_unsmear', 'amp_x_unsmear',"X_IMAGE","Y_IMAGE","X_IMAGE_unsmear"]#, 'l203', 'l214', 'l206']
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
        # print(n)
        if n>1:
            print(t[(t["name"]==line["Internal-count"])])
            # sys.exit()
        elif n==1:
            for col in cols:
                mask_table[col][i] = t[mask][0][col]
    # sys.exit()
    # map  subt['xmask'], subt['ymask'], subt['X_IMAGE'], subt['Y_IMAGE']
    # mask_table["X_IMAGE"] = mask_table["x"]      
    # mask_table["Y_IMAGE"] = mask_table["y"]      
    mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) #&(mask_table["Y_IMAGE"]>540) #& (mask_table["wavelength"]==0.20619)
    print(mask_table[mask_mapping])
    # plt.plot(mask_table[mask_mapping]["Y_IMAGE"],mask_table[mask_mapping]["Y_IMAGE"])
    mappings_linw, centers_linw =  map_mask_detector(mask_table[mask_mapping], bywave=False, deg=[1,3,3])
    mappings_linw.save('/Users/Vincent/Github/FireBallPipe/Calibration/Mappings/2022/mapping-mask-det-w-2022-5-%s.pkl'%(field))
    # plt_mapping(mask_table, mappings_linw,field)
    save_region_file(mask_table,field,mappings_linw)
    mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) & (mask_table["wavelength"]==0.21382)#0.20619)
    fit_magnfication(mask_table[mask_mapping],field)
    return

# curve_fit(guider_to_detector, np.array([mask_table[mask_mapping]["xmm"],mask_table[mask_mapping]["ymm"]]),np.array([mask_table[mask_mapping]["X_IMAGE"],mask_table[mask_mapping]["Y_IMAGE"]]),p0=None)

for field in ["F3"]:#,"F2","F3","F4","QSO"]:
    create_mapping(field)
    # sys.exit()
#%%

#%%





#%%

# sys.exit()



##%% PLOT DISTORTION MAP


# %%CREATE REGION FROM MAPPING

# t = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Targets/2022/targets_%s.txt"%(field),format="ascii")
# t = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Targets/2022/targets_F3.txt",format="csv")
# t = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Targets/2022/targets_QSO.txt",format="csv")

save_region_file(mask_table    )



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
