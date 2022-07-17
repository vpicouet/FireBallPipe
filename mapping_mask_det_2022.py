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
    a = minimize(err_func, x0=(30,15,0,0.89,1), args=(xmm,ymm,xpix,ypix))
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

def create_mapping(field, t=None):
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
    if t is None:
        try:
            # t =Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/%s/%s_2022_6_-106.csv"%(field,field))
            # t =Table.read("/Volumes/ExtremePro/LAM/FIREBALL/2022/DetectorData/220708/xycalib/diffuse/85_%s_-60.csv"%(field))
            # t=Table.read("/Volumes/GoogleDrive-105248178238021002216/.shortcut-targets-by-id/1ZgB7kY-wf7meXrq8v-1vIzor75aRdLDn/FIREBall-2/FB2_2022/Detector_Data/220712/xycalib/stack_nanmean_1-5.csv")
            t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F3_60.csv")
            # t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F4_-20.csv")
            t = t[t["amp_x"]>5]
        except FileNotFoundError:
            t =Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/%s/%s_2022_6_-106_cat.fits"%(field,field))
            t= t[t["FLUX_MAX"]>200]
             
    # mask_table = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Slits/%s_new.csv"%(field))
    
    mask_table = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Targets/2022/targets_%s.csv"%(field),format="ascii")
    try:
        mask_table["xmask"],mask_table["ymask"] = mask_table["x_mask_corrected"],mask_table["_mask_corrected"]
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
    print(mask_table["Y_IMAGE"])
    print(mask_table[mask_mapping])
    # plt.plot(mask_table[mask_mapping]["Y_IMAGE"],mask_table[mask_mapping]["Y_IMAGE"])
    mappings_linw, centers_linw =  map_mask_detector(mask_table[mask_mapping], bywave=False, deg=[1,3,3])
    mappings_linw.save('/Users/Vincent/Github/FireBallPipe/Calibration/Mappings/2022/mapping-mask-det-w-2022-5-%s.pkl'%(field))
    plt_mapping(mask_table, mappings_linw,field)
    save_region_file(mask_table,field,mappings_linw)
    mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) #& (mask_table["wavelength"]==0.21382)#0.20619)
    fit_magnfication(mask_table[mask_mapping],field)
    return

# curve_fit(guider_to_detector, np.array([mask_table[mask_mapping]["xmm"],mask_table[mask_mapping]["ymm"]]),np.array([mask_table[mask_mapping]["X_IMAGE"],mask_table[mask_mapping]["Y_IMAGE"]]),p0=None)


create_mapping("F2",t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F2_19.csv") )
create_mapping("F3",t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F3_60.csv") )
create_mapping("QSO",t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/QSO_-100.5.csv") )
