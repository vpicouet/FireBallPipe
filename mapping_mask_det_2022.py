#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 21:16:07 2018

@author: vpicouet
"""

import matplotlib.pyplot as plt
from astropy.table import Table
from Calibration.mapping_mask_det import map_mask_detector, recompute_mask_pos#, create_DS9regions
from Calibration.focustest import Focus
import numpy as np
import sys 
from pyds9plugin.DS9Utils import create_ds9_regions
import glob, os
from datetime import datetime

def plt_mapping(mask_table, mappings_linw,field):
    wl = 0.21382
    wl = 0.20619
    
    mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) & (mask_table["wavelength"]==wl)
    # fig, (ax0,ax1) = plt.subplots(1,2)
    # ax0.scatter(mask_table['ymask'][mask_mapping], mask_table['xmask'][mask_mapping],c=mask_table[mask_mapping]['wavelength'])
    # ax0.set_xlabel("x_mm")
    # ax0.set_ylabel("y_mm")
    # ax1.set_xlabel("x PIX")
    # ax1.set_ylabel("y PIX")
    # ax1.scatter(mask_table['X_IMAGE'][mask_mapping],-mask_table['Y_IMAGE'][mask_mapping],c=mask_table[mask_mapping]['wavelength'],label="slit center")
    # ax1.plot(mappings_linw.map(wl,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[0],-mappings_linw.map(0.20619,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[1],"+",c="orange",label="mask to det mapping")
    # ax1.legend()
    # plt.savefig("/Users/Vincent/Github/FireBallPipe/Calibration/Mappings/2023/%s.png"%(field))
    # plt.show()
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
        savename="/Users/Vincent/Nextcloud/LAM/FIREBALL/FireBallPipe/Calibration/Mappings/2023/%s.reg"%(field),
        form=["box","box","box"],
        color=["red","yellow","blue"],
        ID=[list(t["Internal-count"]),list(t["Internal-count"]),list(t["Internal-count"])],
    )
    return




def guider_to_detector(xy, x0, y0, theta, mx, my):
    """
    x0,x0 is just center of rotation no 
    """

    x,y = xy
    # print(x)
    xrot = x * np.cos(theta*np.pi/180) + y * np.sin(theta*np.pi/180)
    yrot = y * np.cos(theta*np.pi/180) - x * np.sin(theta*np.pi/180) 
    xpix = (xrot+x0)*mx 
    ypix = (yrot+y0)*my  
    return (xpix,ypix)

def detector_to_guider(xy, x0, y0, theta, mx, my):
    """
    x0,x0 is just center of rotation no 
    """

    x,y = xy
    # print(x)
    xrot = (x/0.013/mx-x0) 
    yrot = (y/0.013/my-y0)  
    
    x_mask = x * np.cos(-theta*np.pi/180) + y * np.sin(-theta*np.pi/180)
    y_mask = y * np.cos(-theta*np.pi/180) - x * np.sin(-theta*np.pi/180) 

    return (x_mask,y_mask)

# def fit_magnfication(mask_table, field):
def fit_magnfication(xmm,ymm,xpix,ypix,id_, field="",x0=(30,15,0,0.89,1),file=None):
    # print(mask_table, field)
    # xmm = mask_table["ymm"]
    # ymm = -mask_table["xmm"]
    # xpix = mask_table["X_IMAGE"]*0.013
    # ypix = mask_table["Y_IMAGE"]*0.013

    
    def get_mean_distance(x1, y1, x2, y2):
        return np.sqrt((x1 - x2)**2 + (y1 - y2)**2).mean()
    
    def err_func(params, x, y, xpix, ypix):
        x0, y0, theta, mx, my = params
        xpred,ypred = guider_to_detector((x,y), x0, y0, theta, mx, my)
        return get_mean_distance(xpix,ypix, xpred,ypred)
    
    from scipy.optimize import leastsq, minimize
    a = minimize(err_func, x0=x0, args=(xmm,ymm,xpix,ypix))
    print(a)
    
    ##%%
    fig, (ax0,ax1, ax2) = plt.subplots(1,3,figsize=(12,5))
    ax0.scatter(xmm,ymm,label="x0,x0=%0.1f, %0.1f \nθ=%0.2f \nmx,my=%0.3f,%0.3f"%(*a["x"],))
    ax1.scatter(xpix,ypix,label="slit center")
    for i, l in enumerate(id_):
        ax0.text(xmm[i],ymm[i],l)#["Internal-count"])
    for i, l in enumerate(id_):
        ax1.text(xpix[i],ypix[i],l)#["Internal-count"])
    ax0.set_xlabel("y_mm")#,c=mask_table['wavelength']
    ax0.set_ylabel("-x_mm")
    ax1.set_xlabel("x PIX * 0.013")
    ax1.set_ylabel("y PIX * 0.013")
    ax2.set_xlabel("x PIX ")
    ax2.set_ylabel("y PIX ")
    x_pred,y_pred = guider_to_detector((xmm,ymm), *a["x"])
    ax1.plot(x_pred,y_pred,"+",c="orange",label="mask to det mapping")
    # ax1.plot(mappings_linw.map(wl,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[0],-mappings_linw.map(0.20619,x=mask_table['xmask'][mask_mapping],y=mask_table['ymask'][mask_mapping])[1],"+",c="orange",label="mask to det mapping")
    ax1.legend()#loc="lower right")
    ax0.legend()#loc="lower right")
    print(np.max(abs(x_pred-xpix)))
    # scale=0.1
    q = ax2.quiver( xpix/0.013, ypix/0.013,(x_pred-xpix)/0.013,(y_pred-ypix)/0.013,scale=0.02, scale_units='xy', angles='xy' )        
    ax2.quiverkey(q, .5, 0.1, 5, 'scale: %0.1f pix'%(5), coordinates='axes', color='r')                                                                  
    ax0.set_title(field)
    ax2.set_title("Residuals")
    fig.tight_layout()
    if file is not None:
        fig.savefig(file.replace(".csv",".png"))
    plt.show()
    return  a["x"]

def create_mapping(field, file=None,x1 = "ymm",y1 = "-xmm",x2="X_IMAGE*0.013",y2="Y_IMAGE*0.013", id_="Internal_count"):
    lambda_lya = 0.121567
    t=Table.read(file)
    # field="grid"
    # field="F4"
    ################ F1 +119
    ######################################################
    
    # path = '/data/FireBall/FTS-06-2018/180612/'
    #filename = 'image-000275-000284-Zinc-with_dark119-stack_table.csv'
    # imagename = 'image-000275-000284-Zinc-with_dark119-stack.fits'
    # t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/XYCalibration/Detector_Mask_mappings/image-000275-000284-Zinc-with_dark119-stack_table.csv")
    # fn = "/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F1/F1_2022_5_-105.fits"

    # if t is None:
    #     try:
    #         t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F3_60.csv")
    #         t = t[t["amp_x"]>5]
    #     except FileNotFoundError:
    #         t =Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/%s/%s_2022_6_-106_cat.fits"%(field,field))
    #         t= t[t["FLUX_MAX"]>200]
             
    # mask_table = Table.read("/Users/Vincent/Github/FireBallPipe/Calibration/Slits/%s_new.csv"%(field))
    
    mask_table = Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/FireBallPipe/Calibration/Targets/2022/targets_%s.csv"%(field),format="ascii")
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
    cols += [c for c in  t.colnames if "_IMAGE_" in c]
    # cols = t.colnames
    # print(np.ptp(mask_table["Y_IMAGE"][np.isfinite(mask_table["Y_IMAGE"])]))

    for col in cols:
        mask_table[col] = np.nan
    t["wavelength"]=np.nan
    t["wavelength"][t["line"]==203] = 0.20255
    t["wavelength"][t["line"]==206] = 0.20619
    t["wavelength"][t["line"]==214] = 0.21382
    
    
    # mask_table["X_IMAGE"] = np.nan
    # mask_table["Y_IMAGE"] = np.nan
    # mask = (mask_table==0.20255) | (mask_table==0.20619) | (mask_table==0.20619)
    # print(t["name"])
    t["name"] = t["name"].astype(str)
    for i,line in enumerate(mask_table):
        print(line["Internal-count"])
        mask = (t["name"]==str(line["Internal-count"])) & (t["wavelength"]==line["wavelength"])

        n=len(t[mask])
        # print(n, line["Internal-count"])
        if (n>1) :#| ((field=="F2")&(line["Internal-count"]=="18")):
            print(t[(str(t["name"])==str(line["Internal-count"]))])
            # sys.exit()
            mask_table[col][i] = np.mean(t[mask][col])
        elif n==1:
            for col in cols:
                mask_table[col][i] = t[mask][0][col]
    # sys.exit()
    # map  subt['xmask'], subt['ymask'], subt['X_IMAGE'], subt['Y_IMAGE']
    # mask_table["X_IMAGE"] = mask_table["x"]      
    # mask_table["Y_IMAGE"] = mask_table["y"]   
    # print(np.ptp(mask_table["Y_IMAGE"][np.isfinite(mask_table["Y_IMAGE"])]))
    mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) #&(mask_table["Y_IMAGE"]>540) #& (mask_table["wavelength"]==0.20619)
    # print(mask_table["Y_IMAGE"])
    # print(mask_table[mask_mapping])
    # plt.plot(mask_table[mask_mapping]["Y_IMAGE"],mask_table[mask_mapping]["Y_IMAGE"])
    print("Creating mapping: ", mask_table[mask_mapping])#, mask_table[mask_mapping])
    mappings_linw, centers_linw =  map_mask_detector(mask_table[mask_mapping], bywave=False, deg=[1,2,2])
    print("saving file: ", )
    mappings_linw.save('/Users/Vincent/Nextcloud/LAM/FIREBALL/FireBallPipe/Calibration/Mappings/2023/mask_to_det_mapping/mapping-mask-det-w-0-%s_%s.pkl'%(field,datetime.now().strftime("%y%m%d")))
    plt_mapping(mask_table, mappings_linw,field)
    save_region_file(mask_table,field,mappings_linw)
    mask_mapping =  np.isfinite(mask_table["Y_IMAGE"]) #& (mask_table["wavelength"]==0.21382)#0.20619)
    
    xmm = mask_table[mask_mapping]["ymm"]
    ymm = -mask_table[mask_mapping]["xmm"]
    xpix = mask_table[mask_mapping]["X_IMAGE"]*0.013
    ypix = mask_table[mask_mapping]["Y_IMAGE"]*0.013
    mask_table.rename_column("Internal-count","Internal_count")
    pd = mask_table[mask_mapping].to_pandas()
    
    # mag = fit_magnfication(xmm,ymm,xpix,ypix,mask_table[mask_mapping]["Internal-count"],field)
    mag = fit_magnfication(pd.eval(x1),pd.eval(y1),pd.eval(x2),pd.eval(y2),pd.eval(id_),field,file=file)
    # mag = fit_magnfication(mask_table[mask_mapping],field)
    # return mag#Table([os.path.basename(os.path.dirname(f)),os.path.basename(f)]+list(mag),names=["iteration","mask","x0","y0","theta","mx","my"])
    return [os.path.basename(os.path.dirname(file)),os.path.basename(file)]+list(mag)#   ,names=["iteration","mask","x0","y0","theta","mx","my"])


# mag= create_mapping("QSObright",file="/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/18_cold_-120_230819/QSObright_image000019.csv")#Table.read(f) )



#%%

if __name__ == "__main__":

    mag= create_mapping("QSO",file="/Users/Vincent/Library/CloudStorage/GoogleDrive-vp2376@columbia.edu/.shortcut-targets-by-id/1ZgB7kY-wf7meXrq8v-1vIzor75aRdLDn/FIREBall-2/FB2_2023/DOBC_data/230620/xy/image000006.csv")#Table.read(f) )
    mag= create_mapping("F3",file="/Users/Vincent/Library/CloudStorage/GoogleDrive-vp2376@columbia.edu/.shortcut-targets-by-id/1ZgB7kY-wf7meXrq8v-1vIzor75aRdLDn/FIREBall-2/FB2_2023/DOBC_data/230620/xy/image000008.csv")#Table.read(f) )
    
    
    sys.exit()
    mag= create_mapping(os.path.basename(f).split("_")[0],file=f)#Table.read(f) )
    
    
    sys.exit()
    
    mag= create_mapping(os.path.basename(f).split("_")[0],t=Table.read(f) )
    #%%
    # fit_magnfication(3600*EL*2,3600*CE*2,X,-Y,Y, field="",x0=(0,0,0,1,1))
    fit_magnfication(X,-Y,3600*EL*2,3600*CE*2,Y, field="",x0=(0,0,0,1,1))
    # fit_magnfication(EL,CE,X,-Y,Y, field="",x0=(0,0,0,1/3600,1/3600))
    
    #%%
    
    # curve_fit(guider_to_detector, np.array([mask_table[mask_mapping]["xmm"],mask_table[mask_mapping]["ymm"]]),np.array([mask_table[mask_mapping]["X_IMAGE"],mask_table[mask_mapping]["Y_IMAGE"]]),p0=None)
    
    
    # create_mapping("F2",t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F2_19.csv") )
    # create_mapping("F3",t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/F3_60.csv") )
    # create_mapping("QSO",t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/all_diffuse_illumination/FocusEvolution/QSO_-100.5_b.csv") )
    
    # create_mapping("F4",t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/10_cold_230320_0.55/F4_image000005.csv") )
    
    # create_mapping("F1",t=Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/11_cold_230407_0.55/F1_image000074.csv") )
    mags = Table(names=["iteration","mask","x0","y0","theta","mx","my"],dtype=["S20","S20",float,float,float,float,float])
    # for f in glob.glob("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/11_cold_230407_0.55/[FQ]*.csv"):
    # for f in glob.glob("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2022/**/[Q]*_6_*.csv"):
    # for f in glob.glob("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2022/*_cold_*/[FQ]*.csv"):
    for f in glob.glob("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/202?/16_cold_*/[Q]*.csv"):
        print(f)
        # for f in glob.glob("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/*_cold_*/F2*.csv"):
        mag= create_mapping(os.path.basename(f).split("_")[0].split(".")[0],file=f)#Table.read(f) )
        # try:
        #     mag= create_mapping(os.path.basename(f).split("_")[0],file=f)#Table.read(f) )
        #     mags.add_row(mag)
        # except ValueError as e:
        #     print("Do not work! ", e)
    mags["iter"] =  [ int(re.findall(r"\d+", m)[0]) for m in mags["iteration"]]
    mags["mask"] =  [ m.split("_")[0] for m in mags["mask"]]
    # mags.write("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/magnifications.csv",overwrite=True)
        # mags.append(create_mapping(os.path.basename(f).split("_")[0],t=Table.read(f) ))
        # sys.exit()    
    
    
    #%%
    for n in [5,10,13,37]:
        mags = Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/magnifications.csv")
        plt.figure(figsize=(6,10))
        for file in glob.glob("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/1?_cold_*/F1*.csv"):
            iteration, mask = os.path.basename(os.path.dirname(f)),os.path.basename(f)
            a= Table.read(file)
            iteration, mask = os.path.basename(os.path.dirname(file)),os.path.basename(file)
            a = a[a["X_IMAGE"]>1600]
            mag = mags[(mags["mask"]==mask)&(mags["iteration"]==iteration)]
            # xmask, ymask = detector_to_guider((a["X_IMAGE"],a["Y_IMAGE"]),mag["x0"],mag["y0"],mag["theta"],mag["mx"],mag["my"])
            # plt.plot(xmask, ymask,".")
            try:
                plt.plot(a["X_IMAGE"]-a["X_IMAGE"][a["name"]==str(n)],a["Y_IMAGE"]-a["Y_IMAGE"][a["name"]==str(n)],".",label=1)
            except ValueError:
                pass
            # plt.plot(a["X_IMAGE"]+mag["y0"],a["Y_IMAGE"]+mag["x0"],".")
        plt.show()
        
        
    #%%
    
    
    mags = Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/magnifications.csv")
    mags = mags[mags["iter"]>10]
    plt.figure(figsize=(10,5))
    for i, (mask,c) in enumerate(zip(["F1","F2","F3","F4","QSO"],["r","g","orange","b","k"])):  
        for it, ms in zip(np.unique(mags["iter"]),["D","X","P","^","s","$\odot$"]*5):  
        # for it, ms in zip(np.arange(0,17),["o","X","P","^","s","v"]*5):  
            submag = mags[(mags["iter"]==it)&(mags["mask"]==mask)]
            plt.scatter(submag["mx"],submag["my"],marker=ms, c=c)
            # if it==16:
            #     plt.text(submag["mx"],submag["my"],mask + ":\n%0.3f\n%0.3f"%(submag["mx"],submag["my"]))
            if i==0:
                # plt.plot(0,0,ms,c="k",label=it)
                submag = mags[(mags["iter"]==it)]
                plt.scatter(None,None,marker=ms,c="k",label="Iter %i: mx~%0.3f±%0.1f%%, my~%0.3f±%0.1f%%"%(it, np.mean(submag["mx"]),100*(submag["mx"].ptp())/submag["mx"].mean()/2,np.mean(submag["my"]),100*(submag["mx"].ptp())/submag["mx"].mean()/2))
    
        submag = mags[(mags["mask"]==mask)]
        # plt.plot(0,0,"o",c=c,label=mask + ": mx~%0.3f±%0.1f%%, my~%0.3f±%0.1f%%"%(np.mean(submag["mx"]),100*(submag["mx"].ptp())/submag["mx"].mean()/2,np.mean(submag["my"]),100*(submag["mx"].ptp())/submag["mx"].mean()/2))
        plt.scatter(None,None,marker="o",c=c,label=mask + ": mx~%0.3f±%0.1f%%, my~%0.3f±%0.1f%%"%( mags[(mags["iter"]==it)&(mags["mask"]==mask)]["mx"],100*(submag["mx"].ptp())/submag["mx"].mean()/2, mags[(mags["iter"]==it)&(mags["mask"]==mask)]["my"],100*(submag["mx"].ptp())/submag["mx"].mean()/2))
    plt.plot(mags[(mags["iter"]==16)]["mx"],mags[(mags["iter"]==16)]["my"],"k:")
         
    plt.legend(title="Maks/iterations",fontsize=8)
    # plt.xlim((0.86,0.9))
    # plt.ylim((1,1.025))
    plt.grid()
    plt.xlabel("Xmm det/Xmm mask magnification")
    plt.ylabel("Ymm det/Ymm mask magnification")
    plt.title("Det/mask magnification evolution with Mask/iteration")
    
    
    
    #%%
    
    
    
    mags = Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/magnifications.csv")
    mags = mags[mags["iter"]>10]
    mags["mx"] *= 1.26 # arcsec/pix
    mags["my"] *= 1.08 # arcsec/pix
    plt.figure(figsize=(10,5))
    for i, (mask,c) in enumerate(zip(["F1","F2","F3","F4","QSO"],["r","g","orange","b","k"])):  
        # for it, ms in zip(np.unique(mags["iter"]),["o","x","+",".","s","^"]*5):  
        for it, ms in zip(np.unique(mags["iter"]),["D","X","P","^","s","$\odot$"]*5):  
            if i==0:
                submag = mags[(mags["iter"]==it)]
                plt.scatter(None,None,marker=ms,c="k",label="Iter %i: mx~%0.3f±%0.1f%%, my~%0.3f±%0.1f%%"%(it, np.mean(submag["mx"]),100*(submag["mx"].ptp())/submag["mx"].mean()/2,np.mean(submag["my"]),100*(submag["mx"].ptp())/submag["mx"].mean()/2))
    
            submag = mags[(mags["iter"]==it)&(mags["mask"]==mask)]
            plt.scatter(submag["mx"],submag["my"],marker=ms, c=c)
        submag = mags[(mags["mask"]==mask)]
        plt.scatter(None,None,marker="o",c=c,label=mask + ": mx~%0.3f±%0.1f%%, my~%0.3f±%0.1f%%"%( mags[(mags["iter"]==it)&(mags["mask"]==mask)]["mx"],100*(submag["mx"].ptp())/submag["mx"].mean()/2, mags[(mags["iter"]==it)&(mags["mask"]==mask)]["my"],100*(submag["mx"].ptp())/submag["mx"].mean()/2))
    plt.plot(mags[(mags["iter"]==16)]["mx"],mags[(mags["iter"]==16)]["my"],"k:")
    plt.legend(title="Maks/iterations",fontsize=8)
    # plt.xlim((1.09,1.125))
    # plt.ylim((1.08,1.105))
    plt.grid()
    plt.xlabel("X sky/mask magnification [asec/pix]")
    plt.ylabel("Y sky/mask  magnification [asec/pix]")
    plt.title("Sky to Mask magnification evolution with Mask/iteration")
    
    
    
    #%%
    
    
    EL = np.array([0.32922	,0.32922	,0.32922	,0.44034	,0.44034])*2*0.9969
    CE = np.array([-0.25557,-0.13057	,0.00832	,0.00832	,-0.25557])*2*1.0180
    X =  np.array([1304.6	,1307.4	    ,1328.8	    ,1965.7	    ,1940.7])
    Y =  np.array([1868.6	,1034.5	    ,107.7	    ,117.5	    ,1882.7])
    
    
    # from scipy.optimize import leastsq, minimize
    # a = minimize(err_func, x0=(30,15,0,0.89,1), args=(EL,CE,X,Y))
    
    # fit_magnfication(3600*EL,3600*CE,X,-Y,Y, field="",x0=(0,0,0,1,1))
    fit_magnfication(X,-Y,3600*EL,3600*CE,Y, field="",x0=(0,0,0,1,1))
        
    #%%
    # Now  I need to derotate the mask to det.
    # So get the of the slits in the det
    # I derotate by the theta I found (and ferify then that theta is 0)
    # Then I compute the difference between 2 far away slits
    # Multiply this pixel value to get degrees on sky.
    # On compare to the difference 
    # or just plot all the slits.
    
    
    # mx,my = 1.105,1.1091
    mx,my = 1.26/1.105, 1.08/1.091
    t=np.deg2rad(-0.13)#-0.13
    rac,decc = 36.9049,	0.65245
    # mx,my = 1.26*0.877, 1.08*1.01
    # mx,my = 1/0.877, 1/1.01
    mx, my = 1.26, 1.08# pixel at detector become sky to detector
    
    Elg = 0.9969 # 2023 #1.0090 2022 # Elg = 1.00379 # 2018
    CEg = 1.018  # 2023 #1.0187 2022 # CEg = 1.02928 # 2018
    
    # mx,my =  1.105/0.877,1.1091/1.01
    # mx,my =  1.105/1.01,1.1091/0.877
    # we want sky to mask so we want to multiply by a mask to det
    f = Table.read("/Users/Vincent/Library/CloudStorage/GoogleDrive-vp2376@columbia.edu/.shortcut-targets-by-id/1ZgB7kY-wf7meXrq8v-1vIzor75aRdLDn/FIREBall-2/FB2_2023/instrument_alignment_focusing/XY_calibration/FireBallPipe/Calibration/Targets/2022/targets_F4.csv")
    p = Table.read("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/16_cold_230612/F4_image000008.csv")
    rotmat = np.array([[np.cos(t), -np.sin(t)], 
                       [np.sin(t),  np.cos(t)]])
    p["X_IMAGE_rot"],p["Y_IMAGE_rot"] = (rotmat @ np.lib.recfunctions.structured_to_unstructured(p["X_IMAGE","Y_IMAGE"].as_array()).T)#.T
    p.write("/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/16_cold_230612/F4_image000008.csv",overwrite=True)
    fig, (ax0, ax1) = plt.subplots(1,2,sharex=True,sharey=True)
    # ax0.plot( (f['RA'] - np.mean(f['RA'])) * np.cos((f["DEC"]-np.mean(f["DEC"]))*np.pi/180),f["DEC"]-np.mean(f["DEC"]),".")
    # ax1.plot(-p["Y_IMAGE"],p["X_IMAGE"],".")
    
    ax0.plot( f["DEC"]-decc,(f['RA'] - rac) * np.cos((f["DEC"]-decc)*np.pi/180),".")
    ax1.plot((p["X_IMAGE_rot"]-p["X_IMAGE_rot"].mean())*mx/3600,-(p["Y_IMAGE_rot"]-p["Y_IMAGE_rot"].mean())*mx/3600,"k.")
    # ax0.plot((p["X_IMAGE_rot"]-p["X_IMAGE_rot"].mean())*mx/3600,-(p["Y_IMAGE_rot"]-p["Y_IMAGE_rot"].mean())*mx/3600,"k.")
    ax0.grid()
    ax1.grid()
    ax0.set_ylim((-0.3,0.2))
        
    
    # f="/Users/Vincent/Nextcloud/LAM/FIREBALL/all_diffuse_illumination/2023/16_cold_230612/F4_image000008.csv"
    # mag= create_mapping(os.path.basename(f).split("_")[0],file=f,x1 = "RA",y1 = "DEC",x2="X_IMAGE_rot*%s/3600"%(mx),y2="Y_IMAGE_rot*%s/3600"%(my), id_="Internal_count" )
    
    
    # fit_magnfication( f["DEC"]-decc,(f['RA'] - rac) * np.cos((f["DEC"]-decc)*np.pi/180),3600*EL*2,3600*CE*2,Y, field="",x0=(0,0,0,1,1))
    
    
        