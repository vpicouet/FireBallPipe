#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 18:59:57 2017

@author: vincent
"""
from __future__ import division, print_function


import os, sys
import numpy as np
import glob

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.patches import Ellipse

from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.integrate import dblquad, quad#, fixed_quad, quadrature
from scipy import ndimage, linalg, special#, signal, misc

from astropy import convolution
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.visualization import AsinhStretch, AsymmetricPercentileInterval, ImageNormalize     # MinMaxInterval, SqrtStretch,
from numpy.polynomial.polynomial import polyval2d#, polyval3d
from mpl_toolkits.mplot3d import axes3d
from pkg_resources  import resource_filename
from pyds9plugin.DS9Utils import create_ds9_regions
try:
    from astroML.crossmatch import crossmatch
except:
    pass

# from .mapping import Mapping    


#from mapping import Mapping    
#
#x = np.array([1683,1854,1866,2036,1849,2021])
#y = np.array([1260,1260,1040,1040,614,614])
#z = np.array([13.143,13.081,13.118,13.068,13.131,13.118])
##plt.figure()
#X,Y,Z, ax, Cdet = fit_quadratic_curve(x,y,z,n=1,order=1)
#plt.show()
#

#loadDS9Image(filename = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018/' + 'stack8103716_pa+078_2018-06-11T06-21-24.fits',Internet =False)
#loadDS9Image(filename = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018/stack8104066_pa+000_2018-06-11T06-23-09.fits',Internet =False)
#loadDS9Image(filename='/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/image-000275-000284-Zinc-with_dark119-stack.fits')
#loadDS9Image(filename='/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018/stack7802709_pa+119_2018-06-10T22-03-58.fits')
#loadDS9Image(filename = '/Users/Vincent/Nextcloud/FIREBALL/Tests-at-FortSumner/170923_SC_GUI03/Detector/170924/image000001.fits')
#


#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018/'
#def PlotFocus2DGuider(path=path, pa=-161,n=22, motors=np.linspace(11.95,14.45,11), dist_max=30, Plot=False, sigma=True,order=1,starsx=None,starsy=None):
#    detfocus = 6
#    if pa==119:
#        focus = 8.5
#    if pa==-121:
#        focus = 6
#    if pa==-161:
#        focus = 9
#    if pa==159:
#        focus = 7.5
#    npa = 'pa%+03d'%(pa)
#    files = glob.glob(path + '/*' + npa + '*.csv')[-11:]
#    stars = []
#    fwhms = []
#    fwhmsvar = []
#    size = []      
#    for i,file in enumerate(files):
#        table = Table.read(file)
#        stars.append(np.array([table['X_IMAGE'],table['Y_IMAGE']])) 
#        fwhms.append(np.sqrt(np.square(table['FWHM_x']) + np.square(table['FWHM_y']))) 
#        fwhmsvar.append(np.sqrt(np.square(table['FWHM_var_x']) + np.square(table['FWHM_var_y']))) 
#        size.append(len(stars[-1][0]))
#    if Plot:
#        plt.figure(figsize=(8,8))  
#        for i in range(len(stars)):
#            print (i,len(stars[i][0]))
#            plt.plot(stars[i][0],stars[i][1],'X',label = '{}'.format(i))
#            plt.legend()
#        plt.show()
#    imageWithMaxStars = np.argmax(size)
#    xc,yc = stars[imageWithMaxStars][0],stars[imageWithMaxStars][1]
#    fwhmx = np.full((len(xc),11),np.nan)
#    fwhmxvar = np.full((len(xc),11),np.nan)
#    fwhmx[:,imageWithMaxStars] = fwhms[imageWithMaxStars][1]
#    fwhmxvar[:,imageWithMaxStars] = fwhmsvar[imageWithMaxStars][1]
#    distance=[]
#    index=[]
#    for i in range(len(stars)):
#        dist,idx = crossmatch(np.array([stars[imageWithMaxStars][0],stars[imageWithMaxStars][1]]).T,np.array([stars[i][0],stars[i][1]]).T,max_distance=dist_max)
#        distance.append(dist)
#        index.append(idx)
#        
#    
#    for i in range(len(stars)):
#        mask = index[i] <  stars[i][0].size
#        fwhmx[mask,i] = fwhms[i][index[i][mask]]
#        fwhmxvar[mask,i] = fwhmsvar[i][index[i][mask]]
#    if Plot:
#        plt.figure()
#        plt.hist(fwhmx)#,bins = np.linspace(-1,15,40))
#        plt.figtext(0.6,0.6,'mean = %0.3f \nsigma = %0.3f'%(np.nanmean(fwhmx),np.sqrt(np.nanvar(fwhmx))))
#        plt.show()
#    n=3
#    cache2sigma = (fwhmx > np.nanmean(fwhmx) + n * np.nanstd(fwhmx)) | (fwhmx < np.nanmean(fwhmx) - n * np.nanstd(fwhmx)) | (fwhmx == 0.0)
#    fwhmx[cache2sigma] = np.nan
#    if Plot:
#        plt.figure()
#        plt.plot(fwhmx.T)
#        plt.show()
#    
#    f = lambda x,a,b,c: a * (x-b)**2 + c#a * np.square(x) + b * x + c
#    center, popt,pcov, center_sig = [], [], [], [] 
#    for star in range (len(fwhmx)):#
#        idx =  np.isfinite(fwhmx[star])
#        if np.nansum(idx)>6:
#           opt, cov = curve_fit(f, np.arange(11)[idx], fwhmx[star][idx], (1., 5,1),  bounds=([0,-10,0],[5,20,20]),sigma = np.sqrt(fwhmxvar[star][idx]))
#                
#        else:
#            opt = [np.nan,np.nan,np.nan]
#            cov = np.full((3,3),np.nan)
#        popt.append(opt)
#        pcov.append(cov)
#        center.append(opt[1])
#        center_sig.append(np.sqrt(cov[1,1]))
#    
#    
#    center = np.array(center)
#    center_sig = np.array(center_sig)
#    center_mask = (center < 0) |  (center > 10) |  (np.isnan(center)) 
#    if sigma:
#        x,y,z, ax, C = fit_quadratic_curve(stars[imageWithMaxStars][0][~center_mask],stars[imageWithMaxStars][1][~center_mask],center[~center_mask],sigma_z=None,n=10,order=order)#center_sig[~center_mask]
#    else:
#        x,y,z, ax, C = fit_quadratic_curve(stars[imageWithMaxStars][0][~center_mask],stars[imageWithMaxStars][1][~center_mask],center[~center_mask],sigma_z=center_sig[~center_mask]+1e-1,n=10,order=order)#center_sig[~center_mask]
#
#    Zfocus = np.ones((z.shape))
#    ax.plot_surface(x, y, Zfocus * focus, rstride=1, cstride=1, alpha=0.2,color='orange',label='Guider focus autocol')
#    ax.plot_surface(x, y, Zfocus * detfocus, rstride=1, cstride=1, alpha=0.2,color='b',label='Detector focus autocol')
#    for xstar, ystar in zip(starsx,starsy):
#        if order==1:
#            defocus = C[0]*xstar + C[1]*ystar + C[2] 
#        if order==2:
#            defocus = np.dot(np.c_[1, xstar, ystar, xstar*ystar, xstar**2, ystar**2], C)
#        ax.plot(np.ones(2)*xstar, np.ones(2)*ystar, np.linspace(focus,defocus,2),linewidth=6,c='black')
#    if order==1:
#        ax.plot(np.ones(2)*903, np.ones(2)*685, np.linspace(focus,C[0]*903 + C[1]*685 + C[2],2),linewidth=6,c='red')
#    if order==2:
#        defocus = np.dot(np.c_[1, 903, 685, 903*685, 903**2, 685**2], C)
#        ax.plot(np.ones(2)*903, np.ones(2)*685, np.linspace(focus,defocus,2),linewidth=6,c='red')
#    plt.title(npa)
#    plt.show()
#    return x,y,z,C
#starsxF2 = [624.1,731.8,496.2,298.1]	
#starsyF2 = [655.7,310.9,226.4,177.7]
#x,y,z,C = PlotFocus2DGuider(path=path, pa=-161,sigma=True,order=2,starsx=starsxF2,starsy=starsyF2)
#files = glob.glob('/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/Autocoll/Detector/1/*.fits')[3:7]    
#PlotFocus2DGuider()

#def giveOffsetPlane(C,value,center):
    
#b = [1150,647]
#a = [1279,222]  
def findpinhole(a,b):
    a, b = np.array(a), np.array(b) 
    if a[0]<b[0]:
        holex = 0.55 * a[0]+ 0.45 *b[0]
        holey = 0.55 * a[1]+ 0.45 *b[1]
    elif a[0]>b[0]:
        holex = 0.55 * b[0]+ 0.45 *a[0]
        holey = 0.55 * b[1]+ 0.45 *a[1] 
    return np.array([holex,holey])
#findpinhole(a,b)    


def PlotFocus2DGuider(path=1, pa=-161,n=22,PlanDetector=None, motors=np.linspace(11.95,14.45,11), dist_max=30, Plot=False, sigma=True,order=1,starsx=None,starsy=None):
    npa = 'pa%+03d'%(pa)
    files = glob.glob(path + '/*' + npa + '*table.fits')[-11:]
    print (files)
    stars = []
    fwhms = []
    fwhmsvar = []
    size = []    
    actuator = []
    for i,file in enumerate(files):
        header = fits.open(file)[0].header
        actuator.append(header['LINAENC'])
        table = Table.read(file)
        stars.append(np.array([table['X_IMAGE'],table['Y_IMAGE']])) 
        fwhms.append(np.sqrt(np.square(table['FWHM_x']) + np.square(table['FWHM_y']))) 
        fwhmsvar.append(np.sqrt(np.square(table['FWHM_var_x']) + np.square(table['FWHM_var_y']))) 
        size.append(len(stars[-1][0]))
        #print(stars)
    print(size)
    if Plot:
        plt.figure(figsize=(8,8))  
        for i in range(len(stars)):
            print (i,len(stars[i][0]))
            plt.plot(stars[i][0],stars[i][1],'X',label = '{}'.format(i))
            plt.legend()
        plt.show()
    imageWithMaxStars = np.argmax(size)
    xc,yc = stars[imageWithMaxStars][0],stars[imageWithMaxStars][1]
    fwhmx = np.full((len(xc),11),np.nan)
    fwhmxvar = np.full((len(xc),11),np.nan)
    fwhmx[:,imageWithMaxStars] = fwhms[imageWithMaxStars][1]
    fwhmxvar[:,imageWithMaxStars] = fwhmsvar[imageWithMaxStars][1]
    distance=[]
    index=[]
    for i in range(len(stars)):
        dist,idx = crossmatch(np.array([stars[imageWithMaxStars][0],stars[imageWithMaxStars][1]]).T,np.array([stars[i][0],stars[i][1]]).T,max_distance=dist_max)
        distance.append(dist)
        index.append(idx)
        
    
    for i in range(len(stars)):
        mask = index[i] <  stars[i][0].size
        fwhmx[mask,i] = fwhms[i][index[i][mask]]
        fwhmxvar[mask,i] = fwhmsvar[i][index[i][mask]]
    if Plot:
        plt.figure()
        plt.hist(fwhmx)#,bins = np.linspace(-1,15,40))
        plt.figtext(0.6,0.6,'mean = %0.3f \nsigma = %0.3f'%(np.nanmean(fwhmx),np.sqrt(np.nanvar(fwhmx))))
        plt.show()
    n=3
    cache2sigma = (fwhmx > np.nanmean(fwhmx) + n * np.nanstd(fwhmx)) | (fwhmx < np.nanmean(fwhmx) - n * np.nanstd(fwhmx)) | (fwhmx == 0.0)
    fwhmx[cache2sigma] = np.nan
    f = lambda x,a,b,c: a * (x-b)**2 + c#a * np.square(x) + b * x + c
    if Plot:
        plt.figure()
        plt.plot(actuator,fwhmx.T)
        plt.plot(actuator, f(np.array(actuator),5,13,7),'-o')
        plt.show()
    
    center, popt,pcov, center_sig = [], [], [], [] 
    for star in range (len(fwhmx)):#
        idx =  np.isfinite(fwhmx[star])
        if np.nansum(idx)>6:
           opt, cov = curve_fit(f, np.array(actuator)[idx], fwhmx[star][idx], (5, 13,7),  bounds=([0,5,0],[10,20,20]),sigma = np.sqrt(fwhmxvar[star][idx]))
#           opt, cov = curve_fit(f, np.arange(11)[idx], fwhmx[star][idx], (1., 5,1),  bounds=([0,-10,0],[5,20,20]),sigma = np.sqrt(fwhmxvar[star][idx]))
                
        else:
            opt = [np.nan,np.nan,np.nan]
            cov = np.full((3,3),np.nan)
        popt.append(opt)
        pcov.append(cov)
        center.append(opt[1])
        center_sig.append(np.sqrt(cov[1,1]))
    
    
    center = np.array(center)
    center_sig = np.array(center_sig)
#    center_mask = (center < 0) |  (center > 10) |  (np.isnan(center)) 
    center_mask = (center < 10) |  (center > 15) |  (np.isnan(center)) 
    if sigma:
        x,y,z, ax, C = fit_quadratic_curve(stars[imageWithMaxStars][0][~center_mask],stars[imageWithMaxStars][1][~center_mask],center[~center_mask],sigma_z=None,n=10,order=order)#center_sig[~center_mask]
    else:
        x,y,z, ax, C = fit_quadratic_curve(stars[imageWithMaxStars][0][~center_mask],stars[imageWithMaxStars][1][~center_mask],center[~center_mask],sigma_z=center_sig[~center_mask]+1e-1,n=10,order=order)#center_sig[~center_mask]

    if pa==119:
        foc = 18.25#8.5
    if pa==-121:
        foc = 18.99#6
    if pa==-161:
        foc = 18.25#9
    if pa==159:
        foc = 18.50#7.5

    if  PlanDetector is not None:  
        offsetDet = (18.996 - PlanDetector[0]*903 + PlanDetector[1]*685)
        offsetGuider = (foc - PlanDetector[0]*903 + PlanDetector[1]*685)#18.25#8.5
        detfocus = PlanDetector[0]*x + PlanDetector[1]*y + offsetDet
        focus = PlanDetector[0]*x + PlanDetector[1]*y + offsetGuider



    Zfocus = np.ones((z.shape))
    if order==1:
        autocf = C[0]*903 + C[1]*685 + C[2]
#        ax.plot(np.ones(2)*903, np.ones(2)*685, np.linspace(focus,focus,2),linewidth=6,c='red')
        ax.scatter(np.ones(2)*903, np.ones(2)*685, np.linspace(autocf,autocf,2), c='b', s=100)
    if order==2:
        autocf = np.dot(np.c_[1, 903, 685, 903*685, 903**2, 685**2], C)
#        ax.plot(np.ones(2)*903, np.ones(2)*685, np.linspace(focus,focus,2),linewidth=6,c='red')
        ax.scatter(np.ones(2)*903, np.ones(2)*685, np.linspace(autocf,autocf,2), c='b', s=100)
    diffInfAutoc = (foc-autocf)
    for xstar, ystar in zip(starsx,starsy):
        if order==1:
            defocus = C[0]*xstar + C[1]*ystar + C[2] 
        if order==2:
            defocus = np.dot(np.c_[1, xstar, ystar, xstar*ystar, xstar**2, ystar**2], C)
#        ax.plot(np.ones(2)*xstar, np.ones(2)*ystar, np.linspace(focus,defocus,2),linewidth=6,c='black')
        focstar = PlanDetector[0]*xstar + PlanDetector[1]*ystar + offsetDet
        ax.plot(np.ones(2)*xstar, np.ones(2)*ystar, np.linspace(focstar-diffInfAutoc,defocus,2),linewidth=6,c='black')
    #if PlanDetector is None:
    ax.plot_surface(x, y, Zfocus * focus - diffInfAutoc, rstride=1, cstride=1, alpha=0.2,color='orange',label='Guider focus autocol')
    ax.plot_surface(x, y, Zfocus * detfocus - diffInfAutoc, rstride=1, cstride=1, alpha=0.2,color='b',label='Detector focus autocol')
#    else:
#        Z = PlanDetector[0]*x + PlanDetector[1]*y + PlanDetector[2]
#        ax.plot_surface(x, y, Z - diffInfAutoc, rstride=1, cstride=1, alpha=0.2,color='orange',label='Guider focus autocol')
#        ax.plot_surface(x, y, Z - diffInfAutoc, rstride=1, cstride=1, alpha=0.2,color='b',label='Detector focus autocol')
    plt.title(npa)
    ax.text(9, 0, 12, 'Infinite to focal: %0.3f actuator mm '%(diffInfAutoc), color='red')
    plt.show()
    return x,y,z,C,Zfocus * focus - diffInfAutoc, Zfocus * detfocus - diffInfAutoc



starsxF1 = [836.6,990.2,954.7]	
starsyF1 = [344.4,476.9,952.1]
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/FBGuider2018_NEW/OpenCluster/F1/'

x = np.array([1683,1854,1866,2036,1849,2021])
y = np.array([1260,1260,1040,1040,614,614])
z = np.array([13.143,13.081,13.118,13.068,13.131,13.118])
#plt.figure()
X,Y,Z, ax, Cdet = fit_quadratic_curve(x,y,z,n=100,order=1)
plt.show()
#
#
#
#starsxF1 = [836.6,990.2,954.7]	
#starsyF1 = [344.4,476.9,952.1]
##x,y,z,C = PlotFocus2DGuider(path=path, pa=119,sigma=True,order=2,starsx=starsxF1,starsy=starsyF1)
#
#x,y,z,C,focus1,focus2 = PlotFocus2DGuider(path=path,PlanDetector=Cdet, pa=119,sigma=True,order=2,starsx=starsxF1,starsy=starsyF1,Plot=True)



def wDetectorCenter(position=1559, wave = 206.20,dispersion=608/13, maskx=0):#disp on mu/nm donc pixel/nm
    print ('dispersion={}\n'.format(dispersion))
    position -=  1067
    
    center_detector = (2124-1067) /2
    print (position - center_detector)
    wavecenter = wave + (position - center_detector)/dispersion
    return wavecenter

#wDetectorCenter()
#PlotFocus2DGuider(path=path, pa=119,sigma=False)
#
#PlotFocus2DGuider(path=path, pa=-161)
#PlotFocus2DGuider(path=path, pa=-121)
#PlotFocus2DGuider(path=path, pa=159)

#    plt.figure()
#    for star in range (len(fwhmx)):#
#        plt.plot(arange(11),f(arange(11),*popt[star]))
#    plt.ylim(0,30)

#        
#    index[6]
#    plot(fwhms[6][0]);plot(fwhms[6][0][index[6]])
#    #ouvrir CSV
#    table = Table.read(path + 'stack8103288_pa+159_2018-06-11T06-19-19_table.csv',format='csv')
#    table1 = Table.read(path + 'stack8103235_pa+159_2018-06-11T06-19-03_table.csv',format='csv')
#    x,y = table['X_IMAGE'],table['Y_IMAGE']
#    x1,y1 = table1['X_IMAGE'],table1['Y_IMAGE']
#    plt.plot(x1,y1,'o')
#    plt.plot(x,y,'+')
#    
#dist,idx = crossmatch(np.array([x1,y1]).T,np.array([x,y]).T,max_distance=10)
#plt.plot(x,y,'x')
#plt.plot(x1[idx],y1[idx],'x')
#plt.plot(x1,x1[idx],'x')
#
#index[dist<10]
#matched = x1[idx]
#
#
#dist, idx = crossmatch(np.array(seen_stars['xcentroid', 'ycentroid']).view((float,2)) + np.array([dx,dy]),
#           np.array(F1_stars['Xguider', 'Yguider'][select]).view((float,2)))
#
##idx = [4,7,24,30,37]
#plt.scatter(F1_stars['Xguider'][matched], F1_stars['Yguider'][matched], 100, 
#            marker="D", color='b')
#
#F1_stars[matched]


#def ConvolveBoxPSF(x, amp=1,r1 = 40, r2 = 40, sigma2 = 40, offset = 0):
#    a = special.erf((r1-x)/np.sqrt(2*sigma2))
#    b = special.erf((r2+x)/np.sqrt(2*sigma2))
#    function = amp * ( a + b )/2*(r1+r2)
#    return offset + function

def ConvolveBoxPSF(x, amp=1, l=40, x0=0, sigma2=40, offset=0):
    a = special.erf((l - (x - x0))/np.sqrt(2*sigma2))
    b = special.erf((l + (x - x0))/np.sqrt(2*sigma2))
    function = amp * ( a + b )/4*l
    return offset + function


def ConvolveBoxPSF2(x, amp=1, l=40, x0=0, sigma2=40, offset=0, k=0):
    a = special.erf((l - (x - x0))/np.sqrt(2*sigma2))
    b = special.erf((l + (x - x0))/np.sqrt(2*sigma2))
    c = k*special.erf((l - (x - x0))/np.sqrt(2*sigma2/4))
    d = k*special.erf((l + (x - x0))/np.sqrt(2*sigma2/4))    
    function = amp * ( a + b - c - d)/4*l
    return offset + function


def ConvolveSlit2D_PSF2(xy, amp=1, l1 = 3, l2=3, L1=9, L2=9, sigma2 = 40):
    x, y = xy
    A1 = special.erf((l1-x)/np.sqrt(2*sigma2))
    A2 = special.erf((l2+x)/np.sqrt(2*sigma2)) 
    B1 = special.erf((L1-y)/np.sqrt(2*sigma2))
    B2 = special.erf((L2+y)/np.sqrt(2*sigma2)) 
    r1, r2 = (l1+l2)/2,(L1+L2)/2 
    function = amp * (1/(16*r1*r2)) * (A1+A2) * (B1 + B2)
    return function.ravel()




#x,y = np.linspace(-20,20,100), np.linspace(-20,20,101)
#xx,yy = np.meshgrid(x,y)
#
#conv = ConvolveSlit2D_PSF2((xx,yy),xo=20,yo=-20)
#plt.imshow(conv.reshape(101,100))
#
#plt.plot(conv.reshape(101,100)[50,:])
#plt.plot(conv.reshape(101,100)[:,50])

#r1,r2 = 25,10
#center = 0
#new = center + (r2-r1)/2
#print(new)
#plot(np.arange(-100,100),ConvolveBoxPSF(np.arange(-100,100),1,r1,r2,0))
#plot([new,new],[0,35])
def fit_quadratic_curve(x,y,z,n,sigma_z=None,order=2,Plot=True):
    if sigma_z is None:
        index = np.isfinite(z) 
#        data = np.array(zip(x[index],y[index],z[index]))  
        data = np.array([x[index],y[index],z[index]]).T
    else:
        index = (np.isfinite(z))  & (np.isfinite(sigma_z)) 
        data = np.array(zip(x[index],y[index],z[index]/sigma_z[index]))
  # regular grid covering the domain of the data
    X,Y = np.meshgrid(np.linspace(x.min(), x.max(), n), np.linspace(y.min(), y.max(), n))
    XX = X.flatten()
    YY = Y.flatten()
    
    order = order  # 1: linear, 2: quadratic
    if order == 1:
        # best-fit linear plane
        if sigma_z is None:
            A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        else:
            A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])] / sigma_z[:,np.newaxis]
            
        C,_,_,_ = linalg.lstsq(A, data[:,2])    # coefficients        
        # evaluate it on grid
        Z = C[0]*X + C[1]*Y + C[2]        
        # or expressed using matrix/vector product
        #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)    
    elif order == 2:
        if sigma_z is None:
        # best-fit quadratic curve
            A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
        else:
        # best-fit quadratic curve
            A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), (data[:,:2]**2)] / sigma_z[:,np.newaxis]
        C,_,_,_ = linalg.lstsq(A, data[:,2])        
        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)
    if Plot==True:
        fig = plt.figure(figsize=(15,10))#(10,8)
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2,color = 'r')
        ax.scatter(data[:,0], data[:,1], z[index], c='r', s=20)
        plt.xlabel('X')
        plt.ylabel('Y')
        ax.set_zlabel('Z')
        ax.axis('auto')
        ax.axis('tight')
    else:
        ax=1
    return X,Y,Z, ax, C

def quicklook(path,ptype='detector', mmr=None,t=1):
    image =[]
    if mmr is None:
        files = glob.glob(path+'*.fits')
    else:
        files = glob.glob(path+'*.fits')[mmr[0]:mmr[1]:mmr[2]]
    for i, file in enumerate(files):#[60:80:2]
        imagefits = fits.open(file)
        if ptype == 'guider':
            image.append(imagefits[0].data[:,:])
            plt.figure(figsize=(8,8))
        if ptype == 'detector':
            image.append(imagefits[0].data[:,1070:2130])
            plt.figure(figsize=(12,6))
        plt.title(imagefits.filename()[-26:-5])
        plt.imshow(image[i].T, vmin = np.percentile(image[i],t), vmax= np.percentile(image[i],100-t))
        plt.colorbar(orientation='horizontal');plt.show()
    return


def ConvolveDiskGaus2D(r, amp = 2 , RR = 4, sig = 4/2.35, offset=0):
    #RR=0.9*3.3/2 #4arcsec
    #RR=3.3/2 #4arcsec
    #RR=1.1*3.3/2 #4arcsec
    #RR=2*3.3/2 #4arcsec
    integrand =  lambda eta,r_ :  special.iv(0,r_ * eta / np.square(sig)) * eta * np.exp(-np.square(eta)/(2*np.square(sig)))
    #def integrand2 (eta, r_):
     #   return special.iv(0,r_ * eta / np.square(sig)) * eta * np.exp(-np.square(eta)/(2*np.square(sig)))
    integ = [quad(integrand,0,RR,args=(r_,))[0] * np.exp(-np.square(r_)/(2*np.square(sig))) / (np.pi*np.square(RR*sig)) for r_ in r]
    #integ = [np.exp(-np.square(r_)/(2*np.square(sig))) / (np.pi*np.square(RR*sig)) * np.nansum(integrand (np.linspace(0,RR,1000),r_))  for r_ in r]    
    #error = [quad(integrand,0,RR,args=(r_,))[1] * np.exp(-np.square(r_/(2*np.square(sig)))) / (np.pi*np.square(RR*sig)) for r_ in r]
    return offset + amp* np.array(integ)#, error



def stackOnCentroid(path, center, DS=0, radius=[20,20], size=100, all=False, function='mean', Plot=False):
    """
    Stack all images contained in a folder, if all=True all images contained in each sub folder
    """
    types = ('*.FIT', '*.fits', '*.fts','*.fit')
    files = []
    for type in types:
        files.extend(glob.glob(path+type))   

    n = len(files)
    image = fits.open(files[1])[0]
    lx, ly = image.data.shape
    stack = np.zeros((lx, ly, n))
    print('\nReading fits files...')
    for i,file in enumerate(files):
        try: 
            stack[:,:,i] = fits.open(file)[0].data.astype(np.int) - DS.astype(np.int)
        except ValueError:
            print("Problem of size, make sure all your images have the same dimension. \nA stacked image already maybe exists.")
            sys.exit()
    stackOnCentroid = np.zeros((2*size, 2*size, n))
    print('Re-centering images...')
    for i in range(n):
        extracted_image = stack[center[0]-radius[0]:center[0]+radius[0], center[1]-radius[1]:center[1]+radius[1], i]
#        imshow(extracted_image)
        x, y = np.where(extracted_image == extracted_image.max())
        x =  center[0] + x[0]-radius[0]
        y =  center[1] + y[0]-radius[0]
        extracted_image_centered = stack[x-size:x+size, y-size:y+size, i]
        if Plot:
            plt.imshow(extracted_image_centered)
            plt.grid()
            plt.show()
        stackOnCentroid[:,:,i] = extracted_image_centered
#        if i ==1:
#            imshow(fits.open(file)[0].data,vmin=np.percentile(fits.open(file)[0].data,1),vmax=np.percentile(fits.open(file)[0].data,99));colorbar(orientation="horizontal");title("Read noise");show()
        if function=='mean':
            image.data = np.nanmean(stackOnCentroid, axis=2)
        if function=='median':
            image.data = np.median(stackOnCentroid, axis=2)
    print('Images stacked')
    image.writeto( fits.open(files[0]).filename()[:-5] + 'stack.fits', overwrite=True)
    print('Stacked image save at: {}'.format(fits.open(files[0]).filename()[:-5] + 'stack.fits'))
    return image.data


def stackImages(path,all=False, DS=0, function = 'mean', numbers=None, save=True, name=""):
    """
    Stack all images contained in a folder, if all=True all images contained in each sub folder
    """
    exts = ('*.FIT', '*.fits', '*.fts','*.fit')
    files=[]
    if all == True:
        folders = os.walk(path).next()[1]#glob.glob(path+'/*')
        for path1 in folders:
            global stackImages
            stackImages(path+'/'+path1,all=False)
    else:
        if numbers is None:
            for ext in exts:
                files.extend(glob.glob(path + ext)) 
        else:
            print('Using files number specified')
            for i in numbers:
                for ext in exts:
                    files.extend(glob.glob("{}/image{:06d}{}".format(path, int(i), ext))) 
        print(print("\n".join(files)))
        n = len(files)
        image = fits.open(files[0])[0]
        lx,ly = image.data.shape
        stack = np.zeros((lx,ly,n))
        print('\nReading fits files...')
        for i,file in enumerate(files):
            with fits.open(file) as f:
                stack[:,:,i] = f[0].data
        if function=='mean':
            image.data = np.nanmean(stack,axis=2) - DS
        if function=='median':
            image.data = np.nanmedian(stack,axis=2) - DS
        print('Images stacked')
        if save:
            fname = path + '/'#os.path.splitext(files[0])[0][:-6] + '-'
            if 'NAXIS3' in image.header:
                image.header.remove('NAXIS3')             

            if numbers is None:
                name = fname + 'stack' + '-' + name + '.fits'                
                image.writeto(name ,overwrite=True)
                print('Stacked image save at: ' + name)
            else:
                name = '{}StackedImage_{}-{}-{}.fits'.format(fname, int(numbers[0]), int(numbers[-1]), name)
#                name = '{}StackedImage_{}-{}-{}.fits'.format(fname, numbers.min(), numbers.max(), name)
                image.writeto(name ,overwrite=True)
                fits.setval(name, 'DARKUSED', value = 0, comment = 'Images subtracted for dark subtraction')
                #add
                print('Stacked image save at: ' + name)                
        #plot = image.data[:,1070:2130]#[30:400,1070+700:2130]
        #plt.figure(figsize=(12,6))
        #plt.title( name)
        #plt.imshow(plot.T, vmin = np.percentile(plot,10), vmax= np.percentile(plot,90))
        #plt.colorbar(orientation='horizontal');plt.show()
    return image.data, name




def radial_profile_normalized(data, center, anisotrope=False, angle=30, radius=40, n=1.5, center_type='barycentre', radius_bg=70,n1=20, stddev=True, size=70):
    """Function that returns the radial profile of a spot
    given an input image + center.
    Use the azimuthal average to compute the profile and determine the encircled energy
    """
    y, x = np.indices((data.shape)) 
    print(data)
    #5#10  
    print('center_type = ',center_type)
    if center_type.lower() == 'maximum':
        image = data[int(center[1])-n1:int(center[1])+n1,int(center[0])-n1:int(center[0])+n1]
        barycentre =  np.array([np.where(image == image.max())[0][0], np.where(image == image.max())[1][0]])#ndimage.measurements.center_of_mass(data[center[1]-n1:center[1]+n1,center[0]-n1:center[0]+n1])
    if center_type.lower() == 'barycentre':
        image = data[int(center[1])-n1:int(center[1])+n1,int(center[0])-n1:int(center[0])+n1]
        #print ('Need to add background substraction')
        background = estimateBackground(data,center,radius,1.8 )
        new_image = image - background
        #print(new_image,new_image.shape)
        index = new_image > 0.5 * np.nanmax(new_image)#.max()
        #print(index)
        new_image[~index] = 0    
        barycentre = ndimage.measurements.center_of_mass(new_image)#background#np.nanmin(image)
    if center_type.lower() == 'user':
        barycentre = [n1,n1]
    else:
        image = data[int(center[1])-n1:int(center[1])+n1,int(center[0])-n1:int(center[0])+n1]
        print('Center type not understood, taking barycenter one')
        background = estimateBackground(data,center,radius,1.8 )
        new_image = image - background
        #print(new_image,new_image.shape)
        index = new_image > 0.5 * np.nanmax(new_image)#.max()
        #print(index)
        new_image[~index] = 0 
        barycentre = ndimage.measurements.center_of_mass(new_image)#background#np.nanmin(image)
    new_center = np.array(center) + barycentre[::-1] - n1
    print('new_center = {}, defined with center type: {}'.format(new_center, center_type))
    
    if radius_bg:
        fond = estimateBackground(data, new_center, radius, n)
    else:
        fond = 0
    image = data - fond#(data - fond).astype(np.int)
    
    r = np.sqrt((x - new_center[0])**2 + (y - new_center[1])**2)#    r = np.around(r)-1
    rint = r.astype(np.int)
    
    
    image_normalized = image #/ np.nansum(image[r<radius])
    if anisotrope == True:
        theta = abs(180*np.arctan((y - new_center[1]) / (x - new_center[0])) / np.pi)#    theta = np.abs(180*np.arctan2(x - new_center[0],y - new_center[1]) / np.pi)
        tbin_spectral = np.bincount(r[theta<angle].ravel(), image_normalized[theta<angle].ravel())
        tbin_spatial = np.bincount(r[theta>90-angle].ravel(), image_normalized[theta>90-angle].ravel())
        nr_spectral = np.bincount(r[theta<angle].ravel())
        nr_spatial = np.bincount(r[theta>90-angle].ravel())
        EE_spatial = 100 * np.nancumsum(tbin_spatial) / np.nanmax(np.nancumsum(tbin_spatial)[:100] + 1e-5)
        EE_spectral = 100 * np.nancumsum(tbin_spectral) / np.nanmax(np.nancumsum(tbin_spectral)[:100] + 1e-5)
        return tbin_spectral / nr_spectral, tbin_spatial / nr_spatial, EE_spectral, EE_spatial
    else:
        tbin = np.bincount(rint.ravel(), image_normalized.ravel())
        nr = np.bincount(rint.ravel())
        rsurf = np.sqrt(np.nancumsum(nr) /np.pi)
        rmean = np.bincount(rint.ravel(), r.ravel())/nr 
        if stddev:
            #import datetime
            #print(datetime.datetime.now())
            dist = np.array(rint[rint<radius].ravel(), dtype=int)
            data = image[rint<radius].ravel()
            stdd = [np.nanstd(data[dist==distance])/ np.sqrt(len(data[dist==distance])) for distance in np.arange(size)]
            #print(datetime.datetime.now())

        radialprofile = tbin / nr
        EE = np.nancumsum(tbin) * 100 / np.nanmax(np.nancumsum(tbin)[:radius] + 1e-5)
        return rsurf[:size], rmean[:size], radialprofile[:size], EE[:size], new_center[:size], stdd[:size]
    #rsurf, rmean, radialprofile, EE, new_center = radial_profile_normalized(fits.open('/Users/Vincent/Nextcloud/Work/FTS2018_FLIGHT/test/InstrumentCentering_180819/CALIBXY180824AF/TF2/stack15069515.fits')[0].data,center=[641,545])


def estimateBackground(data, center, radius=30, n=1.8):
    """Function that estimate the Background behing a source given an inner radius and a factor n to the outter radius
    such as the background is computed on the area which is on C2(n*radius)\C1(radius)
    """
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)
    mask = (r>=radius) & (r<=n*radius)
    fond = np.nanmean(data[mask])
    return fond




def Dirac2Slit(TFarray, TwoDprofile ,size=None):
    from DS9Utils import addAtPos
    TFarray = TFarray.astype(np.int)
    lx,ly = TFarray.shape
    NewSlit = np.zeros((lx,ly))
    cxs,cys = np.where(TFarray==1)
    for cx,cy in zip(cxs,cys):
        NewSlit = addAtPos(NewSlit, TwoDprofile, [cx,cy])
    return NewSlit
    
def radialEE(data, center, radius=30, n=3):
    """Function that only compute the encircled energy.
    Based on the same formula used to compute the radial profile, elypticity is not available
    """
    y, x = np.indices((data.shape))
    n1=5    
    image = data[center[1]-n1:center[1]+n1, center[0]-n1:center[0]+n1]
    barycentre =  np.array([np.where(image == image.max())[1][0],np.where(image == image.max())[0][0]])#ndimage.measurements.center_of_mass(data[center[0]-n1:center[0]+n1,center[1]-n1:center[1]+n1])
    barycentre =  ndimage.measurements.center_of_mass(data[center[1]-n1:center[1]+n1,center[0]-n1:center[0]+n1])
    new_center = np.array(center) + barycentre - n1
    print('new_center = {}'.format(new_center))
    r = np.sqrt((x - new_center[0])**2 + (y - new_center[1])**2)
    r = r.astype(np.int)
    data2 = data.copy()
    data2.ravel().sort()
    fond = estimateBackground(data, new_center, radius, n)
    tbin2 = np.bincount(r.ravel(), (data.ravel()-fond).astype(np.int))
    EE = np.cumsum(tbin2)#*100/np.nanmax(np.cumsum(tbin2)[:radius+20]+1e-5)
    return EE 



def plot_rp2_wo_latex(data, center, size=40, n=1.5, anisotrope=False, angle=30, radius=40, ptype='linear', fit=True, center_type='barycentre', maxplot=0.013, minplot=-1e-5, radius_ext=12, platescale=None):
  """Function used to plot the radial profile and the encircled energy of a spot,
  Latex is not necessary
  """
  if anisotrope == True:
      spectral, spatial, EE_spectral, EE_spatial = radial_profile_normalized(data, center, anisotrope=anisotrope, angle=angle, radius=radius, n=n, center_type=center_type)
      spectral = spectral[~np.isnan(spectral)]
      spatial = spatial[~np.isnan(spatial)]
      #min1 = min(spatial[:size])
      #min2 = min(spectral[:size])
      norm_spatial = spatial[:size]#(spatial[:n] - min(min1,min2)) / np.nansum((spatial[:n] - min(min1,min2) ))
      norm_spectral = spectral[:size]#(spectral[:n] - min(min1,min2)) / np.nansum((spectral[:n] - min(min1,min2) ))              
      if ptype == 'linear':
          popt1, pcov1 = curve_fit(gausexp, np.arange(size), norm_spatial)       
          popt2, pcov2 = curve_fit(gausexp, np.arange(size), norm_spectral) 
          plt.plot(np.arange(size), norm_spectral, label='spectral direction')   
          plt.plot(np.arange(size), norm_spatial, label='spatial direction') 
          if fit==True:
              plt.plot(np.linspace(0,size,10*size), gausexp(np.linspace(0,size,10*size), *popt1), label='Spatial fit')            
              plt.plot(np.linspace(0,size,10*size), gausexp(np.linspace(0,size,10*size), *popt2), label='Spectral fit')                       
              plt.figtext(0.5,0.5,'Sigma = %0.3f-%0.3f pix \nLambda = %0.3f-%0.3f pix \npcGaus = %0.0f-%0.3fpc' % (popt1[2],popt2[2],popt1[3],popt2[3],100*popt1[0]/(popt1[1] + popt1[0]),100*popt2[0]/(popt2[1] + popt2[0])), fontsize=11,bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
      else:
          plt.semilogy(np.arange(size),norm_spectral, label='spectral direction')   
          plt.semilogy(np.arange(size),norm_spatial, label='spatial direction')           
      return popt1, popt2
  else: 
          rsurf, rmean, a, EE, NewCenter = radial_profile_normalized(data, center, anisotrope=anisotrope, angle=angle, radius=radius, n=n, center_type=center_type)
          norm = a[:size]#(a[:n] - min(a[:n]) ) / np.nansum((a[:n] - min(a[:n]) ))
          if ptype == 'linear':
              fig, ax1 = plt.subplots(figsize=(12, 6))
              if platescale:
                  size = platescale * size
                  popt, pcov = curve_fit(gausexp,platescale * rmean[:len(norm)], norm, p0=[0.001,0.001, platescale*1, 1*platescale, 0.3])#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):    
                  ax1.set_xlabel('Distance to center in Diffusing Angle (DA) [arcsec]', fontsize=18)
              else:
                  popt, pcov = curve_fit(gausexp, rmean[:len(norm)], norm,p0=[0.001, 0.001, 1, 1, 1])#[1,1,1,1,1] 
                  ax1.set_xlabel('Distance to center [pix]', fontsize=18)                      

              ax1.plot(rmean[:len(norm)], norm, '+', c='black', label='Normalized isotropic profile')
              if fit==True: 
                  ax1.plot(np.linspace(0,size,10*size), gausexp(np.linspace(0, size, 10*size), *popt), c='royalblue') #)r"$\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!"
                  ax1.fill_between(np.linspace(0, size, len(norm)), norm - 1.5*np.abs(norm - gausexp(np.linspace(0, size, len(norm)), *popt)), 
                                   norm + 1.5*np.abs(norm - gausexp(np.linspace(0, size, len(norm)), *popt)), alpha=0.3, label=r"3*Residuals")
                  ax1.plot(np.linspace(0, size, 10*size), exp(np.linspace(0, size, 10*size), *popt), c='navy')
                  ax1.plot(np.linspace(0, size, 10*size), gaus(np.linspace(0, size, 10*size), *popt), c='blue')
                  ax1.set_ylabel('Radial Profile', color='b', fontsize=18)
                  ax1.tick_params('y', colors='b')
                  ax1.set_ylim((minplot, np.nanmax([np.nanmax(1.1*(norm)), maxplot])))
                  ax2 = ax1.twinx()
                  ax2.plot(rsurf[:len(norm)], EE[:len(norm)], 'r--x')
#                    print(np.linspace(0,size,len(norm)))
#                    mina = min(size[EE[:len(size)]>80])
                  mina = min(np.linspace(0,size, len(norm))[EE[:len(norm)]>77])
                  minb = min(np.linspace(0,size, len(norm))[EE[:len(norm)]>48])
#                    print(mina)
                  ax2.plot(np.linspace(minb, minb, 2), np.linspace(0, 50, 2), 'r-o')                    
                  ax2.plot(np.linspace(minb, size, 2), np.linspace(50, 50, 2), 'r-o')
                  ax2.plot(np.linspace(mina, mina, 2), np.linspace(0, 80, 2), 'r-o')                    
                  ax2.plot(np.linspace(mina, size, 2), np.linspace(80, 80, 2), 'r-o')
                  #EE_gaus = np.cumsum(gaus(np.linspace(0,size,100*size),*popt) *2 * np.pi * np.linspace(0,size,100*size)**1)
                  #EE_exp = np.cumsum(exp(np.linspace(0,size,100*size),*popt) * 2 * np.pi * np.linspace(0,size,100*size)**1)
                  ax2.set_ylabel('Encircled Energy', color='r', fontsize=18)
                  ax2.tick_params('y', colors='r')
                  fig.tight_layout()
                  ax1.xaxis.grid(True)
                  ax1.tick_params(axis='x', labelsize=18)
                  ax1.tick_params(axis='y', labelsize=18)
                  ax2.tick_params(axis='y', labelsize=18)                    
                  e_gaus = np.nansum(gaus(np.linspace(0,size,100*size),*popt) *2 * np.pi * np.linspace(0,size,100*size)**1)
                  e_exp = np.nansum(exp(np.linspace(0,size,100*size),*popt) * 2 * np.pi * np.linspace(0,size,100*size)**1)
                  ax1.legend(loc = (0.54,0.05),fontsize=18)
                  if platescale:
                      plt.figtext(0.63,0.57,"sigma = %0.3f DA[\"] \nLambda = %0.3f DA[\"] \npGaus = %0.2fpc" % (popt[2],popt[3],100*e_gaus/(e_gaus + e_exp)), 
                                  fontsize=19,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})#    norm_gaus = np.pi*sigma    norm_exp = 2*np.pi * lam**2 * gamma(2/alpha)/alpha    
                  else:
                      plt.figtext(0.63,0.57,"sigma = %0.3f pix \nLambda = %0.3f pix \npGaus = %0.2fpc" % (popt[2],popt[3],100*e_gaus/(e_gaus + e_exp)), 
                                  fontsize=19,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})#    norm_gaus = np.pi*sigma    norm_exp = 2*np.pi * lam**2 * gamma(2/alpha)/alpha

  #                plt.figtext(0.74,0.18,r'$\displaystyle\sigma =$ %0.3f pix \n$\displaystyle\lambda =$ %0.3f pix \n$\displaystyle pGaus = \%$%0.2f\n$\displaystyle\alpha = $%0.1f' % (popt[2],popt[3],100*e_gaus/(e_gaus + e_exp),popt[4]), fontsize=18,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})
  #                plt.show()
          else:
              plt.semilogy(np.arange(n), norm, label='isotropic profile')
  return popt, np.linspace(0, size, len(norm)), norm #np.array([np.linspace(0,size,len(norm)),norm])


def plot_rp2_convolved_wo_latex(data, center, size=40, n=1.5, anisotrope=False, angle=30, radius=40, ptype='linear', fit=True, center_type='barycentre', maxplot=0.013, minplot=-1e-5, radius_ext=12, platescale=None,fibersize = 100,SigmaMax=4):
  """Function used to plot the radial profile and the encircled energy of a spot,
  Latex is not necessary
  """
  if anisotrope == True:
      spectral, spatial, EE_spectral, EE_spatial = radial_profile_normalized(data, center, anisotrope=anisotrope, angle=angle, radius=radius, n=n, center_type=center_type)
      spectral = spectral[~np.isnan(spectral)]
      spatial = spatial[~np.isnan(spatial)]
      #min1 = min(spatial[:size])
      #min2 = min(spectral[:size])
      norm_spatial = spatial[:size]#(spatial[:n] - min(min1,min2)) / np.nansum((spatial[:n] - min(min1,min2) ))
      norm_spectral = spectral[:size]#(spectral[:n] - min(min1,min2)) / np.nansum((spectral[:n] - min(min1,min2) ))              
      if ptype == 'linear':
          popt1, pcov1 = curve_fit(gausexp, np.arange(size), norm_spatial)       
          popt2, pcov2 = curve_fit(gausexp, np.arange(size), norm_spectral) 
          plt.plot(np.arange(size), norm_spectral, label='spectral direction')   
          plt.plot(np.arange(size), norm_spatial, label='spatial direction') 
          if fit==True:
              plt.plot(np.linspace(0,size,10*size), gausexp(np.linspace(0,size,10*size), *popt1), label='Spatial fit')            
              plt.plot(np.linspace(0,size,10*size), gausexp(np.linspace(0,size,10*size), *popt2), label='Spectral fit')                       
              plt.figtext(0.5,0.5,'Sigma = %0.3f-%0.3f pix \nLambda = %0.3f-%0.3f pix \npcGaus = %0.0f-%0.3fpc' % (popt1[2],popt2[2],popt1[3],popt2[3],100*popt1[0]/(popt1[1] + popt1[0]),100*popt2[0]/(popt2[1] + popt2[0])), fontsize=11,bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
      else:
          plt.semilogy(np.arange(size),norm_spectral, label='spectral direction')   
          plt.semilogy(np.arange(size),norm_spatial, label='spatial direction')           
      return popt1, popt2
  else: 
          rsurf, rmean, profile, EE, NewCenter = radial_profile_normalized(data, center, anisotrope=anisotrope, angle=angle, radius=radius, n=n, center_type=center_type)
          profile = profile[:size]#(a[:n] - min(a[:n]) ) / np.nansum((a[:n] - min(a[:n]) ))
          if ptype == 'linear':
              fig, ax1 = plt.subplots(figsize=(8, 4))
              if platescale:
                  size = platescale * size
                  #popt, pcov = curve_fit(ConvolveDiskGaus2D, np.linspace(0,size,size), profile, p0=[2, platescale*2, 2*platescale])#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):    
                  popt, pcov = curve_fit(ConvolveDiskGaus2D, platescale * rmean[:size], profile, p0=[2, platescale*2, 2*platescale, np.nanmean(profile)])#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):    
                  ax1.set_xlabel('Distance to center in Diffusing Angle (DA) [arcsec]', fontsize=18)
              else:
                  #popt, pcov = curve_fit(ConvolveDiskGaus2D, np.linspace(0,size,size), profile, p0=[2,2,2])#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):  3.85  
                  fiber = fibersize #/ (2*1.08*(1/0.083))
                  popt, pcov = curve_fit(ConvolveDiskGaus2D, rmean[:size], profile, p0=[1,fiber,2, np.nanmean(profile)],bounds=([0,0.95*fiber-1e-5,1,-1],[2,1.05*fiber+1e-5,SigmaMax,1]))#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):    
                  ax1.set_xlabel('Distance to center [pix]', fontsize=18)                      
              #ax1.plot(rmean[:size], profile, '+', c='black', label='Normalized isotropic profile')
              ax1.plot(rmean[:size], profile, '+', c='black', label='Normalized isotropic profile')
              #ax1.plot(np.linspace(0, size, size), profile, '+', c='black', label='Normalized isotropic profile')
              if fit==True: 
                  ax1.plot(np.linspace(0,size,10*size), ConvolveDiskGaus2D(np.linspace(0, size, 10*size), *popt), c='royalblue') #)r"$\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!"
                  ax1.fill_between(rmean[:size], profile - 1.5*np.abs(profile - ConvolveDiskGaus2D(rmean[:size], *popt)), 
                                   profile + 1.5*np.abs(profile - ConvolveDiskGaus2D(rmean[:size], *popt)), alpha=0.3, label=r"3*Residuals")
                  #ax1.plot(np.linspace(0, size, 10*size), exp(np.linspace(0, size, 10*size), *popt), c='navy')
                  #ax1.plot(np.linspace(0, size, 10*size), gaus(np.linspace(0, size, 10*size), *popt), c='blue')
                  ax1.set_ylabel('Radial Profile', color='b', fontsize=18)
                  ax1.tick_params('y', colors='b')
                  ax1.set_ylim((minplot, np.nanmax([np.nanmax(1.1*(profile)), maxplot])))
                  ax2 = ax1.twinx()
                  #ax2.plot(np.linspace(0, size, len(norm)), EE[:len(norm)], 'r--x')
                  EE_interp = interpolate.interp1d(rsurf[:size], EE[:size],kind='cubic')
                  ninterp = 10
                  xnew = np.linspace(rsurf[:size].min(),rsurf[:size].max(),ninterp*len(rsurf[:size]))
                  ax2.plot(xnew,EE_interp(xnew),linestyle='dotted',c='r')
#                  x = np.linspace(0,size,100*size)
#                  aire = np.pi * np.square(np.linspace(0,size,100*size))
#                  rp2 = ConvolveDiskGaus2D(np.linspace(0, size, 100*size), *popt) - ConvolveDiskGaus2D(np.linspace(0, size, 100*size), *popt).min()
 #                 ee = np.cumsum(aire*rp)
#                  ax2.plot(x,100*(ee/ee.max()),linestyle='dotted')
                  
                  ax2.plot(rsurf[:size], EE[:size], 'rx')
#                    print(np.linspace(0,size,len(norm)))

                  mina = min(xnew[EE_interp(xnew)[:ninterp*size]>79])
                  minb = min(xnew[EE_interp(xnew)[:ninterp*size]>49])

#                    print(mina)
                  ax2.plot(np.linspace(minb, minb, 2), np.linspace(0, 50, 2), 'r-o')                    
                  ax2.plot(np.linspace(minb, size, 2), np.linspace(50, 50, 2), 'r-o')
                  ax2.plot(np.linspace(mina, mina, 2), np.linspace(0, 80, 2), 'r-o')                    
                  ax2.plot(np.linspace(mina, size, 2), np.linspace(80, 80, 2), 'r-o')
                  #EE_gaus = np.cumsum(gaus(np.linspace(0,size,100*size),*popt) *2 * np.pi * np.linspace(0,size,100*size)**1)
                  #EE_exp = np.cumsum(exp(np.linspace(0,size,100*size),*popt) * 2 * np.pi * np.linspace(0,size,100*size)**1)
                  ax2.set_ylim((0, 110))
                  ax2.set_ylabel('Encircled Energy', color='r', fontsize=18)
                  ax2.tick_params('y', colors='r')
                  fig.tight_layout()
                  ax1.xaxis.grid(True)
                  ax1.tick_params(axis='x', labelsize=18)
                  ax1.tick_params(axis='y', labelsize=18)
                  ax2.tick_params(axis='y', labelsize=18)                    
#                  e_gaus = np.nansum(gaus(np.linspace(0,size,100*size),*popt) *2 * np.pi * np.linspace(0,size,100*size)**1)
#                  e_exp = np.nansum(exp(np.linspace(0,size,100*size),*popt) * 2 * np.pi * np.linspace(0,size,100*size)**1)
                  ax1.legend(loc = (0.54,0.05),fontsize=18)
                  if platescale:
                      plt.figtext(0.53,0.53,"Amp = %0.3f\nRadius = %0.3f pix \nSigmaPSF = %0.3f pix \nEE50-80 = %0.2f - %0.2f p" % (popt[0],platescale*popt[1],platescale*popt[2],platescale*minb,platescale*mina), 
                                  fontsize=19,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})#    norm_gaus = np.pi*sigma    norm_exp = 2*np.pi * lam**2 * gamma(2/alpha)/alpha
                  else:
                      plt.figtext(0.53,0.53,"Amp = %0.3f\nRadius = %0.3f pix \nSigmaPSF = %0.3f pix \nEE50-80 = %0.2f - %0.2f p" % (popt[0],popt[1],popt[2],minb,mina), 
                                  fontsize=19,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})#    norm_gaus = np.pi*sigma    norm_exp = 2*np.pi * lam**2 * gamma(2/alpha)/alpha
                  d = {"SizeSource":popt[1],"FWHM":popt[2],"EE50":mina,"EE80":minb,"Platescale":platescale,"Center":NewCenter}
                  print("SizeSource = {}\nFWHM = {} \nEE50 = {}\nEE80 = {}\nPlatescale = {}\nCenter = {}".format(popt[1],popt[2],minb,mina,platescale,NewCenter))
                  return d
  #                plt.figtext(0.74,0.18,r'$\displaystyle\sigma =$ %0.3f pix \n$\displaystyle\lambda =$ %0.3f pix \n$\displaystyle pGaus = \%$%0.2f\n$\displaystyle\alpha = $%0.1f' % (popt[2],popt[3],100*e_gaus/(e_gaus + e_exp),popt[4]), fontsize=18,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})
  #                plt.show()
          else:
              plt.semilogy(rmean[:size], profile, label='isotropic profile')
  return #[2]#popt, np.linspace(0, size, len(norm)), norm #np.array([np.linspace(0,size,len(norm)),norm])



def AnalyzeSpot(data, center, size=40, n=1.5,radius=40, fit=True, center_type='barycentre', radius_ext=12, platescale=None,fibersize = 100,SigmaMax=4):
  """Function used to plot the radial profile and the encircled energy of a spot,
  Latex is not necessary
  """

  rsurf, rmean, profile, EE, NewCenter = radial_profile_normalized(data, center, radius=radius, n=n, center_type=center_type)
  profile = profile[:size]#(a[:n] - min(a[:n]) ) / np.nansum((a[:n] - min(a[:n]) ))
  #popt, pcov = curve_fit(ConvolveDiskGaus2D, np.linspace(0,size,size), profile, p0=[2,2,2])#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):  3.85  
  fiber = fibersize / (2*1.08*(1/0.083))
  if fiber == 0:
      gaus = lambda x, a, sigma: a**2 * np.exp(-np.square(x / sigma) / 2)
      popt, pcov = curve_fit(gaus, rmean[:size], profile, p0=[1, 2])#,bounds=([0,0],[1,5]))#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):    
  else:
      popt, pcov = curve_fit(ConvolveDiskGaus2D, rmean[:size], profile, p0=[1,fiber,2, np.nanmean(profile)],bounds=([0,0.95*fiber-1e-5,1,-1],[2,1.05*fiber+1e-5,SigmaMax,1]))#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):      
  EE_interp = interpolate.interp1d(rsurf[:size], EE[:size],kind='cubic')
  ninterp = 10
  xnew = np.linspace(rsurf[:size].min(),rsurf[:size].max(),ninterp*len(rsurf[:size]))
  mina = min(xnew[EE_interp(xnew)[:ninterp*size]>79])
  minb = min(xnew[EE_interp(xnew)[:ninterp*size]>49])
  if fiber == 0:
      flux = 2*np.pi*np.square(popt[1])*np.square(popt[0])
      d = {"Flux":flux,"SizeSource":0,"Sigma":abs(popt[1]),"EE50":mina,"EE80":minb,"Platescale":platescale,"Center":NewCenter}
      print("Flux = {}\nSizeSource = {}\nSigma = {} \nEE50 = {}\nEE80 = {}\nPlatescale = {}\nCenter = {}".format(flux,0,popt[1],minb,mina,platescale,NewCenter))
  else:
      d = {"Flux":0,"SizeSource":popt[1],"Sigma":abs(popt[2]),"EE50":mina,"EE80":minb,"Platescale":platescale,"Center":NewCenter}
      print("Flux = 0\nSizeSource = {}\nSigma = {} \nEE50 = {}\nEE80 = {}\nPlatescale = {}\nCenter = {}".format(popt[1],popt[2],minb,mina,platescale,NewCenter))

#  if fiber == 0:
#      d = {"SizeSource":0,"Sigma":popt[1],"EE50":mina,"EE80":minb,"Platescale":platescale,"Center":NewCenter}
#      print("SizeSource = {}\nFWHM = {} \nEE50 = {}\nEE80 = {}\nPlatescale = {}\nCenter = {}".format(0,popt[1],minb,mina,platescale,NewCenter))
#  else:
#      d = {"SizeSource":popt[1],"Sigma":popt[2],"EE50":mina,"EE80":minb,"Platescale":platescale,"Center":NewCenter}
#      print("SizeSource = {}\nFWHM = {} \nEE50 = {}\nEE80 = {}\nPlatescale = {}\nCenter = {}".format(popt[1],popt[2],minb,mina,platescale,NewCenter))
  return d


def plot_rp2(data, center, size=40, n=1.5, anisotrope=False, angle=30, radius=40, ptype='linear', fit=True, center_type=None, maxplot=0.013, minplot=-1e-5, radius_ext=12, platescale=None):
  """Function used to plot the radial profile and the encircled energy of a spot,
  Latex is necessary, used plot_rp2_latex if you dont have latexds9
  """
  if anisotrope == True:
      spectral, spatial, EE_spectral, EE_spatial = radial_profile_normalized(data, center, anisotrope=anisotrope, angle=angle, radius=radius, n=n, center_type=center_type)
      spectral = spectral[~np.isnan(spectral)]
      spatial = spatial[~np.isnan(spatial)]
      #min1 = min(spatial[:size])
      #min2 = min(spectral[:size])
      norm_spatial = spatial[:size]#(spatial[:n] - min(min1,min2)) / np.nansum((spatial[:n] - min(min1,min2) ))
      norm_spectral = spectral[:size]#(spectral[:n] - min(min1,min2)) / np.nansum((spectral[:n] - min(min1,min2) ))              
      if ptype == 'linear':
          popt1, pcov1 = curve_fit(gausexp, np.arange(size), norm_spatial)       
          popt2, pcov2 = curve_fit(gausexp, np.arange(size), norm_spectral) 
          plt.plot(np.arange(size), norm_spectral, label='spectral direction')   
          plt.plot(np.arange(size), norm_spatial, label='spatial direction') 
          if fit==True:
              plt.plot(np.linspace(0, size, 10*size), gausexp(np.linspace(0, size, 10*size), *popt1), label='Spatial fit')            
              plt.plot(np.linspace(0, size, 10*size), gausexp(np.linspace(0, size, 10*size), *popt2), label='Spectral fit')                       
              plt.figtext(0.5,0.5,'Sigma = %0.3f-%0.3f pix \nLambda = %0.3f-%0.3f pix \npcGaus = %0.0f-%0.3fpc' % (popt1[2],popt2[2],popt1[3],popt2[3],100*popt1[0]/(popt1[1] + popt1[0]),100*popt2[0]/(popt2[1] + popt2[0])), fontsize=11,bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
      else:
          plt.semilogy(np.arange(size), norm_spectral, label='spectral direction')   
          plt.semilogy(np.arange(size), norm_spatial, label='spatial direction')           
      return popt1, popt2
  else: 
      rsurf, rmean, a, EE = radial_profile_normalized(data, center, anisotrope=anisotrope, angle=angle, radius=radius, n=n, center_type=center_type)
      norm = a[:size]#(a[:n] - min(a[:n]) ) / np.nansum((a[:n] - min(a[:n]) ))
      if ptype == 'linear':
          fig, ax1 = plt.subplots(figsize=(12, 6))
          if platescale:
              size = platescale * size
              popt, pcov = curve_fit(gausexp, platescale * rmean[:len(norm)], norm, p0=[0.001, 0.001, platescale*1, 1*platescale, 0.3])#[1,1,1,1,1] (x,a,b,sigma,lam,alpha):      
              ax1.set_xlabel(r'\textbf{Distance to center in Diffusing Angle (DA) [arcsec]}', fontsize=18)
          else:
              popt, pcov = curve_fit(gausexp, rmean[:len(norm)], norm, p0=[0.001, 0.001, 1, 1, 1])#[1,1,1,1,1] 
              ax1.set_xlabel(r'\textbf{Distance to center [pix]}', fontsize=18)                      
          plt.rc('text', usetex=True)
          plt.rc('font', family='serif')
          ax1.plot(rmean[:len(norm)], norm, '+', c='black', label=r'Normalized isotropic profile')
          if fit==True: 
              ax1.plot(np.linspace(0, size, 10*size), gausexp(np.linspace(0, size, 10*size), *popt), c='royalblue',label=r"Fit:$\displaystyle f(r) = g(r) + h(r) $") #)r"$\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!"
              ax1.fill_between(np.linspace(0, size, len(norm)), norm - 1.5*np.abs(norm-gausexp(np.linspace(0, size, len(norm)), *popt)), norm + 1.5*np.abs(norm-gausexp(np.linspace(0, size, len(norm)), *popt)), alpha=0.3, label=r"3*Residuals")
#                ax1.plot(np.linspace(0,size,10*size),exp(np.linspace(0,size,10*size),*popt), c='navy', label = r"$\displaystyle g(r) = \frac{A_1 \alpha}{2 \pi \lambda^2 \Gamma(2/\alpha)} e^{- (r/\rho)^\alpha}$")
#                ax1.plot(np.linspace(0,size,10*size),gaus(np.linspace(0,size,10*size),*popt), c='blue', label = r"$\displaystyle h(r) = \frac{A_2}{2 \pi \sigma^2} e^{-r^2/\sigma}$")
              ax1.plot(np.linspace(0, size,10*size), exp(np.linspace(0, size, 10*size), *popt), c='navy', label=r"$\displaystyle g(r) = A_1 e^{- (r/\rho)}$")
              ax1.plot(np.linspace(0, size,10*size), gaus(np.linspace(0, size, 10*size), *popt), c='blue', label=r"$\displaystyle h(r) = A_2 e^{-r^2/\sigma}$")
              ax1.set_ylabel(r'\textbf{Radial Profile}', color='b', fontsize=18)
              ax1.tick_params('y', colors='b')
              ax1.set_ylim((minplot,np.nanmax([np.nanmax(1.1*(norm)),maxplot])))
              ax2 = ax1.twinx()
              ax2.plot(rsurf[:len(norm)], EE[:len(norm)], 'r--x')
              mina = min(np.linspace(0,size,len(norm))[EE[:len(norm)]>77])
              minb = min(np.linspace(0,size,len(norm))[EE[:len(norm)]>48])
#                    print(mina)
              ax2.plot(np.linspace(minb,minb,2), np.linspace(0,50,2), 'r-o')                    
              ax2.plot(np.linspace(minb,size,2), np.linspace(50,50,2), 'r-o')
              ax2.plot(np.linspace(mina,mina,2), np.linspace(0,80,2), 'r-o')                    
              ax2.plot(np.linspace(mina,size,2), np.linspace(80,80,2), 'r-o')
              #EE_gaus = np.cumsum(gaus(np.linspace(0,size,100*size),*popt) *2 * np.pi * np.linspace(0,size,100*size)**1)
              #EE_exp = np.cumsum(exp(np.linspace(0,size,100*size),*popt) * 2 * np.pi * np.linspace(0,size,100*size)**1)
              ax2.set_ylabel(r'\textbf{Encircled Energy}', color='r',fontsize=18)
              ax2.tick_params('y', colors='r')
#                fig.tight_layout()
              ax1.xaxis.grid(True)
              ax1.tick_params(axis='x', labelsize=18)
              ax1.tick_params(axis='y', labelsize=18)
              ax2.tick_params(axis='y', labelsize=18)
#                plt.grid()

              e_gaus = np.nansum(gaus(np.linspace(0,size,100*size),*popt) *2 * np.pi * np.linspace(0,size,100*size)**1)
              e_exp = np.nansum(exp(np.linspace(0,size,100*size),*popt) * 2 * np.pi * np.linspace(0,size,100*size)**1)
              ax1.legend(loc = (0.59,0.05),fontsize=18)
#                align_yaxis_np(ax1, ax2)
#                plt.figtext(0.74,0.60,"sigma = %0.3f pix \nLambda = %0.3f pix \npGaus = %0.2fpc\nalpha = %0.1f\nA1 = %0.4f\nA2 = %0.4f" % (popt[2],popt[3],100*e_gaus/(e_gaus + e_exp),popt[4]**2, 2*np.pi*popt[2]**2 * popt[0]**2,popt[1]**2* 2*np.pi * popt[3]**2 * gamma(2/popt[4]**2)/popt[4]**2), fontsize=19,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})#    norm_gaus = np.pi*sigma    norm_exp = 2*np.pi * lam**2 * gamma(2/alpha)/alpha
              if platescale:
                  plt.figtext(0.63,0.62,"sigma = %0.3f DA[\"] \nLambda = %0.3f DA[\"] \npGaus = %0.2fpc" % (popt[2],popt[3],100*e_gaus/(e_gaus + e_exp)), fontsize=19,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})#    norm_gaus = np.pi*sigma    norm_exp = 2*np.pi * lam**2 * gamma(2/alpha)/alpha    
              else:
                  plt.figtext(0.68,0.63,"sigma = %0.3f pix \nLambda = %0.3f pix \npGaus = %0.2fpc" % (popt[2],popt[3],100*e_gaus/(e_gaus + e_exp)), fontsize=19,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})#    norm_gaus = np.pi*sigma    norm_exp = 2*np.pi * lam**2 * gamma(2/alpha)/alpha

#                plt.figtext(0.74,0.18,r'$\displaystyle\sigma =$ %0.3f pix \n$\displaystyle\lambda =$ %0.3f pix \n$\displaystyle pGaus = \%$%0.2f\n$\displaystyle\alpha = $%0.1f' % (popt[2],popt[3],100*e_gaus/(e_gaus + e_exp),popt[4]), fontsize=18,bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})
#                plt.show()
      else:
          plt.semilogy(np.arange(n),norm, label = 'isotropic profile')


def gaus(x, a, b, sigma, lam, alpha):
    """1D gaussian centered on zero
    """
    gaus = a**2 * np.exp(-np.square(x / sigma) / 2) 
    return gaus 


def exp(x, a, b, sigma, lam, alpha):
    """1D exponential centered on zero
    """
    exp =  b**2 * np.exp(-(x/lam)**(1**2))
    return exp


def gausexp(x, a, b, sigma, lam, alpha):
    """1D gaussian + exponential centered on zero
    """
    gaus = a**2 * np.exp(-np.square(x / sigma) / 2) 
    exp =  b**2 * np.exp(-(x/lam)**(1**2))
    return gaus + exp  


def through_slit(path='Detector/170924/', nimages=np.arange(2,15), pos_image=np.arange(2,15), radius=15, center=[933, 1450], n_bg=1.3, sizefig=4):#, center_bg=[500,500]
    print('Sum pixel is used (another estimator may be prefarable)')
    fluxes=[]
    n=radius
    for i in nimages:
        image = fits.open(path + 'image%06d.fits'%(i))[0].data
        plt.figure(figsize=(sizefig,sizefig))
        plt.imshow(image[center[0]-n:center[0]+n, center[1]-n:center[1]+n])#;plt.colorbar();plt.show()
        #flux = np.nansum(image[center[0]-n:center[0]+n,center[1]-n:center[1]+n])-np.nansum(image[center_bg[0]-n:center_bg[0]+n,center_bg[1]-n:center_bg[1]+n])
        flux = np.nanmean(image[center[0]-n:center[0]+n, center[1]-n:center[1]+n]) - estimateBackground(image, center, radius, n_bg)
        fluxes.append(flux)
        fluxesn = (fluxes - min(fluxes)) / max(fluxes - min(fluxes))
        maxf = pos_image[np.where(fluxes==np.nanmax(fluxes))][0]
    plt.figure()
    plt.plot(pos_image, fluxesn,'--*')
    plt.plot(np.linspace(maxf, maxf, len(fluxes)), fluxesn/max(fluxesn))
    plt.grid()
    plt.xlabel('# image')
    plt.ylabel('Sum pixel')
    #plt.savefig(path + 'image%06d.fits'%(i))
    plt.plot()
    return pos_image, fluxesn




#a = create_DS9regions([[1],[2]],[[1],[2]],save=True,color=['red','green'],form='square')#['red','green']
#create_DS9regions(xim=[10,10], yim=[20,10], radius=20, save=True, savename="test", form=['square','square'], color=['green','red'],ID=[3,4])
#create_DS9regions([F2.table['X_IMAGE']], [F2.table['Y_IMAGE']], radius=20, save=True, savename="test", form=['circle'], color=['green'],ID=[F2.internalCount])
#create_DS9regions([F2.xdetpix], [F2.ydetpix], radius=20, save=True, savename="test", form=['circle'], color=['green'],ID=[F2.internalCount])

def plot3d(data, center, n=50):
    """
    Plot a spot in 3D, thresholds should be modified
    """
    data = data[center[1]-n:center[1]+n, center[0]-n:center[0]+n]
    data[data>21000] = 13800
    data[data<13500] = 13800
    n1=2
    image = ndimage.grey_dilation(ndimage.grey_erosion(data, size=(n1,n1)), size=(n1,n1))  
    xx, yy = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(xx, yy, image ,rstride=1, cstride=1, cmap='cubehelix',
            linewidth=0)
    image2 = image.copy()
    image2 = image2.flatten()
    image2.sort()
    print(image2)
    ax.set_zlim(13500, image2[-10]+100)
    plt.show() 

def ComputeSmallRotationOffset(delta,position):
    data = np.concatenate((delta[:,0], delta[:,1]))
    row_x = np.hstack((-position[:,[1]], np.ones((len(position),1)), np.zeros((len(position),1)) )) # xn - x =  x dgamma - y theta + dx
    row_y = np.hstack((position[:,[0]], np.zeros((len(position),1)), np.ones((len(position),1)) )) # yn - y =  y dgamma + x theta      + dy
    mat = np.vstack((row_x, row_y))
    matinv =  np.linalg.pinv(mat)
    sol = matinv.dot(data)
    theta_rad = sol[0]
    deltax = sol[1]
    deltay = sol[2]
    theta = theta_rad*180/np.pi*60 #arcmin
    print("theta: {} arcmin\ndx: {} pixel\ndy: {} pixel".format(theta, deltax, deltay))
    covar = matinv.dot(matinv.T)
    # accuracy, assuming 1 arcsec measurement error
    print("variances: {}\n".format(np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 3600, 3600])) #
    #residual
    data_new = mat.dot(sol)
    print("residuals in arcsec:", (data_new - data).reshape((len(position),2))*3600)
    return theta_rad, deltax, deltay

def GaussianM(x, x0, amp):
    xo = np.array(x0)
    sigma = 2.5
    amplitude = amp
    #A = amplitude/(sigma * np.sqrt(2*np.pi))    
    g = 8600 + amplitude * np.exp( - 0.5*(((x-xo)/sigma)**2))
    return g







    
def delete_doublons(sources, dist):
    """Function that delete doublons detected in a table, 
    the initial table and the minimal distance must be specifies
    """
    try:
        sources['doublons'] = 0
        for i in range(len(sources)):
            a = distance(sources[sources['doublons']==0]['xcentroid'],sources[sources['doublons']==0]['ycentroid'],sources['xcentroid'][i],sources['ycentroid'][i]) > dist
            a = list(1*a)
            a.remove(0)
            if np.nanmean(a)<1:
                sources['doublons'][i]=1
        return sources[sources['doublons']==0]
    except TypeError:
        print('no source detected')
#        quit()
        

def distance(x1,y1,x2,y2):
    """
    Compute distance between 2 points in an euclidian 2D space
    """
    return np.sqrt(np.square(x1-x2)+np.square(y1-y2))


def distance2(a,b):
    """
    Compute distance between 2 points in an euclidian 2D space
    """
    a = np.array(a)
    b = np.array(b)
    return np.sqrt(np.square((a - b)).sum())
    #return np.sqrt(np.square(x1-x2)+np.square(y1-y2))


def Gaussian(x, amplitude, xo, sigma2, offset):
    """Defines a gaussian function with offset
    """
    xo = float(xo)
    #A = amplitude/np.sqrt(2 * np.pi * sigma2)   #attention is not used anymore  
    g = offset + amplitude * np.exp( - 0.5*(np.square(x-xo)/sigma2))
    return g.ravel()





def Normalization(amplitude, xo, sigma):
    """Given the parameters of a gaussian function, returns the amplitude so that it is normalized in area 
    """
    B = dblquad(lambda x: Gaussian((x,), amplitude, xo, np.square(sigma)), -float('inf'), +float('inf'))       
    return (amplitude/B[0], xo, sigma)

        


def interpDispersion(x,y,disp, order=2):
    """Interpolate and map the distortion
    """
    data = np.array([x,y,disp]).T
    a = 1200
    b = 2000#900
    c = 100
    d = 1800#1600wavelength
    # regular grid covering the domain of the data
    X,Y = np.meshgrid(np.arange(a,b,10), np.arange(c,d, 10))
    X,Y = np.meshgrid(np.linspace(a,b,7), np.linspace(c,d,20))
    XX = X.flatten()
    YY = Y.flatten()
    
    if order == 1:
        # best-fit linear plane
        A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        C,_,_,_ = linalg.lstsq(A, data[:,2])    # coefficients
        
        # evaluate it on grid
        Z = C[0]*X + C[1]*Y + C[2]
        
        # or expressed using matrix/vector product
        #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)
    
    elif order == 2:
        # best-fit quadratic curve
        A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
        C,_,_,_ = linalg.lstsq(A, data[:,2])
        
        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)
    
    # plot points and fitted surface
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
    ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlabel('Z')
    ax.axis('equal')
    ax.axis('tight')
    plt.show()
    
    plt.figure(figsize=(10,8))
    plt.title('Mask to det dispersion measured with zinc lamp')
#    plt.xlabel('Detector [mm]' )
#    plt.ylabel('Detector [mm]')
    plt.xlabel('Field of view x-coordinate (arcmin)' )
    plt.ylabel('Field of view y-coordinate (arcmin)')
    plt.imshow(13*Z.T,cmap='jet')#coolwarm        
#    plt.yticks(np.linspace(0,77,5),[-450/60,-225/60,0,225/60,450/60])
#    plt.xticks(np.linspace(0,160,5),[-780/60,-390/60,0,390/60,780/60])
    plt.yticks(np.linspace(0,7-1,5),[-450/60,-225/60,0,225/60,450/60])
    plt.xticks(np.linspace(0,20-1,5),[-780/60,-390/60,0,390/60,780/60])
    cb = plt.colorbar(orientation='horizontal')
    cb.set_label("dispersion [microns/nanometer]")
    #plt.savefig('/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/InstrumentCaracteriscs/dispersion.png')
    return 13*Z
#z = interpDispersion(np.concatenate((x2,x3,x4)),np.concatenate((y2,y3,y4)),np.concatenate((disp2,disp3,disp4)))



def defoc_hyperboloid(xy, cx, cy, a, b):
    x, y = xy
    best2 = np.square(x-cx) + np.square(y-cy)
    
    return b*np.sqrt(1+best2/(a*a))
    

def defoc_paraboloid(xy, cx, cy, a, b):
    x,y = xy
    best = np.sqrt(np.square(x-cx) + np.square(y-cy))
    
    return (a*best)**2 + b


def defoc_parabola(xy, offset, dx, a, b):
    x, y = xy
    best = y + offset + dx*x
    
    return (a*best)**2 + b
    

def make_model_allw(model):
    
    def model_allw(xywi, c0, c1, a0, b0, a1, b1, a2, b2):
        x, y, wi = xywi
        wi = wi.astype(int)
        a = wi.choose([a0, a1, a2]  )
        b = wi.choose([b0, b1, b2]  )
        return model((x,y), c0, c1, a, b)
        
    return model_allw


def best_focus_line(x, offset, dx):
    
    return -(offset + x*dx)


def plot_parabola(focus, w, xy, zdata, model, popt=None, allwave=False):

    fitdone = popt is not None
    
    x, y = xy
    
    if fitdone:
        gx, gy = np.meshgrid(np.linspace(-12, 12, 20), np.linspace(-6, 6, 20))
        if allwave:
            gx_det = np.zeros(gx.shape + (3,))
            gy_det = gx_det.copy()
            gz = gx_det.copy()
            for i,ww in enumerate(focus.w):
                gxy_det = focus.apply_map_mask_detector(gx, gy, np.full_like(gx, ww))
                gx_det[...,i],  gy_det[...,i] = gxy_det[0], gxy_det[1]
                model_allw = make_model_allw(model)
                gz[...,i] = model_allw((gx_det[...,i], gy_det[...,i], np.full(gx_det.shape[:-1], i)), *popt)
        else:
            gxy_det = focus.apply_map_mask_detector(gx, gy, np.full_like(gx, w))
            gx_det, gy_det = gxy_det[0], gxy_det[1]
            gz = model((gx_det, gy_det), *popt)

    if allwave:
        colors = [focus.wcolors[ww] for ww in w]
    else:
        colors = focus.wcolors[w]
        
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(121, projection='3d')
    if allwave:
        ax.set_title('defocus jointly at all wavelength')
    else:        
        ax.set_title('defocus lambda={:.4f}'.format(w))
    ax.set_xlabel('x detector pix')
    ax.set_ylabel('y detector pix')
    ax.set_zlabel('FWHM size pix')
    ax.scatter(x, y, zdata, c=colors)
    ylim = [0, 2069]
    xlim = [focus._DFLT_Yinf, focus._DFLT_Yinf + 1064]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if fitdone:
        if allwave:
            for i, c in enumerate(['blue','green', 'red']):  
                imax = np.where(gx_det[...,i] >= xlim[0])[0].max()
                imin = np.where(gx_det[...,i] <= xlim[1])[0].min()
                ax.plot_wireframe(gx_det[imin:imax,:,i],
                                  gy_det[imin:imax,:,i],
                                  gz[imin:imax,:,i], colors=c, alpha=.15, clip_on=True)
        else:
            imax = np.where(gx_det >= xlim[0])[0].max()
            imin = np.where(gx_det <= xlim[1])[0].min()
            ax.plot_wireframe(gx_det[imin:imax, :],
                              gy_det[imin:imax, :],
                              gz[imin:imax, :], colors=focus.wcolors[w], alpha=.15, clip_on=True)
        if model == defoc_parabola:                           
            ax.set_zlim([0,30])
        else:
            ax.set_zlim([0,10])
    
    # plot on flat scatter
    zmin = max(zdata.min(), 2)
    zmax = min(zdata.max(), 30)
    marker_size = lambda z: 20 + (z-zmin)/(zmax - zmin)*300

    ax = fig.add_subplot(122)
    if allwave:
        ax.set_title('defocus jointly at all wavelength')
    else:  
        ax.set_title('defocus lambda={:.4f}'.format(w))
    ax.set_xlabel('x detector pix')
    ax.set_ylabel('y detector pix')
    print("zdata minmax : {} {}".format(zdata.min(), zdata.max()))
    ax.scatter(x, y, s=marker_size(zdata), c='none', edgecolors=colors)
    #ref_mask = (np.array([-4.6, -3.94, -2.2]), np.array([3, 6, 9]))
    ref_in_mask = (np.array([ -4.50677685,  -4.09257524,  -2.13289232]), np.array([3, 6, 9])) 
    # comes from still pen measurement where tilt mask is crossing the mean science mask plane
    ref_in_mask = focus.pa_rotation(ref_in_mask, focus.pa, reverse=False)
    if model == defoc_parabola:
        if allwave:
            for ww, c in zip(focus.w, ['blue','green', 'red']):  
                ref_det = focus.apply_map_mask_detector(ref_in_mask[0], ref_in_mask[1], np.array([ww, ww, ww]))
                ax.scatter(ref_det[0], ref_det[1], marker='*', c=c, s=100 )
        else:
            ref_det = focus.apply_map_mask_detector(ref_in_mask[0], ref_in_mask[1], np.array([w, w, w]))
            print("ref_det")
            print(ref_det)
            ax.scatter(ref_det[0], ref_det[1], marker='*', c='black', s=100 )
    #ax.plot(gx_det[imin:imax,10], np.full(imax-imin,ref[1]), ':', c=focus.wcolors[w])
    #ax.plot( np.full(20, ref[0]), gy_det[10,:], ':', c=focus.wcolors[w])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

        
    if fitdone:
        # best line plot
        if allwave:
            for i, c in enumerate(['blue','green', 'red']):  
                imax = np.where(gx_det[..., i] >= xlim[0])[0].max()
                imin = np.where(gx_det[..., i] <= xlim[1])[0].min()
                ax.scatter(gx_det[imin:imax, :, i], gy_det[imin:imax, :, i], s=marker_size(gz[imin:imax,:,i]), c=c, alpha=.15)
            if model == defoc_parabola:
                ax.plot(xlim, best_focus_line(np.array(xlim), popt[0], popt[1]), c='black')
            elif model == defoc_paraboloid or model == defoc_hyperboloid:
                ax.scatter(popt[0], popt[1], s=100, c='black')
        else:
            imax = np.where(gx_det >= xlim[0])[0].max()
            imin = np.where(gx_det <= xlim[1])[0].min()
            ax.scatter(gx_det[imin:imax,:], gy_det[imin:imax,:], s=marker_size(gz[imin:imax,:]), c=colors, alpha=.15)
            if model == defoc_parabola:
                ax.plot(gx_det[imin:imax,10], best_focus_line(gx_det[imin:imax,10], popt[0], popt[1]), c=colors)
            elif model == defoc_paraboloid or model == defoc_hyperboloid:
                ax.scatter(popt[0], popt[1], s=100, c='black')

        
def select_focus_crit(table, fwhm_type='quad'):
    if fwhm_type == 'quad':
        zdata = np.sqrt(np.square(table['FWHM_x']) + np.square(table['FWHM_y']))
        sigma = np.sqrt(table['FWHM_var_x'] + table['FWHM_var_y'])
    elif fwhm_type == 'max':
        zdata = np.zeros(len(table))
        sigma = np.zeros(len(table))
        xgty = table['FWHM_x'] > table['FWHM_y']
        zdata[xgty] = table[xgty]['FWHM_x']
        sigma[xgty] = np.sqrt(table[xgty]['FWHM_var_x'])
        zdata[~xgty] = table[~xgty]['FWHM_y']
        sigma[~xgty] = np.sqrt(table[~xgty]['FWHM_var_y'])
        #zdata = np.maximum(table['FWHM_x'], table['FWHM_y'])
    elif fwhm_type == 'fwhmx':
        zdata = table['FWHM_x']
        sigma = np.sqrt(table['FWHM_var_x'])
    elif fwhm_type == 'fwhmy':
        zdata = table['FWHM_y']
        sigma = np.sqrt(table['FWHM_var_y'])
        
    return zdata, sigma
 


def estim_focus_tilted_mask(focus, tilt=2.56,  model=defoc_parabola, 
                            use_mapped_det=False, plot=False, 
                            allwave=False, fwhm_threshold=None, fwhm_type='quad'):
    '''
    Estim detector focus from tilted mask by fitting parabolic cylinder
    
        estim_focus_tilted_mask(focus, tilt=4., xmask0=0., 
                                use_mapped_det=False, plot=False)
        
            focus: is an instance of Focus object (image with spots detected 
                                                   and measured)
            
            tilt: angle of mask in degree (default to 2.56.)
            
            xmask0: position on the mask (in mm) where the tilted mask cross 
                    the focal sphere (default to 0)
            
            
            use_mapped_det: boolean, if True use mapped coordinates of the spots
                            on the detector, else use measured coordinates. 
                            Note: the mapping is the measured one
                            (default False)
            
            plot: boolean (default False), if True do some plots
    '''
    
    # get centers
    try:
        centers = focus.centers
        mapping = focus.mapping
    except AttributeError:
        mapping, centers = focus.map_mask_detector(bywave=True, plot=False)
    
   
    # do it by wavelength
    for w, m in zip(focus.w, mapping):
        subt = focus.table[focus.table['wavelength'] == w]

        if len(subt) < 10:
            print("Not enough spots at wavelength {}".format(w))
            continue
        
        # remove points where FWHM not defined
        subt = subt[(subt['FWHM_x'] > 0) & (subt['FWHM_y'] > 0)]
        if fwhm_threshold is not None:
            subt = subt[(subt['FWHM_x'] < fwhm_threshold) & (subt['FWHM_y'] < fwhm_threshold)]
            
        sid = subt['id_slit']
        zdata, sigma = select_focus_crit(subt, fwhm_type)
               
        if use_mapped_det:
            xy_det = focus.apply_map_mask_detector(focus.x_mask[sid], 
                                                   focus.y_mask[sid], 
                                                   np.full_like(focus.x_mask[sid], w))
            x_det, y_det = xy_det[0], xy_det[1]
        else:
            x_det, y_det = subt['X_IMAGE'], subt['Y_IMAGE'] 
                
        #xydata = np.array([focus.xmask[sid], x_det])
        xydata = (x_det, y_det)
        
        if model == defoc_parabola:
            p0 = [-1000., 0., 1E-1, 2]
        elif model == defoc_paraboloid or model == defoc_hyperboloid:
            p0 = [1500., 1000., 1E-1, 2] # field center

        #bounds = (np.array([-np.inf, -np.inf, 0, 0]), np.inf)
        try:
            popt, pcov = curve_fit(model, xydata, zdata, sigma=sigma , p0=p0, maxfev=5000)#, bounds=bounds)#, maxfev=10000)
        except RuntimeError or ValueError:
            print("Fit at wavelength {} not achieved".format(w))
            popt = None
            print(RuntimeError)
        else:
            if model == defoc_parabola:
                yoffset, xslope, sqrta, b = popt
                yoffset_err, xslope_err, sqrta_err,  b_err = np.sqrt(pcov.diagonal())
            else:
                xcenter, ycenter, sqrta, b = popt
                xcenter_err, ycenter_err, sqrta_err,  b_err = np.sqrt(pcov.diagonal())
                
            if model == defoc_parabola:
                print('''
                detector focus fitted at wavelength {}:
                      y offset: {}  +-{}
                      x slope:  {}  +-{}
                      a: {}  +-{}
                      b: {}  +-{}
                      '''.format(w, yoffset, yoffset_err, xslope, xslope_err, 
                                  sqrta**2, 2*sqrta*sqrta_err, b, b_err))
            else:
                print('''
                detector focus fitted at wavelength {}:
                      x center: {}  +-{}
                      y center: {}  +-{}
                      a: {}  +-{}
                      b: {}  +-{}
                      '''.format(w, xcenter, xcenter_err, ycenter, ycenter_err, 
                                  sqrta**2, 2*sqrta*sqrta_err, b, b_err))
                    
                          
        if plot:
             plot_parabola(focus, w, xydata, zdata, model, popt)
                
                
    # do it for all  wavelength
    if allwave:
        subt = focus.table[focus.table['wavelength'] > 0 ]
    
        # remove points where FWHM not defined
        subt = subt[(subt['FWHM_x'] > 0) & (subt['FWHM_y'] > 0)]
        if fwhm_threshold is not None:
            subt = subt[(subt['FWHM_x'] < fwhm_threshold) & (subt['FWHM_y'] < fwhm_threshold)]
        sid = subt['id_slit']
        w = subt[ 'wavelength']
        wi = np.searchsorted(focus.zinc_lines, w)        

        zdata, sigma = select_focus_crit(subt, fwhm_type)
        
        if use_mapped_det:
            xy_det = focus.apply_map_mask_detector(focus.x_mask[sid], 
                                                   focus.y_mask[sid], 
                                                   subt['wavelength'])        
            x_det, y_det = xy_det[0], xy_det[1]
        else:
            x_det, y_det = subt['X_IMAGE'], subt['Y_IMAGE'] 
                
        xydata = (x_det, y_det, wi)
        
        model_allw = make_model_allw(model)
        if model == defoc_parabola:
            p0_allw = p0 + 2*[1E-1, 2]
        else:
            p0_allw = p0 + 2*[1E-1, 2]
        
        #bounds = (np.array([-np.inf, -np.inf, 0., 0., 0., 0., 0., 0.]), np.inf)                        
        try:
            popt, pcov = curve_fit(model_allw, xydata, zdata, sigma=sigma, p0=p0_allw, maxfev=5000)#, bounds = bounds)# maxfev=5000)
        except RuntimeError or ValueError:
            print("Fit at all wavelength not achieved")
            popt = None
            print(RuntimeError)
        else:
            if model == defoc_parabola:
                yoffset, xslope, sqrta0, b0, sqrta1, b1, sqrta2, b2 = popt
                yoffset_err, xslope_err, sqrta0_err, b0_err, sqrta1_err,  b1_err , sqrta2_err,  b2_err = np.sqrt(pcov.diagonal())
            else:
                xcenter, ycenter, sqrta0, b0, sqrta1, b1, sqrta2, b2 = popt
                xcenter_err, ycenter_err, sqrta0_err, b0_err, sqrta1_err,  b1_err , sqrta2_err,  b2_err = np.sqrt(pcov.diagonal())
        
            if model == defoc_parabola:
                print('''
                detector focus fitted jointly at all wavelength :
                      y offset: {}  +-{}
                      x slope: {}  +-{}
                      a0: {}  +-{}
                      b0: {}  +-{}
                      a1: {}  +-{}
                      b1: {}  +-{}
                      a2: {}  +-{}
                      b2: {}  +-{}
                      '''.format(yoffset, yoffset_err, xslope, xslope_err, 
                              sqrta0**2, 2*sqrta0*sqrta0_err, b0, b0_err, 
                              sqrta1**2, 2*sqrta1*sqrta1_err, b1, b1_err, 
                              sqrta2**2, 2*sqrta2*sqrta2_err, b2, b2_err))
            else:
                print('''
                detector focus fitted jointly at all wavelength :
                      x center: {}  +-{}
                      y center: {}  +-{}
                      a0: {}  +-{}
                      b0: {}  +-{}
                      a1: {}  +-{}
                      b1: {}  +-{}
                      a2: {}  +-{}
                      b2: {}  +-{}
                      '''.format(xcenter, xcenter_err, ycenter, ycenter_err, 
                                  sqrta0**2, 2*sqrta0*sqrta0_err, b0, b0_err, 
                                  sqrta1**2, 2*sqrta1*sqrta1_err, b1, b1_err, 
                                  sqrta2**2, 2*sqrta2*sqrta2_err, b2, b2_err))
                            

        if plot:
             plot_parabola(focus, w, xydata[0:2], zdata, model, popt, allwave=True)


    # convert fitted value to usable quantities
        
    return
     
#a = Focus(filename='/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/Autocoll/Detector/1/image000379.fits',plot=True,HumanSupervision=True)
#a = Focus(filename='/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/Autocoll/Detector/1/image000379.fits',plot=True)
         
    
class Focus(object):

    _DFLT_extension  = '.fits'
    _DFLT_catalog  = 'cat'
    _DFLT_computeImage  = False

#------------------------------------------------------------------------------#
# Special Methods                                                              #
#------------------------------------------------------------------------------#                                                 
    def __init__(self, **kwargs):
        """
        Z calibration (through focus analysis)
        This code helps analyse the focus (using detector images or guider images).
        It includes:
        	- Spot detection in the whole field (with or without human supervision - plot=True and press enter to accept the spot, ¬´ n ¬ª to reject it)
        	- Fitting a gaussian on every detected spot 
        	- Plotting the images + spots + evolution of the FWHM with respect to the field
        The code is easy to use, you can only enter as parameter the file name of an image but there are also other parameters you can use.
        One of the most important to use is the ¬´ quick ¬ª one. As we observe some defocus on the images, the detection code (based on astropy DAO star finder) we have to iterate on the gaussian kernel and the detection threshold we use to detect the maximum of sources with a different focus. quick=False will allow this iteration.
        This codes returns different outputs:
        -A nice image that shows how the image quality change in the field (attached)
        -A csv file with information about all thesubt[['xcentroid','ycentroid']] spot detected (position, FWHM, flux, EE50%, etc)
        -A region (.reg) file that can be open with DS9 that contains a circle that encircled each PSFat its EE80%
        """
        print('''\n\n\n\n      START IMAGE ANALYSIS \n\n\n\n''')
        print('''reading parameters...''')
        
        image = 'image0000{}'
        number = kwargs.get('number', None)
        Type = kwargs.get('Type', 'guider')
        self.reversex = kwargs.get('reversex', False)

        self.centre = kwargs.get('centre', None)
        print('centre = {}'.format(self.centre))
        subDark = kwargs.get('subDark', False)
        self.line = kwargs.get('line', None)
        verbose = kwargs.get('verbose', False)#     
        quick = kwargs.get('quick', True)
        verbose = kwargs.get('verbose', False)
        HumanSupervision = kwargs.get('HumanSupervision', False)        
        plot = kwargs.get('plot', True)
        shape = kwargs.get('shape', 'holes')

        fname = kwargs.get('filename', None)# size = kwargs.get('size', [self._DFLT_Xinf,self._DFLT_Xsup,self._DFLT_Yinf,self._DFLT_Ysup])
        plot = kwargs.get('plot', False)
        print('plot =', plot)
        stack = kwargs.get('stack', True)
        stack_image = kwargs.get('stack_image', False)
        py = kwargs.get('py', True)
        figsize = kwargs.get('figsize', 12)
        
        #cloudpath = kwargs.get('cloudpath', '/Users/Vincent/Nextcloud/FIREBALL/')
        try:
            Target_dir = resource_filename('__package__', 'Targets')
        except:
            Target_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Targets')
            
        self.windowing = kwargs.get('windowing', False) 
        self.mask = kwargs.get('mask', 'grid')
        pa = {'f1': 119., 'f2':-161. , 'f3':-121. , 'f4':159., 'grid':-81., 'tilted':-1 }            
        self.pa0 = pa[self.mask.lower()]
        self.pa = kwargs.get('pa', self.pa0)
        if self.windowing:
            if self.mask.lower() == 'f1':
                self.slit_location = os.path.join(Target_dir, 'targets_F1.txt')
                offset = [1044-1573,25-999] #[1044-1573,25-999]
            if self.mask.lower() == 'grid':
                self.slit_location = os.path.join(Target_dir, 'grid_mask.txt')
                offset = [-490,-905]  #disp-field
            if self.mask.lower() == 'f2':
                self.slit_location = os.path.join(Target_dir, 'targets_F2.txt')
                offset = [-490-13,-905+20]  #disp-field
            if self.mask.lower() == 'f3':
                self.slit_location = os.path.join(Target_dir, 'targets_F3.txt')
                offset = [-490-9,-905+4]  #disp-field
            if self.mask.lower() == 'f4':
                self.slit_location = os.path.join(Target_dir, 'targets_F4.txt')
                offset = [-490-14,-905-81]  #disp-field                        
            self.date = kwargs.get('date', 29)
            print('date given : ', self.date)
            manual_offsets = kwargs.get('manual_offsets', [0,0])            
            if (self.reversex == True) or (self.date>10):
                xo, yo = 0, 0
                x1, y1 = 0, 0
                if self.date < 26: 
                    xo,yo = 7,58
                elif self.date > 26: # 27 = image cold 0601
                    xo,yo = -70+40,+30#-103#-6+32+70,-28-30#-76,50
                if self.date > 27: # 28 = images cold 0606
                    x1,y1 = 4-7,-23
                if self.date == 180612: # images cold 180612
                    offset_180612 = {'F1': (x1+10,y1+6), 'F2': (x1+6,y1+8), 'F3':(x1+5,y1+9), 'F4':(x1+9,y1+9), 'grid':(x1,y1) }
                    x1, y1 = offset_180612[self.mask]
                elif self.date>28:  # images 1808
                    x1,y1 = 4-7+20,-23+6 #testvincent +20 et +6 ajouté le 23 aout 2018
                    x1, y1 = x1+93, y1-93# after flight, because of soem detector shift due to impact most probably
                print('x1,y1', x1,y1)
                if self.mask == 'F1':
                    self.slit_location = os.path.join(Target_dir, 'targets_F1.txt')
                    offset = [-629-xo+x1+manual_offsets[0],-896-yo+y1+manual_offsets[1]] #[1044-1573,25-999]
                if self.mask == 'grid':
                    print('Using manual offsets : ', manual_offsets)
                    manual_offsets = kwargs.get('manual_offsets', [-3,-47]) 
                    self.slit_location = os.path.join(Target_dir, 'grid_mask.txt')
                    offset = [-599-xo+x1+manual_offsets[0],-819-yo-103+y1+manual_offsets[1]] # dispersion - pour aller a droite  ,spatiale + pour descendre,       
                if self.mask == 'F2':
                    self.slit_location = os.path.join(Target_dir, 'targets_F2.txt')
                    offset = [-599-xo+x1+manual_offsets[0],-807-yo+y1+manual_offsets[1]]  #disp-field
                if self.mask == 'F3':
                    self.slit_location = os.path.join(Target_dir, 'targets_F3.txt')
                    offset = [-599-xo+x1+manual_offsets[0],-806-yo+y1+manual_offsets[1]]  #disp-field
                if self.mask == 'F4':
                    manual_offsets = kwargs.get('manual_offsets', [0,100]) 
                    print('Using manual offsets : ', manual_offsets)
                    self.slit_location = os.path.join(Target_dir, 'targets_F4.txt')
                    offset = [-607-xo+x1+manual_offsets[0],-912-yo-103+y1+manual_offsets[1]]  #disp-field  
                if self.mask == 'tilted':
                    self.slit_location = os.path.join(Target_dir, 'targets_Tilted.txt')
                    offset = [-490-13-116+16+4+89-xo+x1+manual_offsets[0],-905+20+28+50+26-yo+y1+manual_offsets[1]]  # like F2
                    

            
            self.source = kwargs.get('source', 'Zn') # can be 'Zn', 'Lya' or 'all'
            IMOversion = kwargs.get('IMOversion', 2015)
            try:
                Mappings_dir = resource_filename('DS9FireBall', 'Mappings')
            except:
                Mappings_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Mappings')

            
            if IMOversion == 2015:
                self.sky2mask_mapping = os.path.join(Mappings_dir,'mapping-sky2mask-Baseline-April-2015.pkl')
                self.sky2det_mapping = os.path.join(Mappings_dir,'mapping-sky2det-Baseline-April-2015.pkl')
            if IMOversion == 2018:
                self.sky2mask_mapping = os.path.join(Mappings_dir,'mapping-sky2mask-Baseline-April-2018.pkl')
                self.sky2det_mapping = os.path.join(Mappings_dir,'mapping-sky2det-Baseline-April-2018.pkl')

        
        self.all = fits.open(fname)
        all = fits.open(fname)[0]
        if self.reversex:
            self.data = all.data[:,::-1]
            self.all[0].data = self.data
            self.all.writeto(fname[:-5]+'reversed.fits',overwrite=True)
        else:
            self.data = all.data
        self.header = all.header

        if self.header['BITPIX']==-32:
            Type='guider'
            self._DFLT_Xinf   = kwargs.get('Xinf',0)
            self._DFLT_Xsup   = kwargs.get('Xsup',1080)  
            self._DFLT_Yinf   = kwargs.get('Yinf',0)
            self._DFLT_Ysup   = kwargs.get('Ysup',1280)            
#        if self.header['BITPIX']==-64:
        else:
            Type='detector'
            self._DFLT_Xinf   = kwargs.get('Xinf',0)
            self._DFLT_Xsup   = kwargs.get('Xsup',2069)  
            self._DFLT_Yinf   = kwargs.get('Yinf',1064)
            self._DFLT_Ysup   = kwargs.get('Ysup',2130) 

        if quick:    
            if Type == 'guider':
                threshold = kwargs.get('threshold', 10)
                fwhm = kwargs.get('fwhm', 8)
            if Type == 'detector':
                threshold = kwargs.get('threshold', 13)
                fwhm = kwargs.get('fwhm', 5)
        else: 
            threshold = kwargs.get('threshold', [15])
            fwhm = kwargs.get('fwhm', [12.5,15])

        print(('To extract the sources: fwhm = {}, threshold = {}'.format(fwhm, threshold)))
            
        ndark = 0#77
        
        if stack_image:                            
            self.image = (self.data.astype('int') + 
                          fits.open( 'image%06d.fits'%(number-1))[0].data.astype('int') + 
                          fits.open('image%06d.fits'%(number+1))[0].data.astype('int') + 
                          fits.open( 'image%06d.fits'%(number+2))[0].data.astype('int') + 
                          fits.open( 'image%06d.fits'%(number-2))[0].data.astype('int'))/5
            self.all[0].data = self.image    
            try:
#                fits.delval(fname,'NAXIS3')
                self.all.writeto(fname[:-5]+'stacked' +self._DFLT_extension, overwrite = True)
            except KeyError:
                fits.delval(fname,'NAXIS3')
                self.all.writeto(fname+'stacked' +self._DFLT_extension, overwrite = True)
        else:
            if (Type == 'detector') or (Type == 'guider'):
                print('''\n\n\n\n     Re-sizing image : [0:%i, 0:%i] --> [%i:%i, %i:%i] \n\n\n\n'''%(self.data.shape[0],self.data.shape[1],self._DFLT_Xinf,self._DFLT_Xsup,self._DFLT_Yinf,self._DFLT_Ysup)) 
                self.image = self.data[self._DFLT_Xinf:self._DFLT_Xsup,self._DFLT_Yinf:self._DFLT_Ysup]#.T  #[1000:1100,1064:2130]  
            #else:        
                #self.image = self.data
        #n=0
        self.data = self.image 
        
        if subDark:
            self.dark = (fits.open( image.format(ndark) + self._DFLT_extension)[0].data + 
                         fits.open( image.format(ndark-2)+ self._DFLT_extension)[0].data + 
                         fits.open( image.format(ndark-1)+ self._DFLT_extension)[0].data + 
                         fits.open( image.format(ndark+1) + self._DFLT_extension)[0].data + 
                         fits.open( image.format(ndark+2)+ self._DFLT_extension)[0].data)/5 #[size[0]:size[1],size[2]:size[3]]
            self.image_dark = self.image - self.dark
            self.all[0].data = self.image_dark
            self.all.writeto(image.format(number)+'stack-dark' +self._DFLT_extension, overwrite = True)
            self.image =  self.image#[400:,1064:1542]  
            self.dark = self.dark#[400:,1064:1542]                         
            self.image = self.image-self.dark


        self.dispersion  = 0
        self.spatialResolution  = 0    
        self.spectralResolution  = 0
        self.pix2sec = 0.93
        self.disp = 0.2167#A/px    # 17.8 A/mm

        self.zinc_lines = np.array([0.20255, 0.20619, 0.21382])
        self.lya_line  = 1215.67*1e-4
        self.wcolors = {}
        for w, c in zip(self.zinc_lines, ['blue', 'green', 'red']):
            self.wcolors[w] = c
        self.wcolors[self.lya_line] = 'white'
        self.wcolors[-1.0] = 'white'
        self.wcolors[0.0] = 'white'

        if self.windowing:
            manual_offsets
            print('''reading parameters...''')
#            from .DS9Utils import returnXY
#            x1, y1, redshift1, slit1, w1 = returnXY(line='202'+self.mask)
#            x2, y2, redshift2, slit2, w2 = returnXY(line='206'+self.mask)
#            x3, y3, redshift3, slit3, w3 = returnXY(line='213'+self.mask)
#            xdetpix, ydetpix, w = np.hstack((x1, x2, x3)), np.hstack((y1, y2, y3)), np.hstack((w1, w2, w3))
            xdetmm, ydetmm, wtab = self.computeLinesLocation()
            xdetpix, ydetpix = self.putOffsetMaskDetect(xdetmm, ydetmm, offset)

            ScienceMasks = ['F1','F2','F3','F4']
 
            if (self.source == 'lya') and (self.mask not in ScienceMasks):
                raise Exception("lya source not compatible with calib mask")
                        
            colors = [self.wcolors[i] for i in self.w]
            self.ID =[self.internalCount]*len(self.w) #[self.internalCount for i in range(len(self.w))] 
            xmaskw =[self.xmask]*len(self.w) #[self.xmask for i in range(len(self.w))] 
            ymaskw =[self.xmask]*len(self.w) #[self.ymask for i in range(len(self.w))] 
            # create_ds9_regions(xdetpix[1], ydetpix[1], form='bpanda', radius=15, save=True, savename=self.all.filename()[:-5]+'predicted', color='blue')

            xdetpix, ydetpix, wtab = np.broadcast_arrays(xdetpix, ydetpix, wtab)
            t = []
            #for x, y, w, count, xmask, ymask in zip(xdetpix, ydetpix, wtab, self.ID, xmaskw, ymaskw):#np.array(self.internalCount)*len(self.w)):
            self.z =[self.z]*len(self.w)
            for x, y, w, count, xmask, ymask , redshift in zip(xdetpix, ydetpix, wtab, self.ID, xmaskw, ymaskw,list(self.z)):#np.array(self.internalCount)*len(self.w)):
                #print (x, y, w,count)
                t.append(Table([x, y, w, np.arange(self.xmask.size), count, xmask, ymask, redshift], names=('x', 'y', 'w','id_slit','Internal-count', 'xmask', 'ymask','redshift')))
            self.predictedSpots = np.hstack(t)
            tab = Table(self.predictedSpots,names=('x','y','w','id_slits','Internal-count','xmask','ymask','redshift'))
            tab.write(self.all.filename()[:-5] + '_table_predicted' + '.csv', format='csv',overwrite = True)
            #print(tab) 
        if py == True:
            print('use of py')
            peak_threshold = kwargs.get('peak_threshold', None)
            self.extractSourcesPy(stack=stack, plot=plot, verbose=verbose, quick=quick, 
                                  threshold=threshold, fwhm=fwhm, centre=self.centre, Type=Type,
                                  peak_threshold=peak_threshold,shape=shape)
            if self.abort:
                print( 'Aborting function, no source detected!')
                return
        else:
            if verbose:
                print('Use of sextractor')
            # self.extractCatalog(**kwargs)            
            self.tranform_catalog(**kwargs)            
            print('Number of sources detected and kept: ', len(self.table))

#        self.Fit2D( **kwargs)
        
        if (shape == 'holes') or (shape == 'slits'):
            print('We fit the profile by the convolution of a box by a gaussian'   )             
            print('Fitting in spatial direction')
            min_fwhm = kwargs.get('min_fwhm', 2.)
            self.popt2, self.pcov2 = self.Fit1D(HumanSupervision, stack, verbose, Type, shape=shape, nstack=1, min_fwhm=min_fwhm, convolved=True, axis=1)
            #self.popt2 = self.Fit1D2_Convolved(stack=stack, HumanSupervision=HumanSupervision, verbose=verbose, quick=quick, Type=Type, shape=shape, min_fwhm=min_fwhm)[0]
            print('Fitting in spectral direction')
            self.popt1, self.pcov1 = self.Fit1D(HumanSupervision, stack, verbose, Type, shape=shape, nstack=3, min_fwhm=min_fwhm, convolved=True, axis=0)
            #self.popt1 = self.Fit1D1_Convolved(stack=stack, HumanSupervision=HumanSupervision, verbose=verbose, quick=quick, Type=Type, shape=shape, min_fwhm=min_fwhm)[0]
        if shape == 'gaussian':
            print('We fit the profile by a gaussian')
            print('Fitting in spectral direction')
            min_fwhm = kwargs.get('min_fwhm', 2.)
            self.popt1, self.pcov1 = self.Fit1D(HumanSupervision, stack, verbose, Type, shape=shape, nstack=3, min_fwhm=min_fwhm, convolved=False, axis=0)
            #self.popt1 = self.Fit1D1(stack = stack, plot = plot, verbose = verbose, quick = quick, Type=Type)[0]
            print('Fitting in spatial direction')
            self.popt2, self.pcov2 = self.Fit1D(HumanSupervision, stack, verbose, Type, shape=shape, nstack=3, min_fwhm=min_fwhm, convolved=False, axis=1)
            #self.popt2 = self.Fit1D2(stack = stack, plot = plot, verbose = verbose, quick = quick, Type=Type)[0]
        if shape == 'slits2D':
            self.popt3 = self.Fit2D_Convolved(stack = stack, plot = plot, verbose = verbose, quick = quick)[0]
            
        self.table['wavelength'] = 0.0
        self.table['id_slit'] = -1
        self.table['redshift'] = -1.0
##########Difference Python DS9 0-1 ########
#        self.table['X_IMAGE'] += 1.0
#        self.table['Y_IMAGE'] += 1.0
#        self.table['xcentroid'] += 1.0
#        self.table['ycentroid'] += 1.0
##########Difference Python DS9 0-1 ########
        
        if self.windowing:
            #print(self.predictedSpots)
            self.dist_min = kwargs.get('dist_min', 20)
            self.AsignWave2Spot(dist_min=self.dist_min )

            centers = np.array(['%0.1f - %0.1f'%(x+1,y+1) for x,y in zip(self.table['X_IMAGE'],self.table['Y_IMAGE'])])
            #  create_DS9regions2(self.table['X_IMAGE'],self.table['Y_IMAGE'], radius=2, form = 'circle',save=True,savename=self.all.filename()[:-5]+'detected',text = centers)
            
            assigned = self.table['id_slit'] != -1
            create_DS9regions([self.table['X_IMAGE'][assigned], self.table['X_IMAGE'][~assigned]],
                              [self.table['Y_IMAGE'][assigned], self.table['Y_IMAGE'][~assigned]], 
                              radius=2, form=['circle']*2, save=True, savename=self.all.filename()[:-5]+'detected', 
                              ID=[centers[assigned], centers[~assigned]], color=['green', 'white'])#vpicouet
#list(np.array((a.table['X_IMAGE'][:]).astype('string')))           
#np.array((a.table['X_IMAGE'][:]).astype('string'))+np.array((a.table['X_IMAGE'][:]).astype('string'))
#        create_DS9regions([self.table['X_IMAGE']],[self.table['Y_IMAGE']],radius=5,save=True,savename=self.all.filename()[:-5],ID=[centers])
        centers = np.array(['%0.1f - %0.1f'%(x+1,y+1) for x,y in zip(self.table['X_IMAGE'],self.table['Y_IMAGE'])])
        # create_ds9_regions(self.table['X_IMAGE'],self.table['Y_IMAGE'], radius=2, form = 'circle',save=False,savename=self.all.filename()[:-5]+'detected',text = centers)

        if plot:
            if not self.abort:
                #self.plotFWHM2D(allw=True)
                #self.plotFWHMEllipse(allw=True)
                #new_table = self.table[(self.table['FWHM_x']>1) &(self.table['FWHM_x']<25) & (self.table['FWHM_y']>1)&(self.table['FWHM_x']<25)]            
                cut = kwargs.get('cut', [50, 99.9])
                if Type =='guider':
                    self.plot_FWHM_guider(name=self.all.filename()[:-5], figsize=figsize, cut=cut)
                if Type =='detector':
                    self.plot_FWHM_det(name=self.all.filename()[:-5], figsize=figsize, cut=cut)
        if fits:
#                if type(self.line)==int:stack
#                    self.table = self.table[(self.table['ycentroid']>self.line-50) & (self.table['ycentroid']>self.line+50)]
            self.table.write(self.all.filename()[:-5] + '_table' + '.csv', format='csv',overwrite = True)
            print('Table saved on: ', self.all.filename()[:-5] + '_table' + '.csv')

        #loadDS9Image(filename = self.all.filename())
#------------------------------------------------------------------------------#
# Methods                                                                      #
#------------------------------------------------------------------------------#                     
    def pa_rotation(self, xy, pa, reverse=False):
  
        x, y = xy
        theta = 3.055*(pa - self.pa0)/4.*np.pi/180
        if reverse:
            theta = -theta
        center = np.array([.3, -101.35])[:,np.newaxis]
        gamma = 1.001
        ct = np.cos(theta)
        st = np.sin(theta)
        R = np.array([[ct, st],[-st, ct]])
        xy_rot = R.dot(np.array([x , y]) - center)*np.array([gamma,1.])[:,np.newaxis] + center
        return xy_rot


    def rotate_pa_mask(self, pa):#, R=125):
        '''
        knowing a pa rotation of the carrsousel, recompute the effetive mask coordinates
        '''
        
        # check already done
        try:
            self._mask_rotated
        except AttributeError:
            self._mask_rotated = False
            
        if not self._mask_rotated:
            print("Rotating the mask of {} degree".format(pa - self.pa0))

            # simple version
            #self.xmask -= np.sin((pa - self.pa)*np.pi/180.) * R 
#            self.xmask +=  (pa - self.pa0)*1.345
#            self._mask_rotated = True
            
            # need to compute actual rotation on sphere 
            # compute rotation
            xymask = self.pa_rotation((self.xmask, self.ymask), pa)
            self.xmask = xymask[0]
            self.ymask = xymask[1]
            self._mask_rotated = True
            
        return self.xmask, self.ymask
            
    def rotate_pa_detector(self, pa, xdet, ydet):#, R=125):
        '''
        knowing a pa rotation of the carrsousel, recompute the effetive mask coordinates
        '''
        def rotate(posx,posy,theta, xc_2 = 10542.805293946425 ,yc_2 = 901.013253389936):
            xnew = xc_2 + (posx - xc_2) * np.cos(theta*np.pi/180) + (posy - yc_2) * np.sin(theta*np.pi/180) 
            ynew = yc_2 - (posx - xc_2) * np.sin(theta*np.pi/180) + (posy - yc_2) * np.cos(theta*np.pi/180) 
            return xnew,ynew
        xnew,ynew = rotate(xdet,ydet,theta = -(pa-self.pa0)*(4/5)/1.6)
        #plt.figure()
        #plt.plot(xdet,ydet,'o')
        #plt.plot(xnew,ynew,'x')
        #plt.axis('equal')
        return xnew, ynew

    def rotate_pa_detector2(self, pa, xdet, ydet):#, R=125):
        '''
        knowing a pa rotation of the carrsousel, recompute the effetive mask coordinates
        '''
        def rotate(posx,posy,theta, xc_2 = 10542.805293946425 ,yc_2 = 901.013253389936):
            xnew = xc_2 + (posx - xc_2) * np.cos(theta*np.pi/180) + (posy - yc_2) * np.sin(theta*np.pi/180) 
            ynew = yc_2 - (posx - xc_2) * np.sin(theta*np.pi/180) + (posy - yc_2) * np.cos(theta*np.pi/180) 
            return xnew,ynew
        xnew,ynew = rotate(xdet,ydet,theta = 0.75*(pa-self.pa0),xc_2=1500,yc_2=1030)
        plt.figure()
        plt.plot(xdet,ydet,'o')
        plt.plot(xnew,ynew,'x')
        plt.axis('equal')
        return xnew, ynew - 100 * (pa-self.pa0)
    
    def plot_mapping(self):
        
        for i, w  in enumerate(self.w):
            if not np.all(self.mapping[i] == 0.): # check mapping was done 
                gx, gy = np.meshgrid(np.linspace(-12, 12, 21), np.linspace(-6, 6, 21))
                gamma = self.mapping[i][:2,:2,:]
                gamma[1,1,:] = 0.
                gxy_det0 = polyval2d(gx, gy, gamma)
                gxy_det = polyval2d(gx, gy, self.mapping[i])
                fig = plt.figure()
                ax = fig.add_subplot(111)
                for ii in range(gxy_det0.shape[1]):
                    plt.plot(gxy_det0[0,ii,:], gxy_det0[1,ii,:],'k:')
                for jj in range(gxy_det0.shape[2]):
                    plt.plot(gxy_det0[0,:,jj], gxy_det0[1,:,jj], 'k:')
                    
                dx = gxy_det[0,:,:] - gxy_det0[0,:,:]
                dy = gxy_det[1,:,:] - gxy_det0[1,:,:]
        
                xc_ = gxy_det0[0,:,:] + 10.*dx
                yc_ = gxy_det0[1,:,:] + 10.*dy
                for ii in range(gxy_det.shape[1]):
                    plt.plot(xc_[ii,:], yc_[ii,:],'b')
                for jj in range(gxy_det.shape[2]):
                    plt.plot(xc_[:,jj], yc_[:,jj], 'b')
                    
                ax.text(0.5, + 0.02, 'difference exagerated 10 times', color='b',
                        horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
                ax.set_title('science mask - detector mapping at waevelength={:.4f}'.format(w))
                ax.set_xlabel('x detector (pix)')
                ax.set_ylabel('y detector (pix)')
    
    
    def map_mask_detector(self, deg=[5, 2], bywave=True, plot=False):
        '''
        From a table of observed mask slits (or holes) on the detector, 
        determine the mapping from science mask coord in mm to detector coord in mm
        and compute as well, the mask center position on the detector for each wavelength

        Focus.map_mask_detector(deg=[5, 2], bywave=True, plot=False)
        
            deg: degree of polynome [deg_x, deg_y]
            
            bywave: if True (default) perform a mapping by wavelength
                    else include w in the fit (deg 2)
                    
            plot: if True plot the residual vectors
            
        return the array of centers  the array of polynomial coeffs.
        '''    
        
        self.centers = np.full((self.w.size, 2), -1.0)

        if bywave:
            
            self._mapping_by_wave = True
            self.mapping = Mapping(wavelength=self.w)
            subt = self.table[self.table['wavelength'] > 0]
            sid = subt['id_slit']
            self.mapping.set(subt['wavelength'], self.xmask[sid], self.ymask[sid], subt['X_IMAGE'], subt['Y_IMAGE'] , deg)

            # compute direct mapping residuals
            pred = self.mapping.map(subt['wavelength'], self.xmask[sid], self.ymask[sid])
            dx = subt['X_IMAGE'] - pred[0]
            dy = subt['Y_IMAGE'] - pred[1]
            for w in self.w:
                flagw = (subt['wavelength'] == w)
                print("at wavelength: {}".format(w))
                print("maximum residual of direct mapping along x, y (in pixels): {}, {}".format(
                        np.abs(dx[flagw]).max(), np.abs(dy[flagw]).max()))
                
            # compute inverse mapping residuals
            pred = self.mapping.inv_map(subt['wavelength'], subt['X_IMAGE'], subt['Y_IMAGE'])
            dx = self.xmask[sid] - pred[0]
            dy = self.ymask[sid] - pred[1]        
            for w in self.w:
                flagw = (subt['wavelength'] == w)
                print("at wavelength: {}".format(w))
                print("maximum residual of inverse mapping along x, y (in mm): {}, {}".format(
                        np.abs(dx[flagw]).max(), np.abs(dy[flagw]).max()))

            # get center
            self.centers = self.mapping.map(self.w, 0., 0.)

            if plot:
                plt.figure()
                plt.quiver(self.xmask[sid], self.ymask[sid], dx, dy)    
                             
            
        return self.mapping, self.centers
    
    
    def apply_map_mask_detector(self, xmask, ymask, wavelength):
        
        xmask = np.array(xmask)
        ymask = np.array(ymask)
        
        # check if mapping by wave or not 
        xy_det = self.mapping.map(wavelength, xmask, ymask)

        return xy_det

        
    def apply_map_detector_mask(self, xdet, ydet, wavelength):
        
        xdet = np.array(xdet)
        ydet = np.array(ydet)
        
        # check if mapping by wave or not 
        xy_mask = self.mapping.inv_map(wavelength, xdet, ydet)

        return xy_mask
        

    def AsignWave2Spot(self, dist_min=20, dist_max=30):
        
        self.table['wavelength'] = 0.0
        self.table['id_slit'] = -1
        self.table['redshift'] = -1.0

        pspot = self.predictedSpots.copy()
        for i in range(len(pspot)):
            distance2spot = distance(pspot['x'][i],pspot['y'][i],pspot['x'],pspot['y'])
            if len(pspot[distance2spot<dist_max])>1:
                pspot[i]['w'] = -1
                pspot[i]['id_slit'] = -1
                pspot[i]['redshift'] = -1
            
        for i in range(len(self.table)):
            distance2spot = distance(self.table['X_IMAGE'][i],self.table['Y_IMAGE'][i],pspot['x'],pspot['y'])
            if len(pspot[distance2spot<dist_min])>0:

                #print('the minimal distance is : ', distance_min)
                id = np.argmin(distance2spot)
                self.table[i]['wavelength'] = pspot[id]['w']
                #print(self.table[i]['wavelength'], pspot[id]['w'])
                self.table[i]['id_slit'] = pspot[id]['id_slit']
                #print (pspot[id]['Internal-count'])
                self.table[i]['Internal-count'] = pspot[id]['Internal-count']
                self.table[i]['xmask'] = pspot[id]['xmask']
                self.table[i]['ymask'] = pspot[id]['ymask']
                self.table[i]['redshift'] = pspot[id]['redshift']
                #self.table[i]['ymask'] = pspot[id]['ymask']
                #print('yes', i, id, pspot[id]['w'], self.table[i]['wavelength'] )

                #print(id)
        print('Number of wavelength asignated spots: ', len(self.table[self.table['wavelength']>0]))
        return 


    def computeLinesLocation(self):
        """
        predict detector position of lines using IMO mapping
        """
        
        mapping_sky2mask = Mapping(self.sky2mask_mapping)
        mapping_sky2det = Mapping(self.sky2det_mapping)
        
        self.internalCount = None
        #### check with slit positions                
        slitfile = self.slit_location#'/Users/Vincent/Nextcloud/FIREBALL/Tests-at-FortSumner/170923_SC_GUI03/Slits/targets_F1.txt'
        if 'targets_F1' in self.slit_location:
            print('Computing predicted location for F1')
            print(slitfile)
            slits = Table.read(slitfile, format='ascii', delimiter='\t')
            # remove multi object slits
            idok = (slits['slit_length_right'] != 0) &  (slits['slit_length_left'] !=0)
            slits = slits[idok]
            self.xmask = slits['xmask'] + (slits['slit_length_right'] - slits['slit_length_left'])/2.
            self.ymask = slits['ymask'] + slits['offset']
            self.z = slits['z'] 
            self.internalCount = slits['Internal-count']
    
        elif 'grid_mask' in self.slit_location:
            print('Computing predicted location for grid mask')
            slits =  Table.read(slitfile, format='ascii', delimiter='\t')
            self.xmask = slits['x']
            self.ymask = slits['y']
            self.z = 0
            #x = xmask
            #y = -ymask
    
        elif ('F2' in self.slit_location) or ('F3' in self.slit_location) or ('F4' in self.slit_location) or ('Tilted' in self.slit_location):
            print('Computing predicted location for science mask')
            slits =  Table.read(slitfile, format='ascii', delimiter='\t')
            self.xmask = slits['xmm']
            self.ymask = slits['ymm']
            self.internalCount = slits['Internal-count']
            #self.mainID = slits['#main_id']
            if not ('Tilted' in self.slit_location):
                self.z = slits['Z'] 

        # eventually rotate mask along pa
        if self.pa != self.pa0:
            self.rotate_pa_mask(self.pa)
            
        # convert mask axis convention to zemax/IMO convention
        xzmx = self.xmask
        yzmx = - self.ymask  
    
        if self.source == 'Zn':
            self.w = self.zinc_lines
            w = self.w[:,np.newaxis]            

        elif self.source == 'Lya':
            print('need to have the redshift with the sources')
            self.w = np.array([self.lya_line])
            w = np.array([(1 + self.z)*self.lya_line])[np.newaxis,:]
            
        elif self.source == 'all':
            self.w = self.zinc_lines
            self.w = np.insert(self.w, 0, self.lya_line)
            w = np.tile(self.w[:,np.newaxis], (1,len(xzmx)))
            w[0,:] = (1 + self.z)*1215.67*1e-4
        print('Wavelength = ',w)    
        xzmx = xzmx[np.newaxis,:]
        yzmx = yzmx[np.newaxis,:]
           
            
        # compute detector position
        xysky = mapping_sky2mask.inv_map(w, xzmx, yzmx)
        xsky = xysky[0]
        ysky = xysky[1]
        xydet = mapping_sky2det.map(w, xsky, ysky)
        xdet_zmx = xydet[0]
        ydet_zmx = xydet[1]
        
        # convert zemax/IMO to fits detector
        xdet, ydet = -ydet_zmx, xdet_zmx
                    
        print('Location of the sources predicted for lambda=', self.w)
        return xdet, ydet, w
    

    def putOffsetMaskDetect(self, xdetmm, ydetmm, offset, pix_size=0.013):
#            new_xdetpix = -ydetmm  * mag - offset + self._DFLT_Yinf - (1522-1217)
#            new_ydetpix = ( xdetmm * mag) - np.nanmin(-xdetmm * mag) - (991-999)
#        else:
        new_xdetpix = xdetmm / pix_size + self._DFLT_Yinf - offset[0]#- (1522-1217)- a
        new_ydetpix = ydetmm / pix_size - offset[1]#- (991-999)- np.nanmin(-xdetmm * mag)

        return new_xdetpix, new_ydetpix 

            
    def extractSourcesPy(self,quick, Type, threshold, fwhm, centre, shape, **kwargs):
        """Extract sources from a given image using DAO star finder
        If there is some defocus in the image, user should use uick=True so that we run the tracker with several instrument PSF and noise threshold
        """ 
        from photutils import DAOStarFinder

        #verbose = kwargs.get('verbose', False)
        n=2
        data = self.image
        data2 = ndimage.grey_dilation(ndimage.grey_erosion(data, size=(n,n)), size=(n,n))       
        self.mean, self.median, self.std = sigma_clipped_stats(data2, sigma=3.0, iters=5)    
        
        if quick:
            #if Type == 'guider':
            if (shape == 'holes') or (shape == 'gaussian'):
                print('Detection algorith is used without any ellipticity')
                daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*self.std, ratio = 1)   #15,10      
            if (shape == 'slits') or (shape == 'slits2D'):
                print('Detection algorith is used without 66% ellipticity')
                daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*self.std,ratio = 0.5, theta = 90)         
            #if Type == 'detector':
            #    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*self.std, ratio = 1) 
            sources0 = daofind(data2 - self.median) 
        if centre:
            sources0.remove_rows(slice(0,-1))
            sources0[0]['xcentroid'] = centre[0]
            sources0[0]['ycentroid'] = centre[1]
            
        if not quick:
            daofind = DAOStarFinder(fwhm=fwhm[0], threshold=threshold[0]*self.std,ratio = 0.5, theta = 90)        
            sources0 = daofind(data2 - self.median)
            for i in fwhm: 
                for j in threshold: 
                    if (shape == 'holes') or (shape == 'gaussian'):
                        print('Detection algorith is used without any ellipticity')
                        daofind = DAOStarFinder(fwhm=i, threshold=j*self.std, ratio = 1)   #15,10      
                    if (shape == 'slits') or (shape == 'slits2D'):
                        print('Detection algorith is used without 66% ellipticity')
                        daofind = DAOStarFinder(fwhm=i, threshold=j*self.std,ratio = 0.5, theta = 90)         
                    
                    sources1 = daofind(data2 - self.median) 
                    print('fwhm = {}, T = {}, len = {}'.format(i,j,len(sources1)))
                    try :
                        sources0 = Table(np.hstack((sources0,sources1)))
                    except TypeError:
                        print('catalog empty')
                        if len(sources0)==0:
                            sources0 = sources1
                        
        self.sources = delete_doublons(sources0, dist = 30)#sources#delete_doublons(sources0, dist = 10)
        
        for value in self.all[0].header[:-2]:
            if (value != 'COMMENT'):# & (value != 'TEMP'):
                try:
                    print(value)
                    self.sources[value] = self.all[0].header[value]
                except TypeError:
                    pass
            
        try:
            self.sources['X_IMAGE'] = self._DFLT_Yinf + self.sources['xcentroid']
            self.sources['Y_IMAGE'] = self._DFLT_Xinf + self.sources['ycentroid']
            self.sources['NUMBER'] = self.sources['id']
            self.sources['xmask'] = float(0)
            self.sources['ymask'] = float(0)
            try:
                self.sources['Image Number'] = float(os.path.basename(self.all.filename())[-11:-4])
            except ValueError:
                pass
                
            self.sources['DATE'] = self.all[0].header['DATE']
                
            self.sources['Internal-count'] = '-'*15
            for xy in ['x', 'y']:
                self.sources['FWHM_' + xy] = float(0)
                self.sources['Amplitude_after_stack_' + xy] = float(0)
                self.sources['Offset_after_stack_' + xy] = float(0)
                self.sources['SumPix_' + xy] = float(0)
                self.sources['MaxPix_' + xy] = float(0)
                self.sources['SlitSize_' + xy] = float(0)
                
                self.sources['Background_' + xy] = float(0)
                self.sources['Background_var_' + xy] = float(0)
                
                self.sources[xy.upper() + '_IMAGE_var'] = float(0)
                self.sources['SlitSize_var_' + xy] = float(0)
                self.sources['FWHM_var_' + xy] = float(0)
                self.sources['Amplitude_after_stack_var_' + xy] = float(0)
                self.sources['Offset_after_stack_var_' + xy] = float(0)
    
                self.sources['FlagFit_' + xy] = int(0)
                
            
            # filter sources with peak
            peak_threshold = kwargs.get('peak_threshold', 'Automatic')
            if peak_threshold is not None:
                if type(peak_threshold) == int:
                #self.sources = self.sources[self.sources['peak'] > peak_threshold]
                    print("Initial number of sources: {}".format(len(self.sources)))                
                    print("Filtering out sources with peak flux below {}".format(peak_threshold))
                    selection = self.sources['peak'] > peak_threshold
                    self.table = self.sources[selection].copy()    #[self.sources['peak'] > peak_threshold]                    
                    print("Number of sources left: {}".format(len(self.table)))
                else:
                    MoreSources = kwargs.get('MoreSources',30)
                    #self.sources = self.sources[self.sources['peak'] > peak_threshold]
                    print("Initial number of sources: {}".format(len(self.sources)))
                    print("Automatic Filtering out sources with peak flux")
                    n = len(self.predictedSpots)
                    print("Using number of predicted sources ({}) + {}: {}".format(n,MoreSources,n+MoreSources))
                    self.sources.sort('peak')
                    self.table = self.sources[-n-MoreSources:]  
            else:
                self.table = self.sources
                
        except TypeError:
            print('no source detected')
            self.abort = True
        else:
            self.abort = False
        return

#        figure()
#        plt.imshow(data2)#, origin='lower', norm=norm)stack#, cmap='Greys', origin='lower', norm=norm)
#        apertures.plot(color='blue', lw=5.5, alpha=1.5)#sources1['X_IMAGE'] = sources1['xcentroid'] 
#        plt.title('Image - median after erosion dilatation' )
#        plt.show()

#        mask = make_source_mask(self.image, snr=1.5, npixels=5, dilate_size=50)#
##        figure();imshow(self.image - median);colorbar();plt.title('Image - median without sources' );show()
#        mean, median, std = sigma_clipped_stats(self.image, sigma=3.0, mask=mask)
#        print('[STATISTICS :] mean = {}, median = {}, std = {} before masking sources'.format(mean, median, std) )
#        self.image = self.image - median
##        figure();imshow(self.image - median);colorbar();show()
        
        
    def extractSources(self, **kwargs):
        """Extract source using Sextractor, 
        should not be used
        """ 
        fname = kwargs.get('filename', False)
        os.system('sex '+ fname + 'm' + self._DFLT_extension + ' -CATALOG_NAME cat'+ fname + '.fits' + ' -CHECKIMAGE_NAME  ' + fname + '.fits.fits')
        
    def extractCatalog(self, **kwargs):
        """Extract sources from a catalog of detection generated with sextractor
        """
        fname = kwargs.get('filename', False)
        if fname:
            self.table = Table.read(self._DFLT_catalog + fname + self._DFLT_extension)
            self.table = self.table[(self.table['FLUX_ISO']<1e5) & (self.table['FWHM_IMAGE']>2) & (self.table['FWHM_IMAGE']<25) & (self.table['ISOAREA_IMAGE']>30) & (self.table['ISOAREA_IMAGE']<300)]

    def transform_catalog(self, **kwargs):
        """Extract sources from a catalog of detection generated with sextractor
        """
        fname = kwargs.get('filename', False)
        if fname:
            self.table = Table.read(fname.replace(".fits","_cat.fits"))
            
            self.sources['X_IMAGE'] = self._DFLT_Yinf + self.sources['xcentroid']
            self.sources['Y_IMAGE'] = self._DFLT_Xinf + self.sources['ycentroid']
            self.sources['NUMBER'] = self.sources['id']
            self.sources['xmask'] = float(0)
            self.sources['ymask'] = float(0)
            try:
                self.sources['Image Number'] = float(os.path.basename(self.all.filename())[-11:-4])
            except ValueError:
                pass
                
            self.sources['DATE'] = self.all[0].header['DATE']
                
            self.sources['Internal-count'] = '-'*15
            for xy in ['x', 'y']:
                self.sources['FWHM_' + xy] = float(0)
                self.sources['Amplitude_after_stack_' + xy] = float(0)
                self.sources['Offset_after_stack_' + xy] = float(0)
                self.sources['SumPix_' + xy] = float(0)
                self.sources['MaxPix_' + xy] = float(0)
                self.sources['SlitSize_' + xy] = float(0)
                
                self.sources['Background_' + xy] = float(0)
                self.sources['Background_var_' + xy] = float(0)
                
                self.sources[xy.upper() + '_IMAGE_var'] = float(0)
                self.sources['SlitSize_var_' + xy] = float(0)
                self.sources['FWHM_var_' + xy] = float(0)
                self.sources['Amplitude_after_stack_var_' + xy] = float(0)
                self.sources['Offset_after_stack_var_' + xy] = float(0)
    
                self.sources['FlagFit_' + xy] = int(0)
            # self.table = self.table[(self.table['FLUX_ISO']<1e5) & (self.table['FWHM_IMAGE']>2) & (self.table['FWHM_IMAGE']<25) & (self.table['ISOAREA_IMAGE']>30) & (self.table['ISOAREA_IMAGE']<300)]

#	    print(len(self.table), 'spots found at positions' , self.table['X_IMAGE'], self.table['Y_IMAGE'])


    def Fit2D(self, HumanSupervision, n=20, **kwargs):
        #quick = kwargs.get('quick', True)
        #verbose = kwargs.get('verbose', False)
        popt = []
        pcov = [] 
        for i in range(len(self.table)):#len(self.table)):
#            print(i)
            image = self.image[int(self.table[i]['ycentroid']-n):int(self.table[i]['ycentroid']+n),int(self.table[i]['xcentroid']-n):int(self.table[i]['xcentroid']+n)]            
            image = ndimage.grey_dilation(ndimage.grey_erosion(image, size=(2,2)), size=(2,2))
            lx,ly = image.shape
            x = np.linspace(0,lx-1,lx)
            y = np.linspace(0,ly-1,ly)
            x, y = np.meshgrid(x,y)
            try:
                a = curve_fit(twoD_Gaussian,(x,y),image.flat,(np.nanmean(image),lx/2,ly/2,1,1,np.nanmin(image)))
                popt.append(a[0])
#            popt.append(Normalization(*a[0]))
                pcov.append(a[1])
            except RuntimeError:
                print('Optimal parameters not found... Trying next image')
            if HumanSupervision:
                fig, axes = plt.subplots(1,2)
                x = np.linspace(0,2*n-1,2*n)
                y = np.linspace(0,2*n-1,2*n)
                x, y = np.meshgrid(x,y)
                im = axes[0].imshow(image, interpolation = 'none')
                cs1 = axes[0].contour(image,popt[i][-1]+abs(popt[i][0]-popt[i][-1])/2,linewidths=3,colors=('black'),levels=[int(5)])
                axes[0].clabel(cs1, inline=0.1, fontsize=12)
                im = axes[1].imshow(twoD_Gaussian((x,y), *popt[i]).reshape(2*n,2*n), interpolation = 'none')
                cs2 = axes[1].contour(image,abs(popt[i][0]+popt[i][-1])/2,linewidths=3,colors=('black'),levels=[int(5)])
                axes[1].clabel(cs2, inline=0.1, fontsize=12)
                cax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
                fig.colorbar(im, cax=cax,ax=axes.ravel().tolist(),orientation='vertical')
                fig.suptitle('FWHMx = {} \n FWHMy = {}'.format(2.355*popt[i][3],2.355*popt[i][4]))
                plt.show()
        return popt,pcov


    def Fit2D_Convolved(self, HumanSupervision, n=20, **kwargs):
        #quick = kwargs.get('quick', True)
        #verbose = kwargs.get('verbose', False)
        from DS9Utils import ConvolveSlit2D_PSF
        popt = []
        pcov = [] 
        for i in range(len(self.table)):#len(self.table)):
#            print(i)
            image = self.image[int(self.table[i]['ycentroid']-n):int(self.table[i]['ycentroid']+n),int(self.table[i]['xcentroid']-n):int(self.table[i]['xcentroid']+n)]            
            image = ndimage.grey_dilation(ndimage.grey_erosion(image, size=(2,2)), size=(2,2))
            image = (image - image.min()) / (image - image.min()).max()
            lx,ly = image.shape
            x = np.linspace(0,lx-1,lx)
            y = np.linspace(0,ly-1,ly)
            x, y = np.meshgrid(x,y)
            try:
                ConvolveSlit2D_PSF_75muWidth = lambda xy , amp, L, xo, yo, sigmax2, sigmay2: ConvolveSlit2D_PSF(xy, amp, 2.5, L, xo, yo, sigmax2, sigmay2)
                popt1, pcov1 = curve_fit(ConvolveSlit2D_PSF_75muWidth,(x,y),image.flat,(1,9,n,n,16,16))
                popt1 = list(popt1)
                popt1.insert(1,3)
                popt.append(popt1)
#            popt.append(Normalization(*a[0]))
                pcov.append(pcov1)
            except RuntimeError:
                print('Optimal parameters not found... Trying next image')
            if HumanSupervision:
                fig, axes = plt.subplots(1,2)
                x = np.linspace(0,2*n-1,2*n)
                y = np.linspace(0,2*n-1,2*n)
                x, y = np.meshgrid(x,y)
                im = axes[0].imshow(image, interpolation = 'none')
                #cs1 = axes[0].contour(image,popt[i][-1]+abs(popt[i][0]-popt[i][-1])/2,linewidths=3,colors=('black'),levels=[int(5)])
                #axes[0].clabel(cs1, inline=0.1, fontsize=12)
                im = axes[1].imshow(ConvolveSlit2D_PSF((x,y), *popt[i]).reshape(2*n,2*n), interpolation = 'none')
                #cs2 = axes[1].contour(image,abs(popt[i][0]+popt[i][-1])/2,linewidths=3,colors=('black'),levels=[int(5)])
                #axes[1].clabel(cs2, inline=0.1, fontsize=12)
                cax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
                fig.colorbar(im, cax=cax,ax=axes.ravel().tolist(),orientation='vertical')
                fig.suptitle('Lx = {} \n Ly = {}\n FWHM_x = {}, FWHMy = {} '.format(2*popt1[1],2*popt1[2],2.355*np.sqrt(popt1[5]),2.355*np.sqrt(popt1[6])))
                plt.show()
                
            self.table['X_IMAGE'][i] += popt1[3]
            self.table['Y_IMAGE'][i] += popt1[4]
            self.table['FWHM_x'][i] = 2.355*np.sqrt(popt1[5])
            self.table['FWHM_y'][i] = 2.355*np.sqrt(popt1[5])
            self.table['SlitSize_X'][i] = 2*popt1[1]
            self.table['SlitSize_Y'][i] = 2*popt1[2]
            
        return popt,pcov


    def Fit1D(self, HumanSupervision, stack, verbose,Type, shape= 'slits', nstack=3, min_fwhm=2., convolved=True, axis=0, nmorph=0):
        """Fit a detected spot with a gaussian or unit box convolved by gaussian  with or without human supervision (plot=True)
        press n if you don't want to keep the detection, enter or any other letter if you want to keep it
        axis = 0 means X axis (spectral), 1 means Y axis (spatial)
        """
        if Type == 'guider':
            npad = 9
        if Type == 'detector':
            npad = 30

        popt = []
        pcov = [] 
        e=0
        if axis==0:
            imagepad = self.image
            cx = self.table['xcentroid'].round().astype(int)
            cy = self.table['ycentroid'].round().astype(int)
            xykey = 'x'
            label = 'Slit spectral profile'
        else:
            imagepad = self.image.T
            cx = self.table['ycentroid'].round().astype(int)
            cy = self.table['xcentroid'].round().astype(int)
            xykey = 'y'
            label = 'Slit spatial profile'
        imagepad = np.pad(imagepad,  ((0,0), (npad,npad)), 'reflect', reflect_type='odd')
        for i in range(len(self.table)):
            image = imagepad[cy[i]-nstack:cy[i]+nstack+1, cx[i]:cx[i]+2*npad+1]
            if nmorph>0:
                image = ndimage.grey_dilation(ndimage.grey_erosion(image, size=(nmorph,nmorph)), size=(nmorph,nmorph))
            image = np.nansum(image, axis=0)
            lx = image.shape[0]
            x = np.linspace(-(lx-1)/2,(lx-1)/2,lx)#x = np.linspace(0,lx-1,lx)
            image = (image-image.min())/(image-image.min()).max()
            min_sigma = min_fwhm/2.35
            if (shape == 'slits') or (shape == 'slits2D') or (shape == 'gaussian'):
                decenter = 7
                if axis==0:
                    size = [2.8, 5]
                    size0 = 3.
                else:
                    size = [8., 12.]
                    size0 = 10.
            if shape == 'holes':
                decenter = 2
                size = [1, 5]
                size0 = 3.
            p0 = [1., 0., 3., np.median(image)]
            bounds = ([0, -decenter , min_sigma, -np.inf], [np.inf, decenter, np.inf, np.inf])
            if convolved:
                model = ConvolveBoxPSF
                p0.insert(1, size0)
                bounds[0].insert(1, size[0])
                bounds[1].insert(1, size[1])
            else:
                model = Gaussian 
            try:
                poptww, pcovww = curve_fit(model, x ,image.flat, p0=p0,  bounds=bounds)
                # ConvolveBoxPSF(x, amp=1, l=40, x0=0, sigma2=40, offset=0)
            except RuntimeError or ValueError:
                e += 1
                popt.append(None)
                pcov.append(None)
                print('Image non analyzed: ', e                )
                self.table['FWHM_' + xykey][i] = 0  
                print(RuntimeError, ValueError)
            except ValueError:
                e += 1
                popt.append(None)
                pcov.append(None)
#                plt.figure()
#                plt.plot(x, image.flat)
#                plt.show()
                print('Value error... Trying next image',i)          
                self.table['FWHM_' + xykey][i] = 0               
            else:
                if convolved:
                    amp, l, x0, sigma2, offset = poptww
                    amp_var, l_var, x0_var, sigma2_var, offset_var = tuple(np.diag(pcovww))
                else:
                    amp, x0, sigma2, offset = poptww
                    amp_var, x0_var, sigma2_var, offset_var = tuple(np.diag(pcovww))
                    l = 0
                backgroundx = (x < x0-l-3*np.sqrt(sigma2)) | (x > x0+l+3*np.sqrt(sigma2))
                background = image[backgroundx]
                popt.append(poptww)
                pcov.append(np.trace(pcovww))
                if HumanSupervision:
                    plt.figure()
                    plt.plot(x, image, label=label)#xykey
                    #plt.plot(x[backgroundx], background, '+')
                    plt.plot(x, model(x, *poptww),label=r'Fit: PSF*Slit') # convolved profile
                    if convolved:
                        max_profile = ConvolveBoxPSF(x0, amp, l, x0, sigma2, 0.)
                        max_gaussian = Gaussian(0., 1., 0., sigma2, 0.)
                        amp_gaussian = max_profile / max_gaussian 
                        plt.plot(x, Gaussian(x, amp_gaussian, x0, sigma2, offset ), ':b',label='Deconvolved PSF') # Gaussian x, amplitude, xo, sigma_x, offset
                        xc = x - x0
                        plt.plot(x, np.piecewise(x, [xc < -l, (xc >=-l) & (xc<=l), xc>l], [offset, max_profile + offset, offset]), ':r', label='Slit size') # slit (UnitBox)
                    plt.plot([x0, x0], [0, 1])
                    #plt.title('%0.0f,   sigma = %0.3f +/- %0.3f pix\nSlitdim = %0.3f +/- %0.3f pix\ncenter = %0.3f +/- %0.3f \nPeak = %0.3f' % (i, np.sqrt(sigma2), np.sqrt(sigma2_var/2.) , 2*l,2*np.sqrt(l_var), x0, np.sqrt(x0_var),self.table[i]['peak'] ))
                    plt.title('Image %0.0f' % (i))
                    plt.ylabel(label)
                    plt.figtext(0.57,0.45,'Sigma = %0.1f +/- %0.1f pix\nSlitdim = %0.1f +/- %0.1f pix\ncenter = %0.1f +/- %0.1f \nPeak = %0.1f' % ( np.sqrt(sigma2), np.sqrt(sigma2_var/2.) , 2*l,2*np.sqrt(l_var), x0, np.sqrt(x0_var),self.table[i]['peak'] ),bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})
                    #plt.figtext(0.60,0.50,'FWHM = %0.1f +/- %0.1f pix\nSlit = %0.1f +/- %0.1f pix' % ( 2.35*np.sqrt(sigma2), 2.35*np.sqrt(sigma2_var/2.) , 2*l,2*np.sqrt(l_var) ),bbox={'facecolor':'blue', 'alpha':0.2, 'pad':10})
                    plt.legend(loc = 'upper right')
                    #vincent
                    plt.show()
                    if sys.version_info.major == 3:
                        var = input("Do you want to keep this fit? : ")
                    if sys.version_info.major == 2:
                        var = raw_input("Do you want to keep this fit? : ")
                    if var != 'n':
                        print("We kept this fit!")
                        self.table['FlagFit_' + xykey][i] = int(1)
                    if var == 'n':
                        print("Fit not taken into account, flag=0")
               
                self.table['Background_' + xykey][i] = background.mean()
                self.table['Background_var_' + xykey][i] = background.var()
                self.table[xykey.upper() + '_IMAGE'][i] += x0
                self.table[xykey.upper() + '_IMAGE_var'][i] = x0_var
                self.table['FWHM_' + xykey][i] = 2.355*np.sqrt(sigma2)      
                self.table['FWHM_var_' + xykey][i] = sigma2_var/2 * 2.35**2
                self.table['Amplitude_after_stack_' + xykey][i] = amp
                self.table['Amplitude_after_stack_var_' + xykey][i] = amp_var
                self.table['Offset_after_stack_' + xykey][i] = offset 
                self.table['Offset_after_stack_var_' + xykey][i] = offset_var 
                if convolved:
                    self.table['SlitSize_' + xykey][i] = 2*l
                    self.table['SlitSize_var_' + xykey][i] = 4*l_var
                n3=50
                self.table['SumPix_' + xykey][i] = np.nansum(self.image[int(cy[i]-n3):int(cy[i]+n3),int(cx[i]-n3):int(cx[i]+n3)]) - np.nansum(self.image[100-n3:100+n3,100-n3:100+n3])
                self.table['MaxPix_' + xykey][i] = np.nanmax(image)
            
        return popt,pcov



    def Fit1D2_CC(self, HumanSupervision, stack, verbose, n1=9, n2=3,  **kwargs):#    def Fit1D2(self, plot, stack, verbose, n1=9, n2=3,  **kwargs):#12
#        popt = []
#        pcov = [] 
        n=0
#        e=0  
        if verbose:
            print(len(self.table['spectral resolution']), 'images to analyze')
        for i in range(len(self.table)):
            if verbose:
                print('Image : ', i)
            image = self.image[int(self.table[i]['ycentroid']-n1):int(self.table[i]['ycentroid']+n1),int(self.table[i]['xcentroid']-n2):int(self.table[i]['xcentroid']+n2)]
#            image = self.image
            #background = np.nansum(ndimage.grey_dilation(ndimage.grey_erosion(self.image[50-n1:50+n1,50-n2:50+n2], size=(n,n)), size=(n,n)), axis=0)
#            image = np.nansum(ndimage.grey_dilation(ndimage.grey_erosion(image, size=(n,n)), size=(n,n)), axis=1)
            image = ndimage.grey_dilation(ndimage.grey_erosion(image, size=(n,n)), size=(n,n))
            lx = image.shape[0]
            x = np.linspace(0,lx-1,lx)

            x = np.arange(-n1,n1)
            im1x = np.nansum(image,axis = 1)

            centre = np.linspace(n1/4,2*n1 - n1/4,10*n1)
            gaussians = GaussianM(x[np.newaxis,:],centre[:,np.newaxis],1)
            crosscorx1 = np.nansum(gaussians * im1x,axis=0)
            centroid1xM = centre[crosscorx1.argmax() ]
            if HumanSupervision:
                plt.figure();plt.title(file)
                plt.plot(x,GaussianM(x,centroid1xM,im1x.max()))
                #plt.plot(x,Gaussian(x,*centroid2x),label = 'cov = {}'.format(cov2xM[0,0]));plt.plot(x,Gaussian(x,*centroid2y),label = 'cov = {}'.format(cov2yM[0,0]))
                plt.legend(loc='upper right')
                plt.plot(x,im1x,'x')
                plt.show()
            self.table[i]['X_IMAGE'] = self.table[i]['xcentroid'] + centroid1xM
        return True


    def createCatalog(self, number, **kwargs):
        """Save the catalog with the detected sources
        """
        fname = kwargs.get('filename', False)
        if 'FLUX_ISO' in self.table.colnames:
            self.table.write(self._DFLT_catalog + fname + self._DFLT_extension,overwrite = True)
        if 'xcentroid' in self.table.colnames:        
            self.table.write(self._DFLT_catalog + 'py{}'.format(number) + '.csv', format = 'csv',overwrite = True)
            
    def CalculateDispersion(self):
        x=[]
        y=[]
        disp=[]
        for i in np.arange(self.table['id_slit'].max()):
            if len(self.table[self.table['id_slit']==i])==2:
        #        print(abs(self.table[self.table['id_slit']==i][0]['X_IMAGE'] - self.table[self.table['id_slit']==i][1]['X_IMAGE'] ))   
        #        print(self.table[self.table['id_slit']==i]['wavelength'])
                distance = abs(self.table[self.table['id_slit']==i][0]['X_IMAGE'] - self.table[self.table['id_slit']==i][1]['X_IMAGE'] ) 
                ldistance = float(abs(self.table[self.table['id_slit']==i][0]['wavelength'] - self.table[self.table['id_slit']==i][1]['wavelength'] )) *1000
                if (float(distance) / ldistance) < 47.1:
                    x.append(self.table[self.table['id_slit']==i]['X_IMAGE'].mean())
                    y.append(self.table[self.table['id_slit']==i]['Y_IMAGE'].mean())
                    disp.append(float(distance) / ldistance)    
            if len(self.table[self.table['id_slit']==i])==3:
                distance = abs(self.table[(self.table['id_slit']==i) & (self.table['wavelength']>0.2050)][0]['X_IMAGE'] - self.table[(self.table['id_slit']==i) & (self.table['wavelength']>0.2050)][1]['X_IMAGE'] ) 
                ldistance = float(abs(self.table[(self.table['id_slit']==i) & (self.table['wavelength']>0.2050)][0]['wavelength'] - self.table[(self.table['id_slit']==i) & (self.table['wavelength']>0.2050)][1]['wavelength'] )) * 1000
                if (float(distance) / ldistance) < 47.1:
                    x.append(self.table[(self.table['id_slit']==i) & (self.table['wavelength']>0.2050) ]['X_IMAGE'].mean())
                    y.append(self.table[(self.table['id_slit']==i) & (self.table['wavelength']>0.2050) ]['Y_IMAGE'].mean())
                    disp.append(float(distance) / ldistance)                
                if (float(distance) / ldistance) < 47.1:
                    distance = abs(self.table[(self.table['id_slit']==i) & (self.table['wavelength']<0.2070)][0]['X_IMAGE'] - self.table[(self.table['id_slit']==i) & (self.table['wavelength']<0.2070)][1]['X_IMAGE'] ) 
                    ldistance = float(abs(self.table[(self.table['id_slit']==i) & (self.table['wavelength']<0.2070)][0]['wavelength'] - self.table[(self.table['id_slit']==i) & (self.table['wavelength']<0.2070)][1]['wavelength'] )) * 1000   
                    x.append(self.table[(self.table['id_slit']==i) & (self.table['wavelength']<0.2070) ]['X_IMAGE'].mean())
                    y.append(self.table[(self.table['id_slit']==i) & (self.table['wavelength']<0.2070) ]['Y_IMAGE'].mean())
                    disp.append(float(distance) / ldistance)    
                #plt.plot(x,y,'x')
        n=12
        a = 1200
        b = 2000#900
        c = 100
        d = 1800#1600
        x1 = np.linspace(a,b,n)
        y1 = np.linspace(c,d,n)
        xx, yy = np.meshgrid(x1,y1)
        dispmap = np.zeros((n,n))
        for k in range(len(disp)):   
            for i in range(n):
                for j in range(n):
                    if (x[k]>x1[i]) &  (x[k]<x1[i] + (b-a)/(n-1)) & (y[k]>y1[j]) &  (y[k]<y1[j] + (d-c)/(n-1)):
                        dispmap[i,j] = disp[k]
        dispb = np.zeros((2*n,n))
        for i in range(n):
            dispb[2*i,:] = dispmap[i,:]
            dispb[2*i+1,:] = dispmap[i,:]
        dispb[dispb==0] = np.nan
        plt.figure()
        plt.title('Mask to det dispersion measured with zinc lamp')
        plt.imshow(13*dispb[:-2,:].T,cmap='coolwarm')#cm.PRGn)#,vmin=45.9, vmax=46.5)cmap='inferno'
        #plt.scatter(np.array(x1),np.array(y1),'o')
        #plt.yticks(np.arange(n),x.astype(intl))
        #plt.xticks(np.arange(n),y.astype(int))
        plt.xlabel('Field of view x-coordinate (arcmin)' )
        plt.ylabel('Field of view y-coordinate (arcmin)')
        plt.yticks(np.linspace(0,9,5),[-450/60,-225/60,0,225/60,450/60])
        plt.xticks(np.linspace(0,18,5),[-780/60,-390/60,0,390/60,780/60])
        cb = plt.colorbar(orientation='horizontal')
        cb.set_label("dispersion [microns/nanometer]")
        #plt.show()
        return x,y,disp

#    def plot_FWHM_guider(self, name, figsize=10, cut=[50, 99.9]):
#        """Function that plot a guider image of the Image class nicely 
#        """
#        n=500
#        table = self.table
#        if type(self.line) is int:
#            image = self.all[0].data[self.line-n:self.line+n, :] 
#            decalage = self.line - n
#        else:
#            image = self.image
#            #self.line = 0.0
#            decalage = 0
#        ly, lx = image.shape
#
#        fs = np.array([lx,ly], dtype=float) 
#        fs *= figsize/fs.max()
#        fig = plt.figure(figsize=fs)
#        axIma = plt.subplot(111)
#        norm = ImageNormalize(image, interval=AsymmetricPercentileInterval(cut[0], cut[1]), stretch=AsinhStretch())   
#        axIma.imshow(image, cmap = 'cubehelix', origin='lower', norm=norm)
#        axIma.set_aspect('equal')
#        axIma.scatter(table['X_IMAGE'], table['Y_IMAGE'] - decalage, 
#                      s=200, c='none', edgecolors='white', alpha=0.50)#;#        plt.colorbar()
#        #axIma.set_xlim((0, lx))
#        axIma.set_ylim((ly, 0))
#        divider = make_axes_locatable(axIma)
#        axx = divider.append_axes("top", size=1.05, pad=0.2, sharex=axIma)
#        axy = divider.append_axes("right", size=1.05, pad=0.2, sharey=axIma)
#        table.sort('xcentroid')
#        axx.plot(table['xcentroid'], table['FWHM_x'], 'or', label='FWHM_x' )
#        axx.plot(table['xcentroid'], table['FWHM_y'], 'og', label='FWHM_y' )
#        x=np.linspace(0,1300,1000)
#        if len(table['xcentroid']>2):
#            a = np.polyfit(np.concatenate((table['xcentroid']  ,table['xcentroid'])), np.concatenate((table['FWHM_x'] ,table['FWHM_y'])),1)
#            b = np.polyfit(np.concatenate((table['xcentroid']  ,table['xcentroid'])), np.concatenate((table['FWHM_x'] ,table['FWHM_y'])),2)
#            axx.plot(x ,a[1] +x* a[0], '--', label='Linear regression' )    
#            axx.plot(x,b[2] +x * b[1]+x**2 * b[0], '--', label='Pol regression' )    
#        table.sort('ycentroid')
#        axy.plot(table['FWHM_x'], table['ycentroid'] - decalage, 'or')
#        axy.plot(table['FWHM_y'], table['ycentroid'] - decalage, 'og')
#        if len(table['xcentroid']>2):
#            c = np.polyfit(np.concatenate((table['ycentroid']  ,table['ycentroid'])),np.concatenate((table['FWHM_x'] ,table['FWHM_y'])),1)
#            d = np.polyfit(np.concatenate((table['ycentroid']  ,table['ycentroid'])),np.concatenate((table['FWHM_x'] ,table['FWHM_y'])),2)
#            axy.plot(c[1] + (x- decalage) * c[0], (x- decalage), '--')    
#            axy.plot(d[2] + (x- decalage) * d[1]+(x- decalage)**2 * d[0], (x- decalage), '--') 
#    #        fig.legend((lx,ly,fx,fy), ('FWHM_x','FWHM_y','Lin reg','Pol reg'), loc = 'upper right')
#    #        fig.legend((lx), ('FWHM_x'), loc = (0.8,0.8))
#        try:
#            plt.figtext(.75, .78,'Tag =  %0.1f\nPos Angle = %0.1f deg \nAxis1 = %0.3f mm \nAxis2 = %0.3f mm \nAxis3 = %0.3f mm\nExposure = %0.0f ms\nAzimuth = %0.1f deg\nElevation = %0.1f deg' % (self.header['IMGTAG'],self.header['ROTENC'],self.header['LINAENC'],self.header['LINBENC'],self.header['LINCENC'],self.header['EXPOSURE'],self.header['AZ'],self.header['AZ']), fontsize=8, color = 'black')
#        except KeyError:
#            print('No header')
#        axx.grid()
#        axy.grid()
#        axx.set_ylim((0, 20))
#        axy.set_xlim((0, 20))
##        axx.set_xlim(axIma.get_xlim())
##        axy.set_ylim(axIma.get_ylim())
#        fig.savefig(self.all.filename()[:-5], dpi=100, figsize=(10,10))
#        print('image saved')




    def plot_FWHM_guider(self, name, figsize=10, cut=[50, 99.9]):
        """Function that plot a guider image of the Image class nicely 
        """
        n=500
        table = self.table
        if type(self.line) is int:
            image = self.all[0].data[self.line-n:self.line+n, :] 
            decalage = self.line - n
        else:
            image = self.image
            #self.line = 0.0
            decalage = 0
        ly, lx = image.shape

        fs = np.array([lx,ly], dtype=float) 
        fs *= figsize/fs.max()
        fig = plt.figure(figsize=fs)
        axIma = plt.subplot(111)
        norm = ImageNormalize(image, interval=AsymmetricPercentileInterval(cut[0], cut[1]), stretch=AsinhStretch())   
        axIma.imshow(image, cmap = 'cubehelix', origin='lower', norm=norm)
        axIma.set_aspect('equal')
        axIma.scatter(table['X_IMAGE'], table['Y_IMAGE'] - decalage, 
                      s=200, c='none', edgecolors='white', alpha=0.50)#;#        plt.colorbar()
        #axIma.set_xlim((0, lx))
        #axIma.set_ylim((0, ly))
        divider = make_axes_locatable(axIma)
        axx = divider.append_axes("top", size=1.05, pad=0.2, sharex=axIma)
        axy = divider.append_axes("right", size=1.05, pad=0.2, sharey=axIma)
        table.sort('xcentroid')
        axx.plot(table['xcentroid'], table['FWHM_x'], 'or', label='FWHM_x' )
        axx.plot(table['xcentroid'], table['FWHM_y'], 'og', label='FWHM_y' )
        x=np.linspace(0,1300,1000)
        if len(table['xcentroid']>2):
            a = np.polyfit(np.concatenate((table['xcentroid']  ,table['xcentroid'])), np.concatenate((table['FWHM_x'] ,table['FWHM_y'])),1)
            b = np.polyfit(np.concatenate((table['xcentroid']  ,table['xcentroid'])), np.concatenate((table['FWHM_x'] ,table['FWHM_y'])),2)
            axx.plot(x ,a[1] +x* a[0], '--', label='Linear regression' )    
            axx.plot(x,b[2] +x * b[1]+x**2 * b[0], '--', label='Pol regression' )    
        table.sort('ycentroid')
        axy.plot(table['FWHM_x'], table['ycentroid'] - decalage, 'or')
        axy.plot(table['FWHM_y'], table['ycentroid'] - decalage, 'og')
        if len(table['xcentroid']>2):
            c = np.polyfit(np.concatenate((table['ycentroid']  ,table['ycentroid'])),np.concatenate((table['FWHM_x'] ,table['FWHM_y'])),1)
            d = np.polyfit(np.concatenate((table['ycentroid']  ,table['ycentroid'])),np.concatenate((table['FWHM_x'] ,table['FWHM_y'])),2)
            axy.plot(c[1] + (x- decalage) * c[0], (x- decalage), '--')    
            axy.plot(d[2] + (x- decalage) * d[1]+(x- decalage)**2 * d[0], (x- decalage), '--') 
    #        fig.legend((lx,ly,fx,fy), ('FWHM_x','FWHM_y','Lin reg','Pol reg'), loc = 'upper right')
    #        fig.legend((lx), ('FWHM_x'), loc = (0.8,0.8))
        try:
            plt.figtext(.75, .78,'Tag =  %0.1f\nPos Angle = %0.1f deg \nAxis1 = %0.3f mm \nAxis2 = %0.3f mm \nAxis3 = %0.3f mm\nExposure = %0.0f ms\nAzimuth = %0.1f deg\nElevation = %0.1f deg' % (self.header['IMGTAG'],self.header['ROTENC'],self.header['LINAENC'],self.header['LINBENC'],self.header['LINCENC'],self.header['EXPOSURE'],self.header['AZ'],self.header['AZ']), fontsize=8, color = 'black')
        except KeyError:
            print('No header')
        axx.grid()
        axy.grid()
        axx.set_ylim((0, 20))
        axy.set_xlim((0, 20))
#        axx.set_xlim(axIma.get_xlim())
#        axy.set_ylim(axIma.get_ylim())
        fig.savefig(self.all.filename()[:-5], dpi=100, figsize=(10,10))
        print('image saved')


    def plot_FWHM_det(self, name=None, figsize=16, cut=[50, 99.9]):
        """Function that plot a detector image of the Image class nicely 
        """
        n =500
        table = self.table
        if self.line is None:
            image = self.image
            #self.line = 0.0
            decalage = 0
        else:
            image = self.all[0].data[self.line-n:self.line+n, 1064:2130] 
            decalage = self.line - n#fsfsfs
        ly,lx = image.shape # ly is vertical size

        if self.windowing:
            try:
                colors = [self.wcolors[w] for w in table['wavelength']]
            except KeyError:
                colors = 'white'
                pass
        else:
            colors = 'white'

        fs = np.array([lx,ly], dtype=float) 
        fs *= figsize/fs.max()
        fig = plt.figure(figsize=fs)

        axIma = plt.subplot(111)
        norm = ImageNormalize(image, interval=AsymmetricPercentileInterval(cut[0], cut[1]), stretch=AsinhStretch())   
        axIma.imshow(image, cmap='cubehelix', norm=norm, origin='lower')
        axIma.set_aspect(1.)
        axIma.scatter(table['X_IMAGE'] - self._DFLT_Yinf, table['Y_IMAGE'] - self._DFLT_Xinf - decalage, 
                      s=200, c='none', edgecolors=colors, alpha=0.70)#;#        plt.colorbar()
        axIma.scatter(table['xcentroid'], table['ycentroid'] - decalage, 
                      s=20, c='none', edgecolors=colors, alpha=0.70)#;#        plt.colorbar()
        axIma.set_xlim((0, lx))
        axIma.set_ylim((0, ly))
        # create new axes on the right and on the top of the current axes.
        divider = make_axes_locatable(axIma)
        axx = divider.append_axes("top", size=1.05, pad=0.2, sharex=axIma)
        axy = divider.append_axes("right", size=1.05, pad=0.2, sharey=axIma)
        try:
            axy.plot(table[table['wavelength']==self.w[0]]['FWHM_x'], table[table['wavelength']==self.w[0]]['Y_IMAGE'] - decalage, 'bo', label='FWHM_x' )
            axy.plot(table[table['wavelength']==self.w[1]]['FWHM_x'], table[table['wavelength']==self.w[1]]['Y_IMAGE'] - decalage, 'go', label='FWHM_x' )
            axy.plot(table[table['wavelength']==self.w[2]]['FWHM_x'], table[table['wavelength']==self.w[2]]['Y_IMAGE'] - decalage, 'ro', label='FWHM_x' )
            axx.plot(table[table['wavelength']==self.w[0]]['X_IMAGE'] - self._DFLT_Yinf, table[table['wavelength']==self.w[0]]['FWHM_y'], 'bo', label='FWHM_y' )
            axx.plot(table[table['wavelength']==self.w[1]]['X_IMAGE'] - self._DFLT_Yinf, table[table['wavelength']==self.w[1]]['FWHM_y'], 'go', label='FWHM_y' )
            axx.plot(table[table['wavelength']==self.w[2]]['X_IMAGE'] - self._DFLT_Yinf, table[table['wavelength']==self.w[2]]['FWHM_y'], 'ro', label='FWHM_y' )
            #axx.set_xlim((0, lx))
            axIma.scatter(self.predictedSpots[self.predictedSpots['w']==self.w[0]]['x'] - self._DFLT_Yinf - decalage, self.predictedSpots[self.predictedSpots['w']==self.w[0]]['y'], s=10, c='blue', marker='+')#;#        plt.colorbar()
            axIma.scatter(self.predictedSpots[self.predictedSpots['w']==self.w[1]]['x'] - self._DFLT_Yinf - decalage, self.predictedSpots[self.predictedSpots['w']==self.w[1]]['y'], s=10, c='green',marker='+')#;#        plt.colorbar()
            axIma.scatter(self.predictedSpots[self.predictedSpots['w']==self.w[2]]['x'] - self._DFLT_Yinf - decalage, self.predictedSpots[self.predictedSpots['w']==self.w[2]]['y'], s=10, c='red', marker='+')#;#        plt.colorbar()
        except AttributeError:
            print('No windowing')           
            axy.plot(table[table['wavelength']>=0]['FWHM_x'], table[table['wavelength']>=0]['Y_IMAGE'] - decalage, 'o', color='black', label='FWHM_x' )
            axx.plot(table[table['wavelength']>=0]['X_IMAGE'] - self._DFLT_Yinf, table[table['wavelength']>=0]['FWHM_y'], 'o', color='black', label='FWHM_y' )
        axx.set_ylabel('FWHM_y')
        axy.set_xlabel('FWHM_x')
            #axx.set_xlim((0, lx))
     #   axx.plot(table['ycentroid'],table['FWHM_y'],'o',label = 'FWHM_y' )
        axx.grid(); axy.grid()
        axx.set_ylim(((table['FWHM_y'].min()), (table['FWHM_x'].max())))
        axy.set_xlim(((table['FWHM_x'].min()), (table['FWHM_x'].max())))
        axx.set_xlim((0,lx))
        axy.set_ylim((0, ly))
        table.sort('xcentroid')
        x = np.linspace(0,2200,1000)
        if len(table['xcentroid']>2):
            a = np.polyfit(np.concatenate((table['ycentroid'], table['ycentroid'])), np.concatenate((table['FWHM_x'], table['FWHM_y'])), 1)
            b = np.polyfit(np.concatenate((table['ycentroid'], table['ycentroid'])), np.concatenate((table['FWHM_x'], table['FWHM_y'])), 2)
            axx.plot(x, a[1] +x* a[0], '--', label = 'Linear regression' )    
            axx.plot(x, b[2] +x * b[1]+x**2 * b[0], '--', label='Pol regression' )    
        table.sort('ycentroid')
        if len(table['xcentroid']>2):
            c = np.polyfit(np.concatenate((table['xcentroid'] , table['xcentroid'])), np.concatenate((table['FWHM_x'], table['FWHM_y'])), 1)
            d = np.polyfit(np.concatenate((table['xcentroid'] , table['xcentroid'])), np.concatenate((table['FWHM_x'], table['FWHM_y'])), 2)
            axy.plot(c[1] +(x- decalage) * c[0], (x- decalage), '--')    
            axy.plot(d[2] +(x- decalage) * d[1] + (x- decalage)**2 * d[0], (x- decalage), '--') 

        if name is None:
            name = self.all.filename()[:-5]
        fig.savefig(name,dpi=100,figsize=(10,10))



    def plotFWHM2D(self,save=True, allw=False, a=7,t=25):
        table = self.table[(self.table['FWHM_x']<20) & (self.table['FWHM_y']<20)].copy()
#        plt.figure()
        fig = plt.figure(figsize = (6,9))
        ax = fig.add_subplot(1, 1, 1)
        ax.axis('equal')
        xmajor_ticks = np.arange(800, 2200, 200)
        xminor_ticks = np.arange(800, 2200, 50)
        ymajor_ticks = np.arange(0, 2000, 200)
        yminor_ticks = np.arange(0, 2000, 50)
        ax.set_xticks(xmajor_ticks)
        ax.set_xticks(xminor_ticks, minor=True)
        ax.set_yticks(ymajor_ticks)
        ax.set_yticks(yminor_ticks, minor=True)
        #ax.grid(which='both')
        # Or if you want different settings for the grids:
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)
        plt.errorbar(1000,0,yerr = a*5,xerr = a*5, ecolor='orange',linewidth=5,label='5pix FWHM PSF')
        plt.legend()
        if allw:
            for w,c in zip([0.20255,0.20619,0.21382,0.0,-1.0],['b','g','r','black','black']):
                plt.errorbar(table[table['wavelength']==w]['X_IMAGE'], table[table['wavelength']==w]['Y_IMAGE'], yerr=a*table[table['wavelength']==w]['FWHM_y'], xerr=a*table[table['wavelength']==w]['FWHM_x'], fmt='+',ecolor=c,linewidth=5)
        else: 
            for w,c in zip([0.20255,0.20619,0.21382,0.0],['b','g','r']):
                plt.errorbar(table[table['wavelength']==w]['X_IMAGE'], table[table['wavelength']==w]['Y_IMAGE'], yerr=a*table[table['wavelength']==w]['FWHM_y'], xerr=a*table[table['wavelength']==w]['FWHM_x'], fmt='+',ecolor=c,linewidth=5)

        plt.title(self.mask + ' ' + np.str(self.pa) + ' FWHM_xy')
#        name = self.all.filename()[:-5] + 'FWHMxy'
        name = self.mask + '_' + np.str(np.int(self.pa)) + '_FWHMxy'
        plt.savefig(os.path.dirname(self.all.filename()) + '/' + name,dpi=100)
        return



    def plotFWHMEllipse(self,save=True, allw=False, a=10,t=25,n=20,name=None):
        if not allw:
            print('Taking only attributed wavelength')
            table = self.table[(self.table['FWHM_x']<t) & (self.table['FWHM_y']<t) & (self.table['wavelength']>0.1)].copy()
        else:
            table = self.table[(self.table['FWHM_x']<t) & (self.table['FWHM_y']<t)].copy()
        X,Y,Z = fit_quadratic_curve(table['X_IMAGE'],table['Y_IMAGE'],table['FWHM_x']**2 + table['FWHM_y']**2,n)
        x = table['Y_IMAGE']
        y = - table['X_IMAGE']
#        ells = [Ellipse(xy=[x[i],y[i]],
##                        width=a*table[i]['FWHM_y'], height=a*table[i]['FWHM_x'],
#                        width=a*np.sqrt(table[i]['FWHM_y']**2 + table[i]['FWHM_x']**2), height=a*np.sqrt(table[i]['FWHM_y']**2 + table[i]['FWHM_x']**2),
#                        angle=0) for i in range(len(table))]
        #for w,c in zip([0.20261371,0.20626604,0.21392365,0.0,-1.0],['b','g','r','black','black']):
        c=[]
        for i in range(len(table)):
            if table[i]['wavelength'] == 0.20255:
                c.append('b')
            if table[i]['wavelength'] == 0.20619:
                c.append('g')
            if table[i]['wavelength'] == 0.21382:
                c.append('r')
            if (table[i]['wavelength']==0.0) or (table[i]['wavelength']==-1.0):
                c.append('black')
        
        fig = plt.figure(figsize=(9,6))
        ax = fig.add_subplot(111, aspect='equal')
        ax.get_yaxis().set_ticklabels([])
        ax.scatter(Y,-X,s=a*Z,alpha=0.15)
        xmajor_ticks = np.arange(-3000,3000, 200)
        xminor_ticks = np.arange(-3000,3000, 50)
        ymajor_ticks = np.arange(-3000,3000, 200)
        yminor_ticks = np.arange(-3000, 3000, 50)
#        xmajor_ticks = np.arange(x.min(),x.max(), 200)
#        xminor_ticks = np.arange(x.min(),x.max(), 50)
#        ymajor_ticks = np.arange(y.min(),y.max(), 200)
#        yminor_ticks = np.arange(y.min(),y.max(), 50)
        ax.set_xticks(xmajor_ticks)
        ax.set_xticks(xminor_ticks, minor=True)
        ax.set_yticks(ymajor_ticks)
        ax.set_yticks(yminor_ticks, minor=True)
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)
        ax.get_xaxis().set_ticklabels([])

        xx,yy = (x.max()+x.min())/2,(y.max()+y.min())/2
        #print((xx,yy))
        #e = Ellipse(xy=[xx,yy],width=a*5, height=a*5,angle=0)       
        #ax.add_patch(e)
        #e.set(clip_box=ax.bbox, alpha=1, facecolor=[0,0,0],label='5pix FWHM PSF')
        #ax.scatter()
        #ax.axis('off')
#        for i,e in enumerate(ells):
#            ax.add_artist(e)
#            e.set_clip_box(ax.bbox)
#            e.set_alpha(0.7)
#            e.set_facecolor(c[i])
        ax.scatter(x.min()-100,y.min()-100,s=a*(3**2 + 3**2),alpha=0.7,c='black',label=r'3 pix FWHM PSF')            
        ax.legend(loc = 'upper right')

        ax.scatter(x,y,s=a*(table['FWHM_y']**2 + table['FWHM_x']**2),  c='none', edgecolors=c ,alpha=0.7,linewidth=3)
        ax.set_ylim(y.min()-50, y.max()+50)
        ax.set_xlim(x.min()-50,x.max()+50)
        if name:
            name = name
        else:
            name = self.mask + '_' + np.str(np.int(self.pa)) + '_FWHMxy_Ellipse'
        plt.title(name)
        plt.savefig(os.path.dirname(self.all.filename()) + '/' +name,dpi=100)
        plt.show()
        
    def plotFWHMEllipse2(self,save=True, allw=False, a=10,t=25):
        table = self.table[(self.table['FWHM_x']<t) & (self.table['FWHM_y']<t)].copy()
        ells = [Ellipse(xy=[table[i]['X_IMAGE'],table[i]['Y_IMAGE']],
                        width=a*table[i]['FWHM_x'], height=a*table[i]['FWHM_y'],
                        angle=0) for i in range(len(table))]
        #for w,c in zip([0.20261371,0.20626604,0.21392365,0.0,-1.0],['b','g','r','black','black']):
        c=[]
        for i in range(len(table)):
            if table[i]['wavelength'] == 0.20255:
                c.append('b')
            if table[i]['wavelength'] == 0.20619:
                c.append('g')
            if table[i]['wavelength'] == 0.21382:
                c.append('r')
            if (table[i]['wavelength']==0.0) or (table[i]['wavelength']==-1.0):
                c.append('black')
        
        fig = plt.figure(figsize=(6,9))
        ax = fig.add_subplot(111, aspect='equal')
        xmajor_ticks = np.arange(800, 2200, 200)
        xminor_ticks = np.arange(800, 2200, 50)
        ymajor_ticks = np.arange(0, 2000, 200)
        yminor_ticks = np.arange(0, 2000, 50)
        ax.set_xticks(xmajor_ticks)
        ax.set_xticks(xminor_ticks, minor=True)
        ax.set_yticks(ymajor_ticks)
        ax.set_yticks(yminor_ticks, minor=True)
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)
        for i,e in enumerate(ells):
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.7)
            e.set_facecolor(c[i])
        ax.set_xlim(table['X_IMAGE'].min()-50, table['X_IMAGE'].max()+50)
        ax.set_ylim(table['Y_IMAGE'].min()-50, table['Y_IMAGE'].max()+50)
        name = self.mask + '_' + np.str(np.int(self.pa)) + '_FWHMxy_Ellipse'
        plt.title(self.mask + ' ' + np.str(np.int(self.pa)) + ' FWHMxy Ellipse')
        plt.savefig(os.path.dirname(self.all.filename()) + '/' +name,dpi=100)
        plt.show()

#grid81.plotFWHMEllipse()


if __name__ == '__main__':
    pass
#    files = [] 
#    path = os.path.dirname(os.path.realpath(__file__))
#    files = glob.glob(path + '*.FIT')[::-1]
#    images = []
#    for i in range(int(len(files))):
#        print(('{} over {}'.format(i,len(files))))
#        im = Image(filename = files[0], quick=True, figsize=12,)#im = Image(filename = files[i], plot = False, stack = True, stack_image = False, py = True, subDark = False, verbose = True, quick = True, Type = 'detector', Fits = True)
#        images.append(im)
#    
#    filename = '/Users/Vincent/Nextcloud/FIREBALL/Tests-at-FortSumner/170901_perf3.0/170901/focus/image000091stack.fits'
#    im = Focus(filename = filename, quick=True, figsize=12,windowing=True, mask='F1')#im = Image(filename = files[i], plot = False, stack = True, stack_image = False, py = True, subDark = False, verbose = True, quick = True, Type = 'detector', Fits = True)

#
#    filename = '/Users/Vincent/Nextcloud/FIREBALL/CalibrationSW_Procedures/Notebook_Test/TestImages/grid.fits'
#    im = Focus(filename = filename, quick=True,reversex=True)#, figsize=12,windowing=True, mask='grid' , plot=False,sources='holes')

#    filename = '/Users/Vincent/Nextcloud/FIREBALL/CalibrationSW_Procedures/Notebook_Test/TestImages/grid.fits'
#    im = Focus(filename = filename, quick=True, figsize=12,windowing=True, mask='grid' , plot=False,sources='holes',date=5)
#    x,y,disp = im.CalculateDispersion()
#   interpDispersion(x,y,disp)

    #filename = '/Users/Vincent/Nextcloud/FIREBALL/CalibrationSW_Procedures/Notebook_Test/TestImages/F1.fits'
    #im = Focus(filename = filename, quick=True, figsize=12,windowing=True, mask='grid', shape='slits2D',plot=True)
#
#    
#    path1 = '/Users/Vincent/NextcloudBackUp/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/'
#    F2 = Focus(filename = path1 + 'image-000025-000034-Zinc-with_dark-161-stack.fits', 
#               quick=False, threshold = [7], fwhm = [9,12.5],
#               HumanSupervision=False,
#    reversex=False, plot=False,source='all',
#    shape='slits',windowing=True, mask='F2', pa=-161,MoreSources=0,peak_threshold=50)


#imshow(xx[:,::-1])
#
#    plot_rp2_convolved_wo_latex(im.image, center=[214,143])
#
#plot_rp2_convolved_wo_latex(newgrid.image, center=[1919-1064,1811])
#plot_rp2_convolved_wo_latex(newgrid.image, center=[1619-1064,1796])

#
#plot_rp2_convolved_wo_latex(newgrid.image, center=[1919-1064,1811])
#
##
##a, b, rp = plot_rp2_convolved_wo_latex(gridFR,G2TDPK center=[1919,1811])
##plot(b,rp)
##plot(b,2*ConvolveDiskGaus2D(b,1/2,2,1))
#
#
#nr = np.bincount(r2.astype(np.int).ravel())
#
#rsurf, rmean,radialprofile, EE = radial_profile_normalized(gridFR, center=[1919,1811])
#plt.plot(rmean[:radius],radialprofile[:radius],'-x')
##plt.plot(arange(40),radialprofile[:radius],'-x')
#
#plt.plot(rsurf[:radius],EE[:radius],'x')
#plot(np.sqrt(cumsum(nr)[:40]/np.pi));plot(rmean[:40])
    
#def defocus(xy, det_offset, det_xslope, det_yslope, a, b):
#x_mask = xy[0]
#x_det = xy[1]
#if fit_yslope:
#    y_det = xy[2]
#
## defocus from mask
#tilt = 4.
#sin_tilt = tilt * np.pi / 180.            
#defoc = (x_mask - xmask0) * sin_tilt
#
## defocus from det
#defoc += det_offset + x_det * det_pix_size * det_xslope # take into account lambda shift
#if fit_yslope:
#    defoc += y_det * det_pix_size * det_yslope # not degenerate with a BECAUSE with distortion y_det not proportional to x_mask 
#
##fwhm from defocus ()
#fwhm = a*defoc**2 + b
#return fwhm
#
#
#
#fig = plt.figure(figsize=(6,12))
#ax = fig.add_subplot(121, projection='3d')
#ax.set_title('defocus lambda={:.4f}'.format(w))
#ax.set_xlabel('X mask mm')
#ax.set_ylabel('y mask mm')
#ax.set_zlabel('FWHM size pix')
#ax.scatter(focus.xmask[sid], focus.ymask[sid], zdata, c=focus.wcolors[w])
#gx, gy = np.meshgrid(np.linspace(-12, 12, 20), np.linspace(-6, 6, 20))
#gxy_det = im.apply_map_mask_detector(gx, gy, np.full_like(gx, 0.214))
#gx_det, gy_det = gxy_det[0], gxy_det[1]
#if fit_yslope:
#    gz = defocus(np.array([gx, gx_det, gy_det]), *popt)
#        
#    path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/Autocoll/Focus_06122018/'
#    for file in glob.glob(path+('*.fits'))[7:]:
#        test = Focus(filename = file, threshold = [9], fwhm = [5,7,9,12.5,15, 17],
#                     quick=False, reversex=False, plot=True, figsize=12,
#                     shape='gaussian', min_fwhm=2.5, cut=[50,99],peak_threshold=None )#, HumanSupervision=True)
#
#
#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/DetectorShim-20180526/'
#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/Autocoll/Focus_06122018/'
#quicklook(path=path,ptype='guider')
#quicklook(path=path,ptype='detector',mmr=[1,10,2])
#
#print(ds9_targets())
#
#    #### DS9 #####
#    path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/'
##        
#    F1 = Focus(filename = path + 'image-000275-000284-Zinc-with_dark119-stack.fits', 
#               threshold = [7], fwhm = [9,12.5],HumanSupervision=False,
#    quick=False, reversex=False, plot=False,source='Zn',
#    shape='slits',windowing=True, mask='F1', pa=119,peak_threshold=50,min_fwhm=1.5,dist_min=30)
#
#
##    
#    F3 = Focus(filename = path + 'image-000075-000084-Zinc-with_dark-121-stack.fits', 
#               threshold = [7], fwhm = [9,12.5],HumanSupervision=False,
#    quick=False, reversex=False, plot=False,source='Zn',
#    shape='slits',windowing=True, mask='F3', pa=-121,peak_threshold=50,min_fwhm=1.5)
#  
#    
#        
#    F2 = Focus(filename = path + 'image-000025-000034-Zinc-with_dark-161-stack.fits', 
#               threshold = [7], fwhm = [9,12.5],HumanSupervision=False,
#    quick=False, reversex=False, plot=False,source='Zn',
#    shape='slits',windowing=True, mask='F2', pa=-161,peak_threshold=50,min_fwhm=1.5)
#  
#    
#        
#    F4 = Focus(filename = path + 'image-000325-000334-Zinc-with_dark159-stack.fits', 
#               threshold = [7], fwhm = [9,12.5],HumanSupervision=False,
#    quick=False, reversex=False, plot=False,source='Zn',
#    shape='slits',windowing=True, mask='F4', pa=159,peak_threshold=50,min_fwhm=1.5)
#  
    
#
#a = Focus(filename = '/Users/Vincent/Nextcloud/FIREBALL/Tests-at-FortSumner/170923_SC_GUI03/Detector/170924/image000001.fits')
#a.table.sort('peak')
#scatter(a.table['X_IMAGE'],a.table['Y_IMAGE'],s=a.table['peak'])
#n1=20
#n2=2
#xcenter,ycenter = int(a.table['xcentroid'][1]), int(a.table['ycentroid'][1])
##plt.plot(a.image[ycenter,xcenter-n:xcenter+n])
#plt.plot(a.image[ycenter-n1:ycenter+n1,xcenter-n2:xcenter+n2].sum(axis=1))
#
#
#
#xcenter,ycenter = int(b.table['X_IMAGE'][0]), int(b.table['Y_IMAGE'][0])
##plt.plot(a.image[ycenter,xcenter-n:xcenter+n])
#plt.plot(b.image[ycenter-n1:ycenter+n1,xcenter-n2:xcenter+n2].sum(axis=1))
#
#
#d.set('plot new stdin ')
#d.set('plot new name profile "Profile" "X image" "profile" xy')
#d.set('plot profile xy [[1],[1]]')
#d.set('plot profile xy 1 1')


#pagesetup scale 
#a = Focus(filename='/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/Autocoll/Detector/1/image000379.fits',plot=True,HumanSupervision=True)

#    path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/AIT-Optical-FTS-2018-Flight/XYCalibration/'
#        
#    F1 = Focus(filename = path + 'StackedImage_24-43-NoDark.fits', 
#               threshold = [7], fwhm = [9,12.5],HumanSupervision=False,
#    quick=False, reversex=False, plot=False,source='Zn',
#    shape='slits',windowing=True, mask='F1', pa=119,peak_threshold=50,min_fwhm=1.5,dist_min=30)
#
#    F4 = Focus(filename = path + 'StackedImage_44-63-NoDark.fits', 
#               threshold = [7], fwhm = [9,12.5],HumanSupervision=False,
#    quick=False, reversex=False, plot=False,source='Zn',
#    shape='slits',windowing=True, mask='F4', pa=159,peak_threshold=50,min_fwhm=1.5,dist_min=30)
#f1 = Focus(filename = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180823/xycalib/image000448.fits',HumanSupervision=False,source='Zn',shape='holes',windowing=False)
#f1 = Focus(filename = '/Users/Vincent/Nextcloud/Work/FTS2018_FLIGHT/test/InstrumentCentering_180819/F3_-121/3_20_12m/stack11398095.fits',HumanSupervision=False,source='Zn',shape='holes',windowing=False)
#files = glob.glob('/Users/Vincent/Nextcloud/Work/FLIGHT_DATA/Diffuse for flight/*.fits')
#
#for file in files:
#    image = fits.open(file)[0]
#    image.data = image.data[:,1053:2133]
#    image.writeto(file[:-5] + '_cut.fits',overwrite=True)
#
#













