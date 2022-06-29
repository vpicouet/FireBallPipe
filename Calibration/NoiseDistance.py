#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 18:59:57 2017

@author: vincent
"""

from __future__ import division
from astropy.io import fits
import astropy
import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import dblquad
from mpl_toolkits.mplot3d import axes3d  ### ne pas ENLEVER : est utilisé
import glob
from scipy import ndimage
from astropy.table import Table
from astropy.table import vstack, Table
from astropy.table import Column
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
import os, sys
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import argparse
from astropy.table import vstack, hstack, Table
from photutils import make_source_mask
from scipy import misc
from skimage import color
from matplotlib.ticker import NullFormatter



def plot_FWHM(table,image,name):
    n=1.4
    nullfmt = NullFormatter()         # no labels
    left, width = 0.1, 1#0.65
    bottom, height = 0.1, 1080./1280#0.65
    bottom_h = left_h = left + width + 0.02
    rect_scatter = [left/n, bottom/n, width/n, height/n]
    rect_histx = [left/n,( bottom_h-0.15)/n, width/n, 0.2/n]
    rect_histy = [left_h/n, bottom/n, 0.2/n, height/n]
    fig = plt.figure(figsize=(7, 7))
    axim = plt.axes(rect_scatter)
    axx = plt.axes(rect_histx)
    axy = plt.axes(rect_histy)

    axx.xaxis.set_major_formatter(nullfmt)
    axy.yaxis.set_major_formatter(nullfmt)
    axim.imshow(im.image - im.median)
    axim.scatter(table['X_IMAGE'],table['Y_IMAGE'] , s=200, c='white', alpha=0.15)#;#        plt.colorbar()
    axim.set_xlim((0, im.image.shape[1]))
    axim.set_ylim((0, im.image.shape[0]))
    table.sort('xcentroid')
    x=np.linspace(0,1300,1000)
    if len(table['xcentroid']>2):
        a = np.polyfit(np.concatenate((table['xcentroid']  ,table['xcentroid'])),np.concatenate((table['FWHMx'] ,table['FWHMy'])),1)
        b = np.polyfit(np.concatenate((table['xcentroid']  ,table['xcentroid'])),np.concatenate((table['FWHMx'] ,table['FWHMy'])),2)
        axx.plot(x ,a[1] +x* a[0], '--',label = 'Linear regression' )    
        axx.plot(x,b[2] +x * b[1]+x**2 * b[0], '--',label = 'Pol regression' )    
    lx = axx.plot(table['xcentroid'],table['FWHMx'],'o',label = 'FWHMx' )
    ly = axx.plot(table['xcentroid'],table['FWHMy'],'o',label = 'FWHMy' )
    table.sort('ycentroid')
    if len(table['xcentroid']>2):
        c = np.polyfit(np.concatenate((table['ycentroid']  ,table['ycentroid'])),np.concatenate((table['FWHMx'] ,table['FWHMy'])),1)
        d = np.polyfit(np.concatenate((table['ycentroid']  ,table['ycentroid'])),np.concatenate((table['FWHMx'] ,table['FWHMy'])),2)
        fx = axy.plot(c[1] +x* c[0],x, '--')    
        fy = axy.plot(d[2] +x * d[1]+x**2 * d[0],x, '--') 
#        fig.legend((lx,ly,fx,fy), ('FWHMx','FWHMy','Lin reg','Pol reg'), loc = 'upper right')
#        fig.legend((lx), ('FWHMx'), loc = (0.8,0.8))
    plt.figtext(0.8, 0.69,'Pos Angle = %0.1f deg \nAxis1 = %0.1f mm \nAxis2 = %0.1f mm \nAxis3 = %0.1f mm\nExposure = %0.0f ms\nAzimuth = %0.1f deg\nElevation = %0.1f deg' % (im.header['ROTENC'],im.header['LINAENC'],im.header['LINBENC'],im.header['LINCENC'],im.header['EXPOSURE'],im.header['AZ'],im.header['AZ']), fontsize=8, color = 'black')

    axy.plot(table['FWHMx'],table['ycentroid'],'o')
    axy.plot(table['FWHMy'],table['ycentroid'],'o')
    axx.grid()
    axy.grid()
    axx.set_ylim((4, 20))
    axy.set_xlim((4, 20))
    axx.set_xlim(axim.get_xlim())
    axy.set_ylim(axim.get_ylim())
    fig.show()
    axim.set_xlabel(name +'    ' +  im.header['Date']) 
    fig.savefig(name +'_'+im.header['Date'][11:13]+'h'+im.header['Date'][14:16]+'m'+im.header['Date'][17:19],dpi=100,figsize=(10,10))
    plt.show()


def commandLine():
    image = 'image0000{}'

     # help
    parser = argparse.ArgumentParser(description="Analyze the spectral and spatial resolution of an image", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # options
    parser.add_argument('number', help='Number of the image to be analysed', type = int)
    parser.add_argument('py',help='True for using photutils from astropy, false for using sextractor', default = False, type = bool)
    parser.add_argument('--plot', help='plot', type = bool, default = False)
    parser.add_argument('--stack', help='True to stack the spot along spatial and spectral resolution',default=True, type = bool)
    parser.add_argument('--stack_image',help='True for stacking the three images +1 & -1', default=False, type = bool)

     # parse
#    args = parser.parse_args()
#    print 'py = ', args.py
#    print 'plot = ', args.plot
#    print 'stack = ', args.stack
#    print 'stack_image = ', args.stack_image

    # start main funciton
    try:
        fits.delval(image.format(args.number) + '.fits','NAXIS3')
    except KeyError:
        print 'ok'
    image_n = Image(number = args.number, plot = args.plot, stack = args.stack, stack_image = args.stack_image, py = args.py)


    
def delete_doublons(sources, dist):
    try:
        sources['doublons'] = 0
        for i in range(len(sources)):
            a = distance(sources[sources['doublons']==0]['xcentroid'],sources[sources['doublons']==0]['ycentroid'],sources['xcentroid'][i],sources['ycentroid'][i]) > dist
            a = list(1*a)
            a.remove(0)
            if np.mean(a)<1:
                sources['doublons'][i]=1
        return sources[sources['doublons']==0]
    except TypeError:
        print 'no source detected'
#        quit()
        


def distance(x1,y1,x2,y2):
    return np.sqrt(np.square(x1-x2)+np.square(y1-y2))


def extraction(a,b, rayon, photo):
    extension = '.fits'
    photo = '{}'.format(photo)
    imgfits = fits.open(photo + extension)
    img = imgfits[0].data
    name = img[b - rayon : b + rayon ,a - rayon :a + rayon]
    return name

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, offset):
    xo = float(xo)
    yo = float(yo)
    A = amplitude/(2*math.pi*sigma_x*sigma_y)    
    g = offset + A * np.exp( - 0.5*(((x-xo)/sigma_x)**2) - 0.5*(((y-yo)/sigma_y)**2))
    return g.ravel()

def Gaussian(x, amplitude, xo, sigma_x, offset):
    xo = float(xo)
    A = amplitude/(sigma_x * np.sqrt(2*math.pi))    
    g = offset + A * np.exp( - 0.5*(((x-xo)/sigma_x)**2))
    return g.ravel()




def Gaussian1(x, amplitude, xo, sigma_x):
    xo = float(xo)
    A = amplitude/(sigma_x * np.sqrt(2*math.pi))    
    g = A * np.exp( - 0.5*(((x-xo)/sigma_x)**2))
    return g.ravel()

def Normalization(amplitude, xo, yo, sigma_x, sigma_y , offset):
    B = dblquad(lambda x, y: twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, 0), -float('inf'), +float('inf'), lambda x: -float('inf'), lambda x: float('inf'))       
    return (amplitude/B[0], xo, yo, sigma_x, sigma_y , 0)


def Normalization(amplitude, xo, sigma_x):
    B = dblquad(lambda x: Gaussian1((x), amplitude, xo, sigma_x), -float('inf'), +float('inf'))       
    return (amplitude/B[0], xo, sigma_x)


class Image(object):

    #static members
    _DFLT_Xinf   = 0#1000
    _DFLT_Xsup   = 2069#2100
    _DFLT_Yinf   = 1064
    _DFLT_Ysup  = 2130#2069
    _DFLT_extension  = '.fits'
    _DFLT_catalog  = 'cat'
    _DFLT_computeImage  = False
    
        
#------------------------------------------------------------------------------#
# Special Methods                                                              #
#------------------------------------------------------------------------------#                                                 
    def __init__(self, **kwargs):
        """
        """
        image = 'image0000{}'
        number = kwargs.get('number', None)

        subDark = kwargs.get('subDark', False)
        quick = kwargs.get('quick', True)
        verbose = kwargs.get('verbose', False)
        plot = kwargs.get('plot', False)

        fname = kwargs.get('filename', None)
        size = kwargs.get('size', [self._DFLT_Xinf,self._DFLT_Xsup,self._DFLT_Yinf,self._DFLT_Ysup])
        plot = kwargs.get('plot', False)
        stack = kwargs.get('stack', True)
        stack_image = kwargs.get('stack_image', True)
        py = kwargs.get('py', False)

	#self.all = pyfits.getdata(fname + self._DFLT_extension)
        all = fits.open(fname)[0]
        self.data = all.data
        self.header = all.header
        ndark = 77
        if stack_image:                            
            self.image = (self.data + fits.open( image.format(number-2)+ self._DFLT_extension)[0].data + fits.open( image.format(number-1)+ self._DFLT_extension)[0].data + fits.open( image.format(number+1) + self._DFLT_extension)[0].data + fits.open( image.format(number+2)+ self._DFLT_extension)[0].data)/5 #[size[0]:size[1],size[2]:size[3]]
            self.all[0].data = self.image                
        else:
            self.image = self.data
        n=0
        self.data = self.image 
        
        if subDark:
            self.dark = (fits.open( image.format(ndark)+ self._DFLT_extension)[0].data + fits.open( image.format(ndark-2)+ self._DFLT_extension)[0].data + fits.open( image.format(ndark-1)+ self._DFLT_extension)[0].data + fits.open( image.format(ndark+1) + self._DFLT_extension)[0].data + fits.open( image.format(ndark+2)+ self._DFLT_extension)[0].data)/5 #[size[0]:size[1],size[2]:size[3]]
            self.image_dark = self.image - self.dark
            self.all[0].data = self.image_dark
            self.all.writeto(image.format(number)+'stack-dark' +self._DFLT_extension, overwrite = True)
            self.image =  self.image[400:,1064:1542]  #[000:,1064:2042] 
            self.dark = self.dark[400:,1064:1542]                         
            self.image = self.image-self.dark


        self.dispersion  = 0
        self.spatialResolution  = 0    
        self.spectralResolution  = 0
        self.pix2sec = 0.93
        self.disp = 0.2167#A/px    # 17.8 A/mm
        
#        plt.figure()
#        plt.imshow(self.image)
#        plt.title(fname)
#        plt.show()
#
#        if 'cat'+ fname + self._DFLT_extension not in glob.glob("*.fits"):
#            try:
#                new_image.writeto(fname + 'm' + self._DFLT_extension)
#                print 'creating new image'
#            except IOError:
#                print 'can not create new image'
#                os.remove(fname + 'm'+ self._DFLT_extension)
#                new_image.writeto(fname + 'm' + self._DFLT_extension) 
#            except astropy.io.fits.verify.VerifyError:
#                print '...'        
#            if self._DFLT_computeImage:
#                print self
#            self.extractSources(**kwargs)
#            
#

        
#        print py
        if py == True:
#            print 'use of py'
#            self.abort = self.extractSourcesPy(**kwargs)
            self.extractBar()
        if py == False:
            print 'Use of sextractor'

#            self.extractSources(**kwargs)

            self.extractCatalog(**kwargs)
#        self.computeDispersionZn(**kwargs)
#        self.computeSpatialResolution()
#        self.computeSpectralResolution()
        
        if plot==3:
            plt.figure()
            plt.imshow(self.image, cmap='hot');plt.colorbar()
            plt.title(fname)
            plt.scatter(self.table['X_IMAGE'],self.table['Y_IMAGE'] , s=60, c='white', alpha=0.4)#;#        plt.colorbar()
#            plt.figtext(.4, .2,'Dispersion = %0.3f A/pix \nFWHM = %0.3f pix \nRes Spatial = %0.3f pix \nRes Spectral = %0.3f' % (self.disp, np.mean(self.table['FWHM_IMAGE']), self.resx, self.resl), fontsize=8, color = 'white')
            plt.show()
            n=60
            
            image = self.image[int(self.table['Y_IMAGE'][-1]-n):int(self.table['Y_IMAGE'][-1]+n),int(self.table['X_IMAGE'][-1]-n):int(self.table['X_IMAGE'][-1]+n)]
            maximum = np.max(self.image)
#            print self.table['X_IMAGE'][-1],self.table['Y_IMAGE'][-1]


#        print ('The FWHM calculated by sextractor are', self.table['FWHM_IMAGE'].data) 
#        self.Fit2D( **kwargs)
        if self.abort==0:
            print 'no source detected'
        else:    
            self.popt1 = self.Fit1D1(**kwargs)[0]
            self.popt2 = self.Fit1D2(**kwargs)[0]

#        self.createCatalog(**kwargs)
#        print ('The average spatiale resolution is ', np.mean(self.resx))
#        print ('The average spectral resolution is ', np.mean(self.resl))

        
        
    def extractSources(self, **kwargs):
        """
        """ 
        fname = kwargs.get('filename', False)
        os.system('sex '+ fname + 'm' + self._DFLT_extension + ' -CATALOG_NAME cat'+ fname + '.fits' + ' -CHECKIMAGE_NAME  ' + fname + '.fits.fits')
        
    def extractCatalog(self, **kwargs):
        fname = kwargs.get('filename', False)
        if fname:
            self.table1 = Table.read(self._DFLT_catalog + fname + self._DFLT_extension)
            self.table = self.table1[(self.table1['FLUX_ISO']<1e5) & (self.table1['FWHM_IMAGE']>2) & (self.table1['FWHM_IMAGE']<25) & (self.table1['ISOAREA_IMAGE']>30) & (self.table1['ISOAREA_IMAGE']<300)]
            self.table1['spatial resolution'] = float(0)
            self.table1['spectral resolution'] = float(0)
            self.table1['spatial resolution'] = self.table1['spatial resolution'].astype(float)
            self.table1['spectral resolution'] = self.table1['spectral resolution'].astype(float)
#	    print len(self.table), 'spots found at positions' , self.table['X_IMAGE'], self.table['Y_IMAGE']

    def extractBar(self, **kwargs):
        image = fits.open('stack310037.fits')[0].data
        n = 20
        n1 = 10
        box1 = [325,452]
        box2 = [755,566]
        im1 = image[box1[1]-n:box1[1]+n,box1[0]-n:box1[0]+n]
        xmax1 = np.sum(im1, axis=0).argmax()
        ymax1 = np.sum(im1, axis=1).argmax()
        im2 = image[box2[1]-n:box2[1]+n,box2[0]-n:box2[0]+n]
        xmax2 = np.sum(im2, axis=0).argmax()
        ymax2 = np.sum(im2, axis=1).argmax()
        im1 = im1[ymax1-n1:ymax1+n1,xmax1-n1:xmax1+n1]
        im2 = im2[ymax2-n1:ymax2+n1,xmax2-n1:xmax2+n1]
        imshow(im2)
        sum1 = np.sum(im1)
        sum2 = np.sum(im2)
        im1y = np.sum(im1,axis = 1)
        im1x = np.sum(im1,axis = 0)
        im2y = np.sum(im2,axis = 1)
        im2x = np.sum(im2,axis = 0)
        centroid1x = np.sum(im1x*np.arange(2*n1))/sum1
        centroid1y = np.sum(im1y*np.arange(2*n1))/sum1
        centroid2x = np.sum(im2x*np.arange(2*n1))/sum2
        centroid2y = np.sum(im2y*np.arange(2*n1))/sum2
        centroid1x_abs = box1[0] +xmax1  + centroid1x
        centroid1y_abs = box1[1] +ymax1  + centroid1y
        centroid2x_abs = box2[0] +xmax2  + centroid2x
        centroid2y_abs = box2[1] +ymax2  + centroid2y
        dist = np.sqrt(np.square(centroid1x_abs - centroid2x_abs) + np.square(centroid1y_abs - centroid2y_abs))
        return dist                         
        
       
    def extractSourcesMax(self, plot, **kwargs):
        """
        """ 

        n=2   
        data = self.image
        data2 = ndimage.grey_dilation(ndimage.grey_erosion(data, size=(n,n)), size=(n,n))       
        x,y = np.where(data2 == data2.max())[0][0] , np.where(data2 == data2.max())[1][0]
        daofind = DAOStarFinder(fwhm=5, threshold=5*std, ratio = 1)    
        sources0 = daofind(data2 - median)
        
    def extractSourcesPy(self, plot, **kwargs):
        """
        """ 
        quick = kwargs.get('quick', True)
        verbose = kwargs.get('verbose', False)
        n=2
        data = self.image
        data2 = ndimage.grey_dilation(ndimage.grey_erosion(data, size=(n,n)), size=(n,n))       
        self.mean, self.median, self.std = sigma_clipped_stats(data2, sigma=3.0, iters=5)    
        
        if quick:
            daofind = DAOStarFinder(fwhm=15, threshold=10*self.std, ratio = 1)         
            sources0 = daofind(data2 - self.median) 
        if not quick:
            daofind = DAOStarFinder(fwhm=15, threshold=10*self.std, ratio = 1)         
            sources0 = daofind(data2 - self.median)
            for i in [10,12.5,17.5,20]:
                data2 = ndimage.grey_dilation(ndimage.grey_erosion(data2, size=(n,n)), size=(n,n))
                mean, median, std = sigma_clipped_stats(data2, sigma=3.0, iters=5)    
                daofind = DAOStarFinder(fwhm=i, threshold=i*self.std, ratio = 1)    
                sources1 = daofind(data2 - median) 
                try :
                    sources0 = Table(np.hstack((sources0,sources1)))
                except TypeError:
                    print 'catalog empty'
        if verbose:
            print '[STATISTICS :] mean = {}, median = {}, std = {} before masking sources'.format(mean, median, std) 
        try:
            sources = delete_doublons(sources0, dist = 25)
        except TypeError:
            print 'No source detected on the image'
            return 0
#        e=0
#        for i in range(len(sources)):
#            if self.image[sources['xcentroid'][i-e],sources['xcentroid'][i-e]]<self.image.max()/2:
#                sources.remove_row(i-e) 
#                e+=1
        self.sources = delete_doublons(sources0, dist = 50)
        try:
            positions = (self.sources['xcentroid'], self.sources['ycentroid'])
            self.sources['X_IMAGE'] = self.sources['xcentroid']
            self.sources['Y_IMAGE'] = self.sources['ycentroid']
            self.sources['NUMBER'] = self.sources['id']        
            self.sources['FWHMx'] = float(0)
            self.sources['FWHMy'] = float(0)
            self.sources['spatial resolution'] = float(0)
            self.sources['spectral resolution'] = float(0)
            self.sources['spatial resolution'] = self.sources['spatial resolution'].astype(float)
            self.sources['spectral resolution'] = self.sources['spectral resolution'].astype(float)
            self.table = self.sources
            self.table1 = self.sources
        except TypeError:
            return 0

#        figure()
#        plt.imshow(data2)#, origin='lower', norm=norm)#, cmap='Greys', origin='lower', norm=norm)
#        apertures.plot(color='blue', lw=5.5, alpha=1.5)#sources1['X_IMAGE'] = sources1['xcentroid'] 
#        plt.title('Image - median after erosion dilatation' )
#        plt.show()

#        mask = make_source_mask(self.image, snr=1.5, npixels=5, dilate_size=50)#
##        figure();imshow(self.image - median);colorbar();plt.title('Image - median without sources' );show()
#        mean, median, std = sigma_clipped_stats(self.image, sigma=3.0, mask=mask)
#        print '[STATISTICS :] mean = {}, median = {}, std = {} before masking sources'.format(mean, median, std) 
#        self.image = self.image - median
##        figure();imshow(self.image - median);colorbar();show()
        
        
    def computeDispersionZn(self, **kwargs):
        """
        """ 
        rays = [2139,2062,2025]
        self.table.sort('X_IMAGE')
        disp = []
        disp.append( (rays[2] - rays[0]) / (self.table[2]['X_IMAGE'] - self.table[0]['X_IMAGE']))
        disp.append( (rays[2] - rays[1]) / (self.table[2]['X_IMAGE'] - self.table[1]['X_IMAGE']))
        disp.append( (rays[1] - rays[0]) / (self.table[1]['X_IMAGE'] - self.table[0]['X_IMAGE']))
        
        print '\033[31mdispersion = \033[31m', disp

    def computeSpatialResolution(self, **kwargs):
        """
        """ 
        self.resx = np.mean(self.table['FWHM_IMAGE'])
        
    def computeSpectralResolution(self, **kwargs):
        """
        """ 
        self.resl = 2040 / (np.mean(self.table['FWHM_IMAGE']) * self.disp)

    def Fit2D(self, plot, n=40, **kwargs):
        popt = []
        pcov = [] 
        for i in range(10):#len(self.table)):
#            print i
            image = self.image[int(self.table[i]['Y_IMAGE']-n):int(self.table[i]['Y_IMAGE']+n),int(self.table[i]['X_IMAGE']-n):int(self.table[i]['X_IMAGE']+n)]            
            image = ndimage.grey_dilation(ndimage.grey_erosion(image, size=(2,2)), size=(2,2))
            lx,ly = image.shape
            x = np.linspace(0,lx-1,lx)
            y = np.linspace(0,ly-1,ly)
            x, y = np.meshgrid(x,y)
            try:
                a = curve_fit(twoD_Gaussian,(x,y),image.flat,(np.mean(image),lx/2,ly/2,15,15,np.min(image)))
                popt.append(a[0])
#            popt.append(Normalization(*a[0]))
                pcov.append(a[1])
            except RuntimeError:
                print ('Optimal parameters not found... Trying next image')
            if plot:
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
                fig.suptitle('spatial resolution = {} \n spectral resolution = {}'.format(2.355*self.pix2sec*popt[i][3],2082/(2.355*popt[i][4]*self.disp)))
                plt.show()
        print '\033[31mfor the {} spots,\n the spatial resolutions are {},\n and the spectral resolutions are {}\033[31m'.format(len(self.table),2.355*0.95*np.array(popt)[:,3],2082/(2.355*np.array(popt)[:,4]*self.disp))
        return popt,pcov

    def Fit1D1(self, plot, stack, n1=3, n2=50,  **kwargs):#12
        popt = []
        pcov = [] 
        n=1
        self.table1['spectral resolution'] = 0
        e=0     
#        print (len(self.table1['spectral resolution']), 'images to analyze')
        for i in range(len(self.table1)):
            image = self.image[int(self.table1[i]['Y_IMAGE']-n1):int(self.table1[i]['Y_IMAGE']+n1),int(self.table1[i]['X_IMAGE']-n2):int(self.table1[i]['X_IMAGE']+n2)]
#            image = self.image
            image = np.sum(ndimage.grey_dilation(ndimage.grey_erosion(image, size=(n,n)), size=(n,n)), axis=0)
            lx = image.shape[0]
            x = np.linspace(0,lx-1,lx)
            try:
                poptww, pcovww = curve_fit(Gaussian,x ,image.flat,(np.mean(image),lx/2,1,np.min(image)))#, sigma = Gaussian(np.linspace(0,lx-1,lx),), absolute_sigma = True)
#                weight = Gaussian1(np.linspace(0,lx-1,lx),np.mean(image),lx/2,8)# / np.max(Gaussian1(np.linspace(0,lx-1,lx),np.mean(image),lx/2,1))
                weight = 1 / (0.2+Gaussian(np.linspace(0,lx-1,lx),*poptww) / np.sum(1/Gaussian(np.linspace(0,lx-1,lx),*poptww)))
                a = curve_fit(Gaussian,x ,image.flat,(np.max(image),lx/2,1,np.min(image)), sigma = weight, absolute_sigma = True)
#                a = curve_fit(Gaussian1,x ,weight*image.flat/np.max(weight),(np.max(image),lx/2,1))#, sigma = weight, absolute_sigma = True)
#                a = curve_fit(Gaussian,x ,weight*image.flat/np.max(weight),(np.max(image),lx/2,1,np.min(image)))#, sigma = weight, absolute_sigma = True)

    #                a = curve_fit(Gaussian,x ,image.flat,(np.mean(image),lx/2,1,np.min(image)), sigma = weight, absolute_sigma = True)
                popt.append(a[0])
                pcov.append(trace(a[1]))
                e+=1
                if plot:
                    plt.figure()
                    plt.plot(image,'x')
#                    plt.plot(Gaussian(x, *popt[i]), label = 'Weighted by first fit')
    #                plt.plot(x, weight, label = 'Weight')
                    plt.plot(Gaussian(x, *poptww), label = 'First fit without weight')
                    plt.title('{},    FWHMx = {} pix'.format(i, 2.355*poptww[2]))#,2082/(2.355*popt[i][3]*self.disp)))
                    plt.legend(loc = 'upper right')
                    plt.show()
                    var = raw_input("Do you want to keep this fit? : ")
                    if var=='y':
                        self.table1['FWHMx'][i] = 2.355*poptww[2]
                    if (var != 'y') & (var != 'n'):
                        print ("Not understood, wee kept this fit!")
                        self.table1['FWHMx'][i] = 2.355*poptww[2]                              
#                if (2082/(2.355*a[0][2]*self.disp))>500 or (2082/(2.355*a[0][2]*self.disp))<3000:
#                    self.table1['spectral resolution'][i] = 2082/(2.355*a[0][2]*self.disp)
                else:
                    self.table1['FWHMx'][i] = 2.355*poptww[2] 
            except RuntimeError:
                popt.append(0)
                print 'Image non analyzed: ', e                
                self.table1['FWHMx'][i] = 0                              
                print RuntimeError               
            except ValueError:
                popt.append(0)
                print ('Value error... Trying next image',i)                
                self.table1['FWHMx'][i] = 0                              
        return popt,pcov

    def Fit1D2(self, plot, stack, n1=50, n2=3,  **kwargs):#12
        popt = []
        pcov = [] 
        n=0
        self.table1['spectral resolution'] = 0
        e=0     
#        print (len(self.table1['spectral resolution']), 'images to analyze')
        for i in range(len(self.table1)):
#            print 'Image : ', i
            image = self.image[int(self.table1[i]['Y_IMAGE']-n1):int(self.table1[i]['Y_IMAGE']+n1),int(self.table1[i]['X_IMAGE']-n2):int(self.table1[i]['X_IMAGE']+n2)]
#            image = self.image
            image = np.sum(ndimage.grey_dilation(ndimage.grey_erosion(image, size=(n,n)), size=(n,n)), axis=1)
            lx = image.shape[0]
            x = np.linspace(0,lx-1,lx)
            try:
                poptww, pcovww = curve_fit(Gaussian,x ,image.flat,(np.max(image),lx/2,1,np.min(image)))#, sigma = Gaussian(np.linspace(0,lx-1,lx),), absolute_sigma = True)
#                weight = Gaussian1(np.linspace(0,lx-1,lx),np.mean(image),lx/2,8)# / np.max(Gaussian1(np.linspace(0,lx-1,lx),np.mean(image),lx/2,1))
                weight = 1 / (10+np.sqrt(Gaussian(np.linspace(0,lx-1,lx),*poptww) / np.sum(1/Gaussian(np.linspace(0,lx-1,lx),*poptww))))
#                plt.plot(weight)
                a = curve_fit(Gaussian,x ,image.flat,(np.max(image),lx/2,1,np.min(image)), sigma = weight, absolute_sigma = True)
#                a = curve_fit(Gaussian1,x ,weight*image.flat/np.max(weight),(np.max(image),lx/2,1))#, sigma = weight, absolute_sigma = True)
#                a = curve_fit(Gaussian,x ,weight*image.flat/np.max(weight),(np.max(image),lx/2,1,np.min(image)))#, sigma = weight, absolute_sigma = True)

    #                a = curve_fit(Gaussian,x ,image.flat,(np.mean(image),lx/2,1,np.min(image)), sigma = weight, absolute_sigma = True)
                popt.append(a[0])
                pcov.append(trace(a[1]))
                e+=1
                if plot:
                    plt.figure()
                    plt.plot(image,'x')
#                    plt.plot(Gaussian(x, *popt[i]), label = 'Weighted by first fit')
    #                plt.plot(x, weight, label = 'Weight')
                    plt.plot(Gaussian(x, *poptww), label = 'First fit without weight')
                    plt.title('{},    FWHM = {} pix'.format(i, 2.355*poptww[2]))#,2082/(2.355*popt[i][3]*self.disp)))
                    plt.legend(loc = 'upper right')
                    plt.show()
                    var = raw_input("Do you want to keep this fit? : ")
                    if var=='y':
                        self.table1['FWHMy'][i] = 2.355*poptww[2]
                    if (var != 'y') & (var != 'n'):
                        print ("Not understood, wee kept this fit!")
                        self.table1['FWHMy'][i] = 2.355*poptww[2]      
                else:
                    self.table1['FWHMy'][i] = 2.355*poptww[2]      
#                if (2082/(2.355*a[0][2]*self.disp))>500 or (2082/(2.355*a[0][2]*self.disp))<3000:
#                    self.table1['spectral resolution'][i] = 2082/(2.355*a[0][2]*self.disp)
            except RuntimeError:
                popt.append(0)
                print 'Image non analyzed: ', e                
                self.table1['FWHMy'][i] = 0                              
                print RuntimeError               
            except ValueError:
                popt.append(0)
                print ('Value error... Trying next image',i)                
                self.table1['FWHMy'][i] = 0                              
        return popt,pcov



    def createCatalog(self, number, **kwargs):
        fname = kwargs.get('filename', False)
        if 'FLUX_ISO' in self.table1.colnames:
            self.table1.write(self._DFLT_catalog + fname + self._DFLT_extension,overwrite = True)
        if 'xcentroid' in self.table1.colnames:        
            self.table1.write(self._DFLT_catalog + 'py{}'.format(number) + self._DFLT_extension,overwrite = True)
            

            
if __name__ == '__main__':
#    commandLine()
#    
#    number = 17
#    try:
#        fits.delval(image.format(number) + '.fits','NAXIS3')
#    except KeyError:
#        print 'ok'
#    for i in [55,60,65,70,75]:
#    number1 = []
#    number2 = []
#    number3 = []
#    xx = []
#    yy = []
#
#    filel = []
#    files = glob.glob('*.fits')[::-1]
#    images = []
#    for i in range(int(len(files))):
#        im = Image(filename = files[i], plot = False, stack = True, stack_image = False, py = True, subDark = False, verbose = False, quick = True, Type = 'guider')
##        im = Image(filename = 'stack557607.fits', plot = False, stack = True, stack_image = False, py = True, subDark = False, verbose = False, quick = False, Type = 'guider')
#        images.append(im)
#        if im.abort == 0 :
#            print 'no sources detected'
#        else:
#            table = im.table[(im.table['FWHMx']>1) &(im.table['FWHMx']<25) & (im.table['FWHMy']>1)&(im.table['FWHMx']<25)]
#            print files[i]
#            plot_FWHM(table,im.image,name=files[i][:-5])


#        print im.table['FWHMx'],im.table['FWHMy']
#        table = im.table[(im.table['FWHMx']>1) &(im.table['FWHMx']<20) & (im.table['FWHMy']>1)&(im.table['FWHMx']<20)]
#        f, ax = plt.subplots(1,2, sharey=True)
#        table.sort('xcentroid')
#        ax[0].plot(table['xcentroid'],table['FWHMx'],'--o')
#        ax[0].plot(table['xcentroid'],table['FWHMy'],'--o')
#        ax[0].grid()
#        ax[0].set_title('FWHM size in x dir')
#        table.sort('ycentroid')
#        ax[1].plot(table['ycentroid'],table['FWHMx'],'--o')
#        ax[1].plot(table['ycentroid'],table['FWHMy'],'--o')
#        ax[1].grid()
#        ax[1].set_ylim((4,20))
#        ax[1].set_title('FWHM size in y dir')
#        plt.show()



#        var = raw_input("next? :")
#        if var != 'n':
#            print '\n \n \n \n \n'
        
#        im.table.sort('xcentroid')
#        plt.figure();plt.plot(im.table['xcentroid'],im.table['FWHMx'],'--o'),plt.plot(im.table['xcentroid'],im.table['FWHMy'],'--o');plt.grid();plt.show
#        im.table.sort('ycentroid')
#        plt.figure();plt.plot(im.table['ycentroid'],im.table['FWHMx'],'--o'),plt.plot(im.table['ycentroid'],im.table['FWHMy'],'--o');plt.grid();plt.show
        
    

#il faut que le code soit generique:
#--ouvre n'importe quel format que ce soit fits ou autre
#-- verbose ou non verbose
#-- plot ou pas plot (deja pas mal) avec demande davis
#--vite ou pas vite avec la boucle sur la recherche 
#--guider ccd ou autre
# dark ou pas dark
# une ou plsusieurs sources

    def extractStarBPondere(file):
        image = fits.open(file)[0].data
        n = 30
        n1 = 8
        box1 = [325,452]
        box2 = [755,566]
        im1 = image[box1[1]-n:box1[1]+n,box1[0]-n:box1[0]+n]
        xmax1 = np.sum(im1, axis=0).argmax()
        ymax1 = np.sum(im1, axis=1).argmax()
        im2 = image[box2[1]-n:box2[1]+n,box2[0]-n:box2[0]+n]
        xmax2 = np.sum(im2, axis=0).argmax()
        ymax2 = np.sum(im2, axis=1).argmax()
        im1 = im1[ymax1-n1:ymax1+n1,xmax1-n1:xmax1+n1]
        im2 = im2[ymax2-n1:ymax2+n1,xmax2-n1:xmax2+n1]
        sum1 = np.sum(im1)
        sum2 = np.sum(im2)
        try:
            im1y = np.sum(im1,axis = 1)
            im1x = np.sum(im1,axis = 0)
            im2y = np.sum(im2,axis = 1)
            im2x = np.sum(im2,axis = 0)
            centroid1x = np.sum(im1x*np.arange(2*n1))/sum1
            centroid1y = np.sum(im1y*np.arange(2*n1))/sum1
            centroid2x = np.sum(im2x*np.arange(2*n1))/sum2
            centroid2y = np.sum(im2y*np.arange(2*n1))/sum2
            centroid1x_abs = box1[0] +xmax1  + centroid1x
            centroid1y_abs = box1[1] +ymax1  + centroid1y
            centroid2x_abs = box2[0] +xmax2  + centroid2x
            centroid2y_abs = box2[1] +ymax2  + centroid2y
            dist = np.sqrt(np.square(centroid1x_abs - centroid2x_abs) + np.square(centroid1y_abs - centroid2y_abs))
            return dist 
        except ValueError:
            print 'Center is not at the center'
            return nan


    
    def extractStarB(file):
        image = fits.open(file)[0].data
        n = 30
        n1 = 8
        box1 = [325,452]
        box2 = [755,566]
        im1 = image[box1[1]-n:box1[1]+n,box1[0]-n:box1[0]+n]
        xmax1 = np.sum(im1, axis=0).argmax()
        ymax1 = np.sum(im1, axis=1).argmax()
        im2 = image[box2[1]-n:box2[1]+n,box2[0]-n:box2[0]+n]
        xmax2 = np.sum(im2, axis=0).argmax()
        ymax2 = np.sum(im2, axis=1).argmax()
        im1 = im1[ymax1-n1:ymax1+n1,xmax1-n1:xmax1+n1]
        im2 = im2[ymax2-n1:ymax2+n1,xmax2-n1:xmax2+n1]
        sum1 = np.sum(im1)
        sum2 = np.sum(im2)
        try:
            im1y = np.sum(im1,axis = 1)
            im1x = np.sum(im1,axis = 0)
            im2y = np.sum(im2,axis = 1)
            im2x = np.sum(im2,axis = 0)
            centroid1x = np.sum(im1x*np.arange(2*n1))/sum1
            centroid1y = np.sum(im1y*np.arange(2*n1))/sum1
            centroid2x = np.sum(im2x*np.arange(2*n1))/sum2
            centroid2y = np.sum(im2y*np.arange(2*n1))/sum2
            centroid1x_abs = box1[0] +xmax1  + centroid1x
            centroid1y_abs = box1[1] +ymax1  + centroid1y
            centroid2x_abs = box2[0] +xmax2  + centroid2x
            centroid2y_abs = box2[1] +ymax2  + centroid2y
            dist = np.sqrt(np.square(centroid1x_abs - centroid2x_abs) + np.square(centroid1y_abs - centroid2y_abs))
            return dist 
        except ValueError:
            print 'Center is not at the center'
            return nan
    
    def GaussianM(x, x0,amp):
        xo = np.array(x0)
        amplitude = amp
        sigma = 2.5
        A = amplitude/(sigma * np.sqrt(2*math.pi))    
        g = 8600 + amplitude * np.exp( - 0.5*(((x-xo)/sigma)**2))
        return g

    def GaussianM1(x, x0):
        xo = float(x0)
        amplitude = 1000
        sigma = 2.5
        A = amplitude/(sigma * np.sqrt(2*math.pi))    
        g = 8600 + amplitude * np.exp( - 0.5*(((x-xo)/sigma)**2))
        return g.ravel()
        
    def extractStarF(file, plot = False):
        image = fits.open(file)[0].data
        n = 30
        n1 = 12
    
        box1 = [320,457]
        box2 = [755,566]
        box3 = [521,494]
        im1 = image[box1[1]-n:box1[1]+n,box1[0]-n:box1[0]+n]
        im2 = image[box2[1]-n:box2[1]+n,box2[0]-n:box2[0]+n]
        im3 = image[box3[1]-n:box3[1]+n,box3[0]-n:box3[0]+n]
        xmax1 = np.sum(im1, axis=0).argmax()
        ymax1 = np.sum(im1, axis=1).argmax()
        xmax2 = np.sum(im2, axis=0).argmax()
        ymax2 = np.sum(im2, axis=1).argmax()
        xmax3 = np.sum(im3, axis=0).argmax()
        ymax3 = np.sum(im3, axis=1).argmax()
        im1 = im1[ymax1-n1:ymax1+n1,xmax1-n1:xmax1+n1]
        im2 = im2[ymax2-n1:ymax2+n1,xmax2-n1:xmax2+n1]
        im3 = im3[ymax3-n1:ymax3+n1,xmax3-n1:xmax3+n1]
        sum1 = np.sum(im1)
        sum2 = np.sum(im2)
        sum3 = np.sum(im3)
        x = np.arange(2*n1)
        im1y = np.sum(im1,axis = 1)
        im1x = np.sum(im1,axis = 0)
        im2y = np.sum(im2,axis = 1)
        im2x = np.sum(im2,axis = 0)
        try:
            centroid1x, cov1x = curve_fit(Gaussian,x ,im1x,(np.max(im1x),n1,1,np.min(im1x)))
            centroid1y, cov1y = curve_fit(Gaussian,x ,im1y,(np.max(im1y),n1,1,np.min(im1y)))
            centroid2x, cov2x = curve_fit(Gaussian,x ,im2x,(np.max(im2x),n1,1,np.min(im2x)))
            centroid2y, cov2y = curve_fit(Gaussian,x ,im2y,(np.max(im2y),n1,1,np.min(im2y)))
            centroid2xM, cov2xM = curve_fit(GaussianM1,x ,im2x,n1)
            centroid2yM, cov2yM = curve_fit(GaussianM1,x ,im2y,n1)
#            centroid1x, cov1x = curve_fit(Gaussian,x ,im1x,(np.max(im1x),n1,1,np.min(im1x)))
#            centroid1y, cov1y = curve_fit(Gaussian,x ,im1y,(np.max(im1y),n1,1,np.min(im1y)))
#            centroid2xa.append(curve_fit(Gaussian,x ,im2x,(np.max(im2x),n1,1,np.min(im2x)))[0])
#            centroid2ya.append(curve_fit(Gaussian,x ,im2y,(np.max(im2y),n1,1,np.min(im2y)))[0])
#            print centroid2x,centroid2y
#            plt.figure();plt.plot(x,Gaussian(x,*centroid1x),label = 'cov = {}'.format(np.log10(cov1x.trace())));plt.plot(x,Gaussian(x,*centroid1y),label = 'cov = {}'.format(np.log10(cov1y.trace())));plt.plot(x,Gaussian(x,*centroid2x),label = 'cov = {}'.format(np.log10(cov2x.trace())));plt.plot(x,Gaussian(x,*centroid2y),label = 'cov = {}'.format(np.log10(cov2y.trace())));plt.legend(loc='upper right')
            if plot:
                plt.figure();plt.title(file)
                plt.plot(x,Gaussian(x,*centroid1x),label = 'cov = {}'.format(cov1x[1,1]));plt.plot(x,Gaussian(x,*centroid1y),label = 'cov = {}'.format(cov1y[1,1]));plt.plot(x,Gaussian(x,*centroid2x),label = 'cov = {}'.format(cov2x[1,1]));plt.plot(x,Gaussian(x,*centroid2y),label = 'cov = {}'.format(cov2y[1,1]));plt.legend(loc='upper right')
                plt.plot(im1x,'x');plt.plot(im1y,'x');plt.plot(im2x,'x');plt.plot(im2y,'x')
                plt.plot(x,GaussianM1(x,*centroid2xM));plt.plot(x,GaussianM1(x,*centroid2yM));plt.show()
            centroid1x_abs = box1[0] +xmax1  + centroid1x[1]
            centroid1y_abs = box1[1] +ymax1  + centroid1y[1]
            centroid2x_abs = box2[0] +xmax2  + centroid2x[1]
            centroid2y_abs = box2[1] +ymax2  + centroid2y[1]
            centroid2x_abs = box2[0] +xmax2  + centroid2x[1]
            centroid2y_abs = box2[1] +ymax2  + centroid2y[1]
            error=0.1
            if (cov1x[1,1]<error) & (cov1y[1,1]<error) & (cov2x[1,1]<error) & (cov2y[1,1]<error):
                dist = np.sqrt(np.square(centroid1x_abs - centroid2x_abs) + np.square(centroid1y_abs - centroid2y_abs))
            else:
                dist = nan
            print 'ok'
            return dist 
        except ValueError:
            print 'Center is not at the center'
            return nan
        except RuntimeError:
            print ('did not converge')
#            plt.plot(im1x,'x');plt.plot(im1y,'x');plt.plot(im2x,'x');plt.plot(im2y,'x');plt.show()
            return nan



    def extractStarFW(file):
        image = fits.open(file)[0].data
        n = 24
        n1 = 12    
        box1 = [320,457]
        box2 = [755,566]
        box3 = [521,494]
        im1 = image[box1[1]-n:box1[1]+n,box1[0]-n:box1[0]+n]
        im2 = image[box2[1]-n:box2[1]+n,box2[0]-n:box2[0]+n]
        im3 = image[box3[1]-n:box3[1]+n,box3[0]-n:box3[0]+n]
        xmax1 = np.sum(im1, axis=0).argmax()
        ymax1 = np.sum(im1, axis=1).argmax()
        xmax2 = np.sum(im2, axis=0).argmax()
        ymax2 = np.sum(im2, axis=1).argmax()
        xmax3 = np.sum(im3, axis=0).argmax()
        ymax3 = np.sum(im3, axis=1).argmax()
        im1 = im1[ymax1-n1:ymax1+n1,xmax1-n1:xmax1+n1]
        fond2y = np.sum(im2[ymax2-n1:ymax2+n1,:xmax2-n1], axis=1) + np.sum(im2[ymax2-n1:ymax2+n1,xmax2+n1:], axis=1)
        fond2x = np.sum(im2[ymax2-n1:,xmax2-n1:xmax2+n1], axis=0) + np.sum(im2[ymax2-n1:ymax2+n1,xmax2-n1:xmax2+n1], axis=0)
        im2 = im2[ymax2-n1:ymax2+n1,xmax2-n1:xmax2+n1]
        im3 = im3[ymax3-n1:ymax3+n1,xmax3-n1:xmax3+n1]
        x = np.arange(2*n1)
        try:
            im1y = np.sum(im1,axis = 1)
            im1x = np.sum(im1,axis = 0)
            im2y = np.sum(im2,axis = 1) #- fond2
            im2x = np.sum(im2,axis = 0)
            weight = Gaussian(np.arange(24),-100,n1,4,10)
#            plot(weight)
            centroid1x = curve_fit(Gaussian,x ,im1x,(np.max(im1x),n1,1,np.min(im1x)), sigma = weight, absolute_sigma = True)[0][1]
            centroid1y = curve_fit(Gaussian,x ,im1y,(np.max(im1y),n1,1,np.min(im1y)), sigma = weight, absolute_sigma = True)[0][1]
            centroid2x = curve_fit(Gaussian,x ,im2x,(np.max(im2x),n1,1,np.min(im2x)), sigma = weight, absolute_sigma = True)[0][1]
            centroid2y = curve_fit(Gaussian,x ,im2y,(np.max(im2y),n1,1,np.min(im2y)), sigma = weight, absolute_sigma = True)[0][1]
            centroid1x_abs = box1[0] +xmax1  + centroid1x
            centroid1y_abs = box1[1] +ymax1  + centroid1y
            centroid2x_abs = box2[0] +xmax2  + centroid2x
            centroid2y_abs = box2[1] +ymax2  + centroid2y
            dist = np.sqrt(np.square(centroid1x_abs - centroid2x_abs) + np.square(centroid1y_abs - centroid2y_abs))
            print 'ok'
            return dist 
        except ValueError:
            print 'Center is not at the center'
            return nan
        except RuntimeError:
            print 'Fit did not converge'
            return nan

    def extractStarF2(file, plot = False):
        image = fits.open(file)[0].data
        n = 30
        n1 = 12
    
        box1 = [320,457]
        box2 = [755,566]
        box3 = [521,494]
        im1 = image[box1[1]-n:box1[1]+n,box1[0]-n:box1[0]+n]
        im2 = image[box2[1]-n:box2[1]+n,box2[0]-n:box2[0]+n]
        im3 = image[box3[1]-n:box3[1]+n,box3[0]-n:box3[0]+n]
        xmax1 = np.sum(im1, axis=0).argmax()
        ymax1 = np.sum(im1, axis=1).argmax()
        xmax2 = np.sum(im2, axis=0).argmax()
        ymax2 = np.sum(im2, axis=1).argmax()
        xmax3 = np.sum(im3, axis=0).argmax()
        ymax3 = np.sum(im3, axis=1).argmax()
        im1 = im1[ymax1-n1:ymax1+n1,xmax1-n1:xmax1+n1]
        im2 = im2[ymax2-n1:ymax2+n1,xmax2-n1:xmax2+n1]
        im3 = im3[ymax3-n1:ymax3+n1,xmax3-n1:xmax3+n1]
        sum1 = np.sum(im1)
        sum2 = np.sum(im2)
        sum3 = np.sum(im3)
        x = np.arange(2*n1)
        im1y = np.sum(im1,axis = 1)
        im1x = np.sum(im1,axis = 0)
        im2y = np.sum(im2,axis = 1)
        im2x = np.sum(im2,axis = 0)
        try:
            centroid1x, cov1x = curve_fit(Gaussian,x ,im1x,(np.max(im1x),n1,1,np.min(im1x)))
            centroid1y, cov1y = curve_fit(Gaussian,x ,im1y,(np.max(im1y),n1,1,np.min(im1y)))
            centroid2x, cov2x = curve_fit(Gaussian,x ,im2x,(np.max(im2x),n1,1,np.min(im2x)))
            centroid2y, cov2y = curve_fit(Gaussian,x ,im2y,(np.max(im2y),n1,1,np.min(im2y)))
            centroid2xM, cov2xM = curve_fit(GaussianM,x ,im2x,(n1,max(im2x)))
            centroid2yM, cov2yM = curve_fit(GaussianM,x ,im2y,(n1,max(im2x)))
#            centroid1x, cov1x = curve_fit(Gaussian,x ,im1x,(np.max(im1x),n1,1,np.min(im1x)))
#            centroid1y, cov1y = curve_fit(Gaussian,x ,im1y,(np.max(im1y),n1,1,np.min(im1y)))
#            centroid2xa.append(curve_fit(Gaussian,x ,im2x,(np.max(im2x),n1,1,np.min(im2x)))[0])
#            centroid2ya.append(curve_fit(Gaussian,x ,im2y,(np.max(im2y),n1,1,np.min(im2y)))[0])
#            print centroid2x,centroid2y
#            plt.figure();plt.plot(x,Gaussian(x,*centroid1x),label = 'cov = {}'.format(np.log10(cov1x.trace())));plt.plot(x,Gaussian(x,*centroid1y),label = 'cov = {}'.format(np.log10(cov1y.trace())));plt.plot(x,Gaussian(x,*centroid2x),label = 'cov = {}'.format(np.log10(cov2x.trace())));plt.plot(x,Gaussian(x,*centroid2y),label = 'cov = {}'.format(np.log10(cov2y.trace())));plt.legend(loc='upper right')
            if plot:
                plt.figure();plt.title(file)
                plt.plot(x,Gaussian(x,*centroid1x),label = 'cov = {}'.format(cov1x[1,1]));plt.plot(x,Gaussian(x,*centroid1y),label = 'cov = {}'.format(cov1y[1,1]))
                #plt.plot(x,Gaussian(x,*centroid2x),label = 'cov = {}'.format(cov2xM[0,0]));plt.plot(x,Gaussian(x,*centroid2y),label = 'cov = {}'.format(cov2yM[0,0]))
                plt.legend(loc='upper right')
                plt.plot(im1x,'x');plt.plot(im1y,'x');plt.plot(im2x,'x');plt.plot(im2y,'x')
                plt.plot(x,GaussianM(x,*centroid2xM));plt.plot(x,GaussianM(x,*centroid2yM));plt.show()
            centroid1x_abs = box1[0] +xmax1  + centroid1x[1]
            centroid1y_abs = box1[1] +ymax1  + centroid1y[1]
            centroid2x_abs = box2[0] +xmax2  + centroid2x[1]
            centroid2y_abs = box2[1] +ymax2  + centroid2y[1]
            centroid2x_abs = box2[0] +xmax2  + centroid2xM[0]
            centroid2y_abs = box2[1] +ymax2  + centroid2yM[0]
            error=0.4
#            if (cov1x[1,1]<error) & (cov1y[1,1]<error) & (cov2x[1,1]<error) & (cov2y[1,1]<error):
            if (cov2xM[0,0]<error) & (cov2yM[0,0]<error):
                dist = np.sqrt(np.square(centroid1x_abs - centroid2x_abs) + np.square(centroid1y_abs - centroid2y_abs))
            else:
                dist = nan
            print 'ok'
            return dist 
        except ValueError:
            print 'Center is not at the center'
            return nan
        except RuntimeError:
            print ('did not converge')
#            plt.plot(im1x,'x');plt.plot(im1y,'x');plt.plot(im2x,'x');plt.plot(im2y,'x');plt.show()
            return nan

    def extractStarCC(file, plot = False):
        image = fits.open(file)[0].data
        n = 30
        n1 = 12
    
        box1 = [320,457]
        box2 = [755,566]
        box3 = [521,494]
        im1 = image[box1[1]-n:box1[1]+n,box1[0]-n:box1[0]+n]
        im2 = image[box2[1]-n:box2[1]+n,box2[0]-n:box2[0]+n]
        im3 = image[box3[1]-n:box3[1]+n,box3[0]-n:box3[0]+n]
        xmax1 = np.sum(im1, axis=0).argmax()
        ymax1 = np.sum(im1, axis=1).argmax()
        xmax2 = np.sum(im2, axis=0).argmax()
        ymax2 = np.sum(im2, axis=1).argmax()
        xmax3 = np.sum(im3, axis=0).argmax()
        ymax3 = np.sum(im3, axis=1).argmax()
        im1 = im1[ymax1-n1:ymax1+n1,xmax1-n1:xmax1+n1]
        im2 = im2[ymax2-n1:ymax2+n1,xmax2-n1:xmax2+n1]
        im3 = im3[ymax3-n1:ymax3+n1,xmax3-n1:xmax3+n1]
        sum1 = np.sum(im1)
        sum2 = np.sum(im2)
        sum3 = np.sum(im3)
        x = np.arange(2*n1)
        im1y = np.sum(im1,axis = 1)
        im1x = np.sum(im1,axis = 0)
        im2y = np.sum(im2,axis = 1)
        im2x = np.sum(im2,axis = 0)
        centre = np.linspace(n1/4,2*n1 - n1/4,10*n1)
        gaussians = GaussianM(x[np.newaxis,:],centre[:,np.newaxis],1)
        try:
            crosscorx1 = np.sum(gaussians * im1x,axis=1)
            crosscory1 = np.sum(gaussians * im1y,axis=1)            
            centroid1xM = centre[crosscorx1.argmax() ]
            centroid1yM = centre[crosscory1.argmax()]
            crosscorx2 = np.sum(gaussians * im2x,axis=1)
            crosscory2 = np.sum(gaussians * im2y,axis=1)  
            cxmax2 = crosscorx2.argmax()
            cymax2 = crosscory2.argmax()            
            cxmin2 = crosscorx2.argmin()
            cymin2 = crosscory2.argmin()
            centroid2xM = centre[cxmax2]
            centroid2yM = centre[cymax2]
#            centroid1x, cov1x = curve_fit(Gaussian,x ,im1x,(np.max(im1x),n1,1,np.min(im1x)))
#            centroid1y, cov1y = curve_fit(Gaussian,x ,im1y,(np.max(im1y),n1,1,np.min(im1y)))
#            centroid2xa.append(curve_fit(Gaussian,x ,im2x,(np.max(im2x),n1,1,np.min(im2x)))[0])
#            centroid2ya.append(curve_fit(Gaussian,x ,im2y,(np.max(im2y),n1,1,np.min(im2y)))[0])
#            print centroid2x,centroid2y
#            plt.figure();plt.plot(x,Gaussian(x,*centroid1x),label = 'cov = {}'.format(np.log10(cov1x.trace())));plt.plot(x,Gaussian(x,*centroid1y),label = 'cov = {}'.format(np.log10(cov1y.trace())));plt.plot(x,Gaussian(x,*centroid2x),label = 'cov = {}'.format(np.log10(cov2x.trace())));plt.plot(x,Gaussian(x,*centroid2y),label = 'cov = {}'.format(np.log10(cov2y.trace())));plt.legend(loc='upper right')
            
            errorx = crosscorx2.max() - crosscorx2.min()
            errory = crosscory2.max() - crosscory2.min()
#            print errory
            if plot:
                plt.figure();plt.title(file)
                plt.plot(x,Gaussian(x,*centroid1x),label = 'cov = {}'.format(errorx));plt.plot(x,Gaussian(x,*centroid1y),label = 'cov = {}'.format(errory))
                #plt.plot(x,Gaussian(x,*centroid2x),label = 'cov = {}'.format(cov2xM[0,0]));plt.plot(x,Gaussian(x,*centroid2y),label = 'cov = {}'.format(cov2yM[0,0]))
                plt.legend(loc='upper right')
                plt.plot(im1x,'x');plt.plot(im1y,'x');plt.plot(im2x,'x',label = 'cov = {}'.format(errorx));plt.plot(im2y,'x',label = 'cov = {}'.format(errory))
                plt.show()
            centroid1x_abs = box1[0] +xmax1  + centroid1xM
            centroid1y_abs = box1[1] +ymax1  + centroid1yM
            centroid2x_abs = box2[0] +xmax2  + centroid2xM
            centroid2y_abs = box2[1] +ymax2  + centroid2yM
            error=9000
            if (errorx>error) &(errory>error):
                dist = np.sqrt(np.square(centroid1x_abs - centroid2x_abs) + np.square(centroid1y_abs - centroid2y_abs))
            else:
                return nan, 0 
            return dist , errorx
        except ValueError:
            print 'Center is not at the center'
            return nan, 0
        except RuntimeError:
            print ('did not converge')
#            plt.plot(im1x,'x');plt.plot(im1y,'x');plt.plot(im2x,'x');plt.plot(im2y,'x');plt.show()
            return nan  , 0  
    c=[]
    d=[]
    e=[]
    f = []
    g = []
    files = glob.glob('stack*7.fits')
    for i in files:
        c.append(extractStarB(i))
#        d.append(extractStarF(i,plot=False))
        e.append(extractStarF2(i,plot=False))
        f.append(extractStarCC(i,plot=False)[0])
        g.append(extractStarCC(i,plot=False)[1])
    c= np.array(c - np.nanmean(c))    
    d= np.array(d - np.nanmean(d)) 
    e= np.array(e - np.nanmean(e)) 
    f= np.array(f - np.nanmean(f)) 
#    d = d[(d<1.5) & (d>-1.5)]
#    e = e[(e<1.5) & (e>-1.5)]
#    plt.plot(c,'o', label = 'barycentre, std = {}'.format(np.nanstd(c)))
#    plt.plot(d,'o', label = 'fit gaussien, std = {}'.format(np.nanstd(d)))
    plt.plot(f,'x', label = 'cross cor, std = {}'.format(np.nanstd(f)))
    plt.plot(e,'p', label = 'fit gaussien weighted, std = {}'.format(np.nanstd(e)))
#    plt.plot(np.array(g)/1000,'p', label = 'error of cross cor, std = {}'.format(np.nanstd(e)))
    plt.grid()
    plt.legend()
    plt.xlabel('number image')
    plt.ylabel('distance Stars (pix)')
    print np.nanstd(d)
    print np.nanmean(c)

#
#mean = np.mean(image[:200,1000:])
#var = np.var(image[:200,1000:])
#n = (mean/var)**2
#k = mean/n  
#
#
#snr = (image - mean)/np.sqrt((image + var ))
#
#flux - fond sur meme surface / racine(flux+var)
