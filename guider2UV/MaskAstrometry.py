#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 10:20:33 2017

@author: dvibert
"""
from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy.table import Table
from astropy import units as u

class LocalGuiderMaskProjector(object):
    
    def __init__(self, wcs):
        self.w = wcs
        self._build_local_frame()
        
    def __str__(self):
        str_ = '''
GuiderMaskProjector object:
    LocalFrame: {}
    wcs: {}
'''.format(self.localframe, self.w)
        return str_
        
    def _build_local_frame(self):
        center = self.w.wcs.crval*u.deg
        center = coordinates.SkyCoord(center[0], center[1])
        self.localframe = coordinates.SkyOffsetFrame(origin=center)
                                                    
    def pix2local(self, pix_coord, distortion=True):
        if distortion:
            radec  = self.w.all_pix2world(pix_coord, 1)
        else:
            radec  = self.w.wcs_pix2world(pix_coord, 1)
        radec_coord = coordinates.SkyCoord(radec, unit='deg')
        return radec_coord.transform_to(self.localframe)
    
    def local2pix(self, local_coord, distortion=True):
        radec_coord = local_coord.transform_to(coordinates.ICRS)
        if distortion:
            pix = self.w.all_world2pix(radec_coord.ra, radec_coord.dec, 1)
        else:
            pix = self.w.wcs_world2pix(radec_coord.ra, radec_coord.dec, 1)            
        return pix




class LocalScienceMaskProjector(object):
    
    def __init__(self, center, rotation=0*u.deg, gamma=1.):
        self.center = center
        self.rotation = rotation
        self.gamma = gamma
        self._build_local_frame()

    def __str__(self):
        str_ = '''
ScienceMaskProjector object:
    LocalFrame: {}
    gamma: {}
'''.format(self.localframe, self.gamma)
        return str_

    @property
    def radial_mag_polynomial(self):
        try:
            return self._radial_mag_polynomial
        except AttributeError:
            self._radial_mag_polynomial = (42.26134, -3.154411e-3, -1.117322) # 2018 masks
            return self._radial_mag_polynomial
            
    @property
    def radial_mag_inv_polynomial(self):
        try:
            return self._radial_mag_inv_polynomial
        except AttributeError:
            # self._radial_mag_inv_polynomial np.roots([c,b,a,-r])[-1]/r # last root seems the right one....
            self._radial_mag_inv_polynomial = (2.366233E-2, -3.610313e-9, 3.566165e-7) # 2018 masks
            return self._radial_mag_polynomial
    
    @radial_mag_polynomial.setter
    def radial_mag_polynomial(self, coeffs):
        self._radial_mag_polynomial = coeffs
    
    @radial_mag_inv_polynomial.setter
    def radial_mag_inv_polynomial(self, coeffs):
        self.radial_mag_inv_polynomial = coeffs

    # THIS SHOULD BE UPDATED for 2022 flight... 
    #  at least replace with a constant near the expected magnification.
    # (the model update method in guider2UV is only able to add a small 1+delta factor to this)/
    # alternatively: return the constant 1.0 here and set self.gamma to the expected value (in deg/mm eg 2.36E-2)  
    # when instanciating the class (eg when instanciating a Guider2UV object).
    def radial_magnification(self, r, reverse=False):
        if not reverse:
            a, b, c = self.radial_mag_polynomial
            return a + b * r + c * np.square(r) # mm/deg
        else:
            a, b, c = self.radial_mag_inv_polynomial
            return c*np.square(r) + b*r + a # deg/mm
        
        
    def  _build_local_frame(self):
        # revert lon axis, swap lat & lon (=> rotation -90deg)
        # to have same axes than guider ones
        self.localframe = coordinates.SkyOffsetFrame(origin=self.center, 
                                                     rotation=self.rotation-90*u.deg)

        
    def pix2local(self, pix_coord, distortion=True): #pix_coord in mm 
        r = np.sqrt(pix_coord[0]**2 + pix_coord[1]**2)
        inv_magnification = self.radial_magnification(r, reverse=True)
         
        pix_coord = pix_coord * inv_magnification * self.gamma
        
        # -90 deg rotatation to match guider x,y axis
        mrot = np.array([[0, 1],
                         [-1, 0]])
        local_coord = mrot.dot(pix_coord)
        local_coord = coordinates.SkyCoord(local_coord[0]*u.deg, local_coord[1]*u.deg, 
                                           frame=self.localframe)
        return local_coord
        
    
    def local2pix(self, local_coord, distortion=True):
        r = np.sqrt(local_coord.lon.deg**2 + local_coord.lat.deg**2)
        magnification = self.radial_magnification(r)

        xy = np.vstack((local_coord.lon.deg, local_coord.lat.deg))
        xy *= magnification / self.gamma
        
        # 90 deg rotation tp match x,y science mask axis
        mrot = np.array([[0, -1],
                         [1, 0]])    
        xy = mrot.dot(xy)        
        return xy


    def world2local(self, world_coord):
        return world_coord.transform_to(self.localframe)

    
    def local2world(self, local_coord):
        return local_coord.transform_to('icrs')
    

    def world2pix(self, world_coord, distortion=True):
        local_coord = self.world2local(world_coord)
        return self.local2pix(local_coord, distortion=distortion)

    def pix2world(self, pix_coord, distortion=True):
        local_coord = self.pix2local(pix_coord, distortion=True)
        return self.local2world(local_coord)
    
    
    
# test the LocalMaskProjector
if __name__ == '__main__':
    
    #from focustest2 import Image
    import matplotlib.pyplot as plt

    f3_astrometry_filename = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/perf3.0-09-01/sky_00_astrometry/stack475788.new"
    f3_astrometry_filename = "/Users/Vincent/Downloads/wcs-4.fits"

    #m13tab = Table.read('stack475788.starpos.fits')

    hdu = fits.open(f3_astrometry_filename)
    w  = wcs.WCS(hdu[0].header)

    F3_converter = LocalGuiderMaskProjector(w)
    


    local = F3_converter.pix2local(np.array([m13tab['xcentroid'],m13tab['ycentroid']]).T)
    plt.figure()
    plt.plot(local.lon, local.lat, 'o')

    pix = F3_converter.local2pix(local)
    plt.figure()
    plt.plot(m13tab['xcentroid'], m13tab['ycentroid'],'+r')
    plt.plot(pix[0], pix[1], 'x')
    
    local = F3_converter.pix2local(np.array([m13tab['xcentroid'],m13tab['ycentroid']]).T, distortion=False)
    plt.figure()
    plt.plot(local.lon, local.lat, 'og')

    pix = F3_converter.local2pix(local, distortion=False)
    plt.figure()
    plt.plot(m13tab['xcentroid'], m13tab['ycentroid'],'+r')
    plt.plot(pix[0], pix[1], 'x')
