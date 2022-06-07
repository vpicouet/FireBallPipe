#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 10:20:08 2018

@author: dvibert
"""

from __future__ import division, print_function

import os, sys
import numpy as np

from .polyfit import PolyFit

try:
   import cPickle as pickle # for python 2 only CPickle is faster
except:
   import pickle


def load_pickle(pickle_file):
    '''
    interface to pickle.load for mitigating python 2 & 3 incompatbilities
    '''

    if __package__ != '':
       package_pref = __package__ + '.'
    else:
       package_pref = __package__
       
    if __package__ != 'Calibration': # needed for DSutils: package name is changed to DS9FireBall
       sys.modules['Calibration.polyfit'] = sys.modules[package_pref + 'polyfit']
       sys.modules['Calibration.polynomial'] = sys.modules[package_pref + 'polynomial']

    if __package__ != '':  # needed for pkls generated from FireBallIMO code 
       sys.modules['polyfit'] = sys.modules[package_pref + 'polyfit']
       sys.modules['polynomial'] = sys.modules[package_pref + 'polynomial']
    
    try:
        with open(pickle_file, 'rb') as f:
            pickle_data = pickle.load(f)
    except UnicodeDecodeError as e:
        with open(pickle_file, 'rb') as f:
            pickle_data = pickle.load(f, encoding='latin1')
    except Exception as e:
        print('Unable to load data ', pickle_file, ':', e)
        raise

    if __package__ != 'Calibration':
       del sys.modules['Calibration.polyfit'], sys.modules['Calibration.polynomial']

    if __package__ != '':
       del sys.modules['polyfit'], sys.modules['polynomial']
    
    return pickle_data


class Mapping(object):
    '''
    Object for storing polynonmial mapping coefficients and reverse mapping as well
    
    The mapping can be symmetric, with a polynomial in r**2, or non symmetric with a polynomial in (x,y) 
    it can be by wavelength (a polynomial is build for each wavelength) 
    or the wavelength is another variable: r**2,w or x,y,w
    
    The methods are: 
        
        set: to compute the polynomial and reverse polynomial coeffs
       
        map: to map coordinates by computing the polynomial values
        
        inv_map: to map in reverse direction by computing the inverse polynomial values
        
        save, restore: to save and restore the object on file 
    '''
    
    def __init__(self, filename=None, symmetric=False, wavelength=None):
        '''
        Initialize the Mapping object.
        if file is given the object is restored from file.
        else
            the object is  configured with symmetric and wavelength arguments:
                
            symmetric: default to False, if true the polynome will be in r**2 instead of x,y
            wavelength: value or list or array of wavelength, 
                        a polynomial will be computed for each wavelength.
                        If not set: the polynomial will be in x,y,w 
        '''
        if filename is not None:
            self.restore(filename)
        else:
            self.symmetric = symmetric # symmetruc fit in r2
            if wavelength is not None:
                self.w = np.array(wavelength)
            else:
                self.w = None


    def set(self, w, xin, yin, xout, yout, deg=[2,2]):
        
        xin = np.array(xin).ravel()
        yin = np.array(yin).ravel()
        xout = np.array(xout).ravel()
        yout = np.array(yout).ravel()
        w = np.array(w).ravel()
        
        if self.symmetric:
            r2 = xin*xin + yin*yin
            rho2 = xout*xout + yout*yout

        # mapping x,y by wavelength
        if self.w is not None:
            if self.symmetric:
                raise Exception("symmetric fit by wavelength not implemented") 
            else:
                self.mapping = []
                self.inv_mapping = []
                xy_out = np.array((xout, yout)) .reshape(2, -1).T
                xy_in = np.array((xin, yin)) .reshape(2, -1).T
                for wi in self.w:
                    idx = (w == wi)
                    if np.any(idx):
                        self.mapping.append(PolyFit((xin[idx], yin[idx]), xy_out[idx], deg))
                        self.inv_mapping.append(PolyFit((xout[idx], yout[idx]), xy_in[idx], deg))
                    else:
                        self.mapping.append(None)
                        self.inv_mapping.append(None)
        # mapping x,y,w    
        else:
            if self.symmetric:
                self.mapping = PolyFit((w, r2), rho2, deg, "2d_woct")
                self.inv_mapping = PolyFit((w, rho2), r2, deg,"2d_woct")
                
            else:
                xy_out = np.array((xout, yout)) .reshape(2, -1).T
                self.mapping = PolyFit((w, xin, yin), xy_out, deg, "3d")
        
                xy_in = np.array((xin, yin)) .reshape(2, -1).T
                self.inv_mapping = PolyFit((w, xout, yout), xy_in, deg, "3d")
        
        
    def map(self, w, x, y, inverse=False):
        '''
        map the input coordinates x,y,w to the output x,y
        '''
        x = np.array(x)
        y = np.array(y)
        w = np.array(w)

        if self.symmetric:
            r2 = x**2 + y**2
            w, r2 = np.broadcast_arrays(w, r2)
        else:
            w, x, y = np.broadcast_arrays(w, x, y)

        if inverse:
            mapping = self.inv_mapping
        else:
            mapping = self.mapping
        

        # mapping by wavelength
        if self.w is not None:
            if not np.all(np.in1d(w, self.w)):
                raise Exception("wavelength given for mapping not valid")
            
            if self.symmetric:
                        raise Exception("symmetric fit by wavelength not implemented")
            
#            if self.w.size == 1:
#                if not np.all(w == self.w):
#                    raise Exception("wavelength given for mapping not valid")
#                xy_mapped = mapping[0]((x, y))
            
#            else:
#                # interpol between 1st 2 wavelength used for mapping
#                if inverse:
#                    raise Exception("Not implemented")
#                    
#                else:
#                    xy_w0 = mapping[0]((x, y))
#                    xy_w1 = mapping[1]((x, y))
#                    alpha = (w - self.w[0])/(self.w[1] - self.w[0])
#                    xy_mapped = xy_w0*(1-alpha) + xy_w1*alpha

                
            # group by wavelength
            xy_mapped = np.empty((2,) + w.shape)
            for wi, m in zip(self.w, mapping):
                idx = (w == wi)
                if np.any(idx):
                    if self.symmetric:
                        raise Exception("symmetric fit by wavelength not implemented")
#                        scale2 = m(r2[idx], overy=True)
#                        scale = np.sqrt(scale2)
#                        xy_mapped[0, idx] = -x * scale
#                        xy_mapped[1, idx] = -y * scale
                    else:
                        xy_mapped[:, idx] = m((x[idx], y[idx]))
    
        # mapping x,y,w
        else:
            if self.symmetric:
                scale2 = mapping((w, r2), overy=True)
                scale = np.sqrt(scale2)
                xy_mapped = np.empty((2,) + w.shape)
                xy_mapped[0] = -x * scale
                xy_mapped[1] = -y * scale                
            else:
                xy_mapped = mapping((w, x, y))
            
        return xy_mapped


    def inv_map(self, w, x, y):
        return self.map(w, x, y, inverse=True)
    
    
    def save(self, filename='mapping.pkl'):

        data = (self.mapping, self.inv_mapping, self.w, self.symmetric)
        with  open(filename,'wb') as file_handle:            
            pickle.dump(data, file_handle, protocol=2) # to be compatible with python 2


    def restore(self, filename='mapping.pkl'):
        
        data = load_pickle(filename)
            
        # unpack    
        self.mapping, self.inv_mapping, self.w, self.symmetric = data


    def linear(self, center=(0., 0.), w=None, inverse=False):
        '''
        return a linear Mapping at given position using local platescale at position
        '''    
        lin_map = Mapping(symmetric=self.symmetric, wavelength=self.w)

        if inverse:
            mapping = self.inv_mapping
        else:
            mapping = self.mapping
            
        if self.w is None:
            mapping = [mapping]
            ww = [w]
        else:
            ww = self.w
        
        lin_map.mapping = []
        for m, w in zip(mapping, ww):
            cm = self.map(w, center[0], center[1]) 
            if m.type == "2d":
                J = m.jacobian((center[0], center[1]))
            elif m.type == "2d_woct":
                r2 = center[0]**2 + center[1]**2
                xy = center[0]*center[1]
                dscale2 = m.jacobian((w, r2))
                scale = cm[0]/center[0]
                dscale = dscale2[1] / scale # 2 d/dr2 scale(w,r2) 
                J = np.empty((2,2))
                J[0, 0] = center[0]**2 * dscale + scale
                J[1, 0] = xy   * dscale
                J[0, 1] = xy   * dscale
                J[1, 1] = center[1]**2 * dscale - scale
            elif m.type == "3d":
                J = m.jacobian((w, center[0], center[1]), axis=[1,2])

            coeffs = np.array([[[cm[0], J[0,0]], [J[0,1], 0.]], 
                               [[cm[1], J[1,0]], [J[1,1], 0.]]])
            p = PolyFit(deg=[1,1], coeffs=coeffs, nbset=2)
            lin_map.mapping.append(p)
                        
        if self.w is None:
            lin_map.mapping = lin_map.mapping[0]
            
        return lin_map

        
    def plot(self, x=np.linspace(-12, 12, 21), y=np.linspace(-6, 6, 21), 
             center=(0., 0.), w=None, inverse=False):

        from matplotlib import pyplot as plt
             
        linear = self.linear(center, w, inverse)

#        if inverse:
#            mapping = self.inv_mapping
#            lin_mapping = linear.inv_mapping
#        else:
#            mapping = self.mapping
#            lin_mapping = linear.mapping

        if self.w is None:
#            mapping = [mapping]
#            lin_mapping = [lin_mapping]
            ww = [w]
        else:
            ww = self.w
                
        for w in ww:
            gx, gy = np.meshgrid(x,y)        
            # mapped coords:
            xy = self.map(w, gx, gy, inverse)        
            # ref coords (using platescale at center)
            xy0 = linear.map(w, gx, gy)
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for ii in range(xy0.shape[1]):
                plt.plot(xy0[0,ii,:], xy0[1,ii,:],'k:')
            for jj in range(xy0.shape[2]):
                plt.plot(xy0[0,:,jj], xy0[1,:,jj], 'k:')
                
#            dx = xy[0,:,:] - xy0[0,:,:]
#            dy = xy[1,:,:] - xy0[1,:,:]
#            
#            xc = xy0[0,:,:] + 10.*dx
#            yc = xy0[1,:,:] + 10.*dy
            xc = xy[0]
            yc = xy[1]
            for ii in range(xy.shape[1]):
                plt.plot(xc[ii,:], yc[ii,:],'b')
            for jj in range(xy.shape[2]):
                plt.plot(xc[:,jj], yc[:,jj], 'b')
                
            ax.text(0.5, + 0.02, 'difference exagerated 10 times', color='b',
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            ax.set_title('mapping at waevelength={:.4f}'.format(w))
            ax.set_xlabel('x')
            ax.set_ylabel('y')

        plt.show()
            
if __name__ == '__main__':
    pass
