#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 10:20:08 2018

@author: dvibert
"""

import numpy as np
from numpy.polynomial import polynomial as poly

from .polynomial import  polyfit3d, polyfit2d, polyfit2d_wocty

polyfit_func = {"2d": polyfit2d, "3d": polyfit3d, "2d_woct": polyfit2d_wocty}

class PolyFit(object):
        
    def __init__(self, xyz=None, val=None, deg=[2, 2], type='2d', 
                 coeffs=None, nbset=None):

        self.deg = deg
        self.dim = len(deg)
        self.type = type
        #self.fit = polyfit_func[type]

        if xyz is not None:
            self.set(xyz, val)

        if coeffs is not None:
            self.coeffs = coeffs
            self.nbset = nbset
            # check dims
            if self.coeffs.ndim != (self.dim  + (self.nbset > 1)):
                raise Exception("coeffs argument has wrong dimension")
                                
    def __call__(self, xyz, overy=False):
        if (self.type ==  "2d"):
            x, y = xyz
            return poly.polyval2d(x, y, self.coeffs)
        elif (self.type ==  "2d_woct"):
            x, y = xyz
            if not overy:
                return y * poly.polyval2d(x, y, self.coeffs)
            else:
                # return P(x,y)/y 
                return poly.polyval2d(x, y, self.coeffs) 
        else:
            x, y, z = xyz
            return poly.polyval3d(x, y, z, self.coeffs)
    

    def set(self, xyz, val):
        val = np.asarray(val, dtype=float)
        if (self.type ==  "2d") or (self.type == "2d_woct"):
            x, y = xyz
            self.coeffs = polyfit_func[self.type](x, y, val, self.deg)
        else:
            x, y, z = xyz
            self.coeffs = polyfit_func[self.type](x, y, z, val, self.deg)
        if val.ndim == 1:
            self.nbset = 1
        else: # val.ndim=2
            self.nbset = val.shape[1]
        return self

        
    def _compute_polyder(self):
        self._Dcoeffs = []
        for j in range(self.dim):
            if self.nbset == 1:
                d = poly.polyder(self.coeffs, axis=j)
            else:
                d = []
                for i in range(self.nbset):
                    d.append(poly.polyder(self.coeffs[...,i], axis=j))
                d = np.array(d)
                d = np.rollaxis(d, 0, d.ndim)
            self._Dcoeffs.append(d)
        return self._Dcoeffs


    @property
    def Dcoeffs(self):
        try:
            return self._Dcoeffs
        except AttributeError:
            return self._compute_polyder()

            
    def jacobian(self, xyz, axis=None):
        if axis is None:
            axis = np.arange(self.dim)
        else:
            axis = np.array(axis, copy=False)
        J = []
        for j in axis.flat:
            if (self.type ==  polyfit2d) or (self.type ==  polyfit2d_wocty):
                x, y = xyz
                d = poly.polyval2d(x, y, self.Dcoeffs[j])
            else:
                x, y, z = xyz
                d = poly.polyval3d(x, y, z, self.Dcoeffs[j])
            if self.Dcoeffs[j].ndim > self.dim: d = np.rollaxis(d, 0, d.ndim)
            J.append(d)
        J = np.array(J, copy=False)
        J = np.rollaxis(J, 0, J.ndim)
        if axis.ndim == 0: # axis is scalar 
            J = np.squeeze(J, axis=-1)
        return J   
    
    
    
    
