from __future__ import division, print_function
#from builtins import range, zip

import numpy as np
from polynomialD import polyfit2d
import numpy.polynomial.polynomial as poly

class PolyFit2d(object):
 
    def __init__(self, x, y, val, deg=[2, 2], w=None):
        self.deg = deg
        val = np.array(val, dtype=float)
        self.coeffs, self.stat = polyfit2d(x, y, val, deg, full = True, w=w)
        if val.ndim == 1:
            self.nbset = 1
        else: # val.ndim=2
            self.nbset = val.shape[1]
                                
    def __call__(self, x, y):
        return poly.polyval2d(x, y, self.coeffs)
           
    def _compute_polyder(self):
        self._Dcoeffs = []
        for j in range(2):
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
            
    def jacobian(self, x, y, axis=np.arange(2)):
        axis = np.array(axis, copy=False)
        J = []
        for j in axis.flat:
            d = poly.polyval2d(x, y, self.Dcoeffs[j])
            if self.Dcoeffs[j].ndim > 2: d = np.rollaxis(d, 0, d.ndim)
            J.append(d)
        J = np.array(J, copy=False)
        J = np.rollaxis(J, 0, J.ndim)
        if axis.ndim == 0: # axis is scalar 
            J = np.squeeze(J, axis=-1)
        return J   
