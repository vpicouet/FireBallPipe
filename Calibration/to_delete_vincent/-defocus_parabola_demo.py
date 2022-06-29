#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  2 07:03:11 2018

@author: dvibert
"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from focustest import build_defocus

gx, gy = np.meshgrid(np.linspace(-12, 12, 20), np.linspace(-6, 6, 20))

gx_det = -0.87 * gy
gy_det = -gx

defocus = build_defocus()
a, b = 11, 6
det_offset, det_xslope, det_yslope = 0., 0., 0. 
fwhm =  defocus(np.array([gx, gx_det, gy_det]), det_offset, det_xslope, det_yslope, a, b)

def plot_defocus():
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(121, projection='3d')
    ax.set_title('defocus')
    ax.set_xlabel('x detector pix')
    ax.set_ylabel('y detector pix')
    ax.set_zlabel('FWHM size pix')
    ax.plot_wireframe(gx_det, gy_det, fwhm, alpha=.15)
    
    ax = fig.add_subplot(122)
    ax.set_title('defocus')
    ax.set_xlabel('x detector pix')
    ax.set_ylabel('y detector pix')
    ax.scatter(gx_det, gy_det, s=fwhm*10)