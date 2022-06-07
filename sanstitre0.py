#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 14:34:49 2019

@author: dvibert
"""

import numpy as np
from Calibration.multiple_threshold import process_image

from matplotlib import pyplot as plt
from astropy.convolution import Gaussian2DKernel

# check background estimation bias from averaging the erlang distribution
rnoise = 109 # electrons
gain = 470 # em amplification gain

kernel = Gaussian2DKernel(3)
Nsamp = kernel.array.size
input_mean = np.linspace(0.1, 5, 100)
bias = np.zeros(100)

Nreal = 3000 # to estim bias

for i in range(100):
    input_poisson =  np.random.poisson(input_mean[i], (Nsamp, Nreal))
    read_noise_real = np.random.normal(scale=rnoise, size=(Nsamp, Nreal) )
    amplification_real = np.random.gamma(input_poisson, gain)
    real = amplification_real + read_noise_real
    #estim = real.mean(axis=0)
    estim = np.average(real, weights=kernel.array.ravel(), axis=0)
    bias[i] = estim.mean()/gain - input_mean[i]

#plt.figure()
#plt.hist(estim/gain)
#print(estim.mean()/gain)
#print(np.sqrt(estim.var())/gain)

plt.figure()
plt.plot(input_mean, bias)