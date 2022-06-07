#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 12:13:16 2018

@author: Vincent
"""

from __future__ import division


import os, sys

path = os.getcwd()
sys.path.append(path + '/Calibration_SW/')
sys.path.append('/Users/Vincent/Documents/FireBallIMO/')
from IPython.display import Image as imdisplay


from focustest import *#PlotFocus2DGuider

FUV = Table.read('/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180821/thrufocus/TotalThroughfocus.csv')
F1 = Table.read('/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/F1_119/TotalThroughfocus.csv')
F2 = Table.read('/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/F2_-161/TotalThroughfocus.csv')
F3 = Table.read('/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/F3_-121/TotalThroughfocus.csv')
F4 = Table.read('/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/F4_159/TotalThroughfocus.csv')

def DefineBestActuator(table):
    mean = np.nanmean(np.array(np.array(table['Best sigma','Best EE50','Best EE80','Best Maxpix', 'Best Varpix']).tolist()),axis=1)
    var = np.nanstd(np.array(np.array(table['Best sigma','Best EE50','Best EE80','Best Maxpix', 'Best Varpix']).tolist()),axis=1)
    table['MeanBestActuator'] = mean
    table['VarBestActuator'] = var
    return table
FUV = DefineBestActuator(FUV)
F1 = DefineBestActuator(F1)
F2 = DefineBestActuator(F2)
F3 = DefineBestActuator(F3)
F4 = DefineBestActuator(F4)

table=FUV
plt.plot()
#plt.plot(table['x'],table['Best EE80'],'x',label = 'EE80')
#plt.plot(table['x'],table['Best EE50'],'+',label = 'EE50')
#plt.plot(table['x'],table['Best sigma'],'.',label = 'sigma')
plt.errorbar(table['y'],table['MeanBestActuator'],  fmt='o',yerr = table['VarBestActuator'],label = 'Mean')
#plt.plot(table['x'],table['Best Maxpix'],'o',label = 'sigma')
#plt.plot(table['x'],table['Best Varpix'],'o',label = 'sigma')
plt.legend()
plt.show()

table = FUV
X161,Y161,Z161, ax161, Cguider161 = fit_quadratic_curve(table['x'],table['y'],table['MeanBestActuator'],sigma_z = table['VarBestActuator'], n=100,order=1)
plt.show()