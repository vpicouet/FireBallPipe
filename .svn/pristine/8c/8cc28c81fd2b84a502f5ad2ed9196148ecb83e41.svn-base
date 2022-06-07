#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 16:43:20 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy.table import Table
from guider2UV import Guider2UV

path_SCGUI01 = '/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/'

G2UV= Guider2UV(filename=path_SCGUI01  + 'Guider2UV_F1.pkl')

path_SCGUI03 = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI03/Pattern/"

    
    
# F1 Field, predict star pattern pos
##########################################
new_F1_star_pattern_filename = path_SCGUI03  + 'F1.csv'

new_F1_star_pattern_tab = Table.read(new_F1_star_pattern_filename)  
new_F1_star_pattern_tab.sort('xcentroid')
new_F1_star_pattern_pos = np.array([new_F1_star_pattern_tab['xcentroid'], new_F1_star_pattern_tab['ycentroid']]).T
starid = 2
slit_pos = np.array([0.5525186, -4.7601449]) # slit 30
pred = G2UV.pattern_in_slit(starid, slit_pos, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred)
#[array([ 380.11130579,  429.91799389,  894.13619392,  939.14423429]), array([ 505.27525037,   -3.75743313,  549.3900237 ,   51.38926488])]

# slit 11
starid = 2
target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/targets_F1.txt'
F1 = Table.read(target_filename, format='ascii')
slit_pos2 =   np.array([F1[F1['Internal-count']=='11']['xmm'][0], F1[F1['Internal-count']=='11']['ymm'][0]])
pred2 = G2UV.pattern_in_slit(starid, slit_pos2, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
print(pred2)
#[array([  832.04175069,   860.62293966,  1329.57788462,  1355.35321079]), array([ 1313.57801462,   817.23204211,  1334.44249052,   848.2060185 ])]#




# F2 Field, predict star pattern pos
##########################################
new_F2_star_pattern_filename = path_SCGUI03  + 'F2.csv'

new_F2_star_pattern_tab = Table.read(new_F2_star_pattern_filename)  
new_F2_star_pattern_tab.sort('xcentroid')
new_F2_star_pattern_pos = np.array([new_F2_star_pattern_tab['xcentroid'], new_F2_star_pattern_tab['ycentroid']]).T
starid = 2

# slit 16
target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/targets_F2.txt'
F2 = Table.read(target_filename, format='ascii')
slit_pos1 =   np.array([F2[F2['Internal-count']=='16']['xmm'][0], F2[F2['Internal-count']=='16']['ymm'][0]])
pred = G2UV.pattern_in_slit(starid, slit_pos1, world=False, FourStarsGuider_pos=new_F2_star_pattern_pos)
print(pred)
# [array([  807.2533361 ,   840.70961313,  1305.02947189,  1336.53788126]), 
#  array([ 1085.48705905,   588.8923821 ,  1112.52765071,   626.14949714])]

# slit 30 
starid = 2
slit_pos2 =   np.array([F2[F2['Internal-count']=='30']['xmm'][0], F2[F2['Internal-count']=='30']['ymm'][0]])
pred2 = G2UV.pattern_in_slit(starid, slit_pos2, world=False, FourStarsGuider_pos=new_F2_star_pattern_pos)
print(pred2)
#[array([  447.67578982,   496.02458309,   958.53905731,  1003.55546651]),
# array([ 527.87252833,   21.28120254,  570.99855443,   75.39609185])]



# F3 Field, predict star pattern pos
##########################################
new_F3_star_pattern_filename = path_SCGUI03  + 'F3.csv'

new_F3_star_pattern_tab = Table.read(new_F3_star_pattern_filename)  
new_F3_star_pattern_tab.sort('xcentroid')
new_F3_star_pattern_pos = np.array([new_F3_star_pattern_tab['xcentroid'], new_F3_star_pattern_tab['ycentroid']]).T
starid = 2

# slit OVI1
target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/targets_F3.txt'
F3 = Table.read(target_filename, format='ascii')
slit_pos1 =   np.array([F3[F3['Internal-count']=='OVI1']['xmm'][0], F3[F3['Internal-count']=='OVI1']['ymm'][0]])
pred = G2UV.pattern_in_slit(starid, slit_pos1, world=False, FourStarsGuider_pos=new_F3_star_pattern_pos)
print(pred)
#[array([  754.18076508,   783.27160284,  1251.4357524 ,  1277.47992686]),
# array([ 1226.71999722,   727.71885158,  1247.56760069,   759.58669904])]



# F4 Field, predict star pattern pos
##########################################
new_F4_star_pattern_filename = path_SCGUI03  + 'F4.csv'

new_F4_star_pattern_tab = Table.read(new_F4_star_pattern_filename)  
new_F4_star_pattern_tab.sort('xcentroid')
new_F4_star_pattern_pos = np.array([new_F4_star_pattern_tab['xcentroid'], new_F4_star_pattern_tab['ycentroid']]).T
starid = 2

# slit #8
target_filename = '/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/targets_F4.txt'
F4 = Table.read(target_filename, format='ascii')
slit_pos1 =   np.array([F4[8]['xmm'], F4[8]['ymm']])
pred = G2UV.pattern_in_slit(starid, slit_pos1, world=False, FourStarsGuider_pos=new_F4_star_pattern_pos)
print(pred)
#[array([  985.19057   ,  1014.69816358,  1473.68704015,  1501.67235076]),
# array([ 1202.02227071,   708.56125385,  1225.05486884,   741.9465758 ])]



