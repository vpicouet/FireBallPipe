#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 09:22:09 2018

@author: dvibert
"""

from __future__ import division, print_function

import os
import numpy as np
from astropy.table import Table
from pkg_resources  import resource_filename
import matplotlib.pyplot as plt
from .mapping import Mapping

    #%%
def map_mask_detector(table, deg=[5, 2], bywave=True):
    '''
    From a table of observed mask slits (or holes) on the detector, 
    determine the mapping from science mask coord in mm to detector coord in mm
    and compute as well, the mask center position on the detector for each wavelength
            
    return a Mapping object and the array of centers  .
    '''    

    # list of wavelength
    subt = table[table['wavelength'] > 0]
    
    wset = np.unique(subt['wavelength'])
    
    wset.sort()
    
    if bywave:
        mapping = Mapping(wavelength=wset)
    else:
        mapping = Mapping()
    
    
    mapping.set(subt['wavelength'], subt['xmask'], subt['ymask'], subt['X_IMAGE'], subt['Y_IMAGE'] , deg)
        
    # compute residuals
    pred = mapping.map(subt['wavelength'], subt['xmask'], subt['ymask'])
    dx = subt['X_IMAGE'] - pred[0]
    dy = subt['Y_IMAGE'] - pred[1]
    for w in wset:
        flagw = (subt['wavelength'] == w)
        print("residuals of DIRECT mapping at wavelength: {}".format(w))
        print("mean: {}, {}".format(
                dx[flagw].mean(), dy[flagw].mean()))
        print("max:  {}, {}".format(
                np.abs(dx[flagw]).max(), np.abs(dy[flagw]).max()))
        print("rms:  {}, {}".format(
                np.sqrt(np.square(dx[flagw]).mean()), np.sqrt(np.square(dy[flagw]).mean())))

    pred = mapping.inv_map(subt['wavelength'], subt['X_IMAGE'], subt['Y_IMAGE'])
    dx = subt['xmask'] - pred[0]
    dy = subt['ymask'] - pred[1]        
    for w in wset:
        flagw = (subt['wavelength'] == w)
        print("residuals of INVERSE mapping at wavelength: {}".format(w))
        print("mean: {}, {}".format(
                dx[flagw].mean(), dy[flagw].mean()))
        print("max:  {}, {}".format(
                np.abs(dx[flagw]).max(), np.abs(dy[flagw]).max()))
        print("rms:  {}, {}".format(
                np.sqrt(np.square(dx[flagw]).mean()), np.sqrt(np.square(dy[flagw]).mean())))

        #print("at wavelength: {}".format(w))
        #print("maximum residual of inverse mapping along x, y (in mm): {}, {}".format(
        #        np.abs(dx[flagw]).max(), np.abs(dy[flagw]).max()))
        
    # get center
    centers = mapping.map(wset, 0., 0.)
                                       
    return mapping, centers



#%%
def recompute_mask_pos(table, mask):
    '''
    Add the xmask and ymask missing column in the table
    for old version of focustest, when windowing was effective the xmask and ymask 
    were not saved in the table
    '''
    __file__ = '/Users/Vincent/Github/FireBallPipe/Calibration/mapping_mask_det.py'

    try:
        Target_dir = resource_filename('Calibration', 'Targets')
    except:
        Target_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Targets')


    if mask == 'F1':
        slit_location = os.path.join(Target_dir, 'targets_F1.txt')
    if mask == 'grid':
        slit_location = os.path.join(Target_dir, 'grid_mask.txt')
    if mask == 'F2':
        raise(TypeError,"Need to have this catalogs for F2 and F3")#slit_location = os.path.join(Target_dir, 'targets_F2.txt')
    if mask == 'F3':
        raise(TypeError,"Need to have this catalogs for F2 and F3")#slit_location = os.path.join(Target_dir, 'targets_F2.txt')
    if mask == 'F4':
        slit_location = os.path.join(Target_dir, 'targets_F4.txt')
    
    slitfile = slit_location
    if 'targets_F1' in slit_location:
        print('Computing predicted location for F1')
        slits = Table.read(slitfile, format='ascii', delimiter='\t')
        # remove multi object slits
        idok = (slits['slit_length_right'] != 0) &  (slits['slit_length_left'] !=0)
        slits = slits[idok]
        xmask = slits['xmask'] + (slits['slit_length_right'] - slits['slit_length_left'])/2.
        ymask = slits['ymask'] + slits['offset']
        z = slits['z'] 
        internalCount = slits['Internal-count']
    
    elif 'grid_mask' in slit_location:
        print('Computing predicted location for grid mask')
        slits =  Table.read(slitfile, format='ascii', delimiter='\t')
        xmask = slits['x']
        ymask = slits['y']
        #x = xmask
        #y = -ymask
    
    elif ('F2' in slit_location) or ('F3' in slit_location) or ('F4' in slit_location) or ('Tilted' in slit_location):
        print('Computing predicted location for science mask')
        slits =  Table.read(slitfile, format='ascii', delimiter='\t')
        xmask = slits['xmm']
        ymask = slits['ymm']
        internalCount = slits['Internal-count']
        #mainID = slits['#main_id']
        if not ('Tilted' in slit_location):
            z = slits['Z'] 

    table['xmask'] = 0.
    table['ymask'] = 0.
    flag = table['id_slit'] >= 0 
    sid = table['id_slit'][flag]             
    table['xmask'][flag] = xmask[sid]
    table['ymask'][flag] = ymask[sid]
    
    return xmask, ymask, internalCount, z


def create_DS9regions(xim, yim, radius=20, save=True, savename="test", form=['circle'], color=['green'], ID=None):#of boxe
    """Returns and possibly save DS9 region (circles) around sources with a given radius
    """
    
    regions = """# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image
"""
    r = radius    
    for i in range(len(xim)):
        if form[i] == 'box':
            rest = '{:.2f},{:.2f})'.format(600, r)
        elif form[i]=='circle':
            rest = '{:.2f})'.format(r)
        elif form[i] == 'bpanda':
            rest = '0,360,4,0.1,0.1,{:.2f},{:.2f},1,0)'.format(r, 2*r)
        elif form[i] == 'annulus':
            rtab = np.linspace(0, r, 10)
            rest = '{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f})'.format(*list(rtab))  
        rest += ' # color={}'.format(color[i])
        for j, (x, y) in enumerate(np.nditer([xim[i], yim[i]])):
            regions += '{}({:.2f},{:.2f},'.format(form[i], x+1, y+1) + rest
            if ID is not None:
                regions += ' text={{{}}}'.format(ID[i][j])
            regions +=  '\n'   

    if save:
        with open(savename+'.reg', "w") as text_file:
            text_file.write(regions)        
        print(('Mapped region file saved at: ' +  savename + '.reg'))
        return 

    #%%

if __name__ == '__main__':
    path = '/data/FireBall/FTS-06-2018/180605/'
    path="/Volumes/ExtremePro/LAM/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180605/"
    
    filename = 'image-120-124-Zinc-Dark-substarcted119-stack_table.csv'

    t=Table.read(path+filename, format='csv')

    xmask, ymask, internalCount, z = recompute_mask_pos(t, 'F1')
    
    mappings, centers = map_mask_detector(t)
    print(centers)
    
    #mappings.plot()

    mappings.save('mapping-mask-det-F1.pkl')
    
    # create mapped region file
    wcolors = {0.20255:'blue', 0.20619:'green', 0.21382:'red'}
    colors = [wcolors[i] for i in mappings.w]
    ID =[internalCount for i in range(len(mappings.w))] 
    xdetpix = []
    ydetpix = []
    for w in mappings.w:
        xydetpix = mappings.map(w, xmask, ymask)
        xdetpix.append(xydetpix[0])
        ydetpix.append(xydetpix[1])
    create_DS9regions(xdetpix, ydetpix, form=['bpanda']*3, radius=10, save=True, 
                          savename=path + 'image-120-124-Zinc-Dark-substarcted119-stack_mapped', color = colors, ID=ID)
