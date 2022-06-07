# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates
from astropy import units as u
from astropy.table import Table

#from guider2UV import Guider2UV
from MaskAstrometry import LocalScienceMaskProjector
### check science mask magnif/dist versus JG calc



target_file = "/home/dvibert/ownCloud/FIREBALL/Target_selection_meeting_NY_20170405/targets_F1.txt"
target_file = '/Users/Vincent/Github/DS9functions/DS9FireBall/Targets/targets_F1.txt'
target_tab = Table.read(target_file, format='ascii')

plt.figure()
plt.plot(target_tab['ra'], target_tab[ 'dec'], '+ ')

a =42.26134
b=-3.154411e-3
c=-1.117322

Field_center = coordinates.SkyCoord(32.19*u.deg, -5.688*u.deg)
FieldP = LocalScienceMaskProjector(Field_center)

r = np.linspace(0, 1/6., 100) # 0 -10arcmin
gr = FieldP.radial_magnification(r)

plt.figure()
plt.plot(r, gr)

def correctFORplatescale(radius):
    #X=MX
    #Y=MY
    # if we adopt plate scale from projection above 
    # which is 0.0235  --> 42.5532 mm/deg = 11.8203 um/"
    # we need to divide out the assumed plate scale
    ps0=11.8203
    #radius=np.sqrt((X**2)+(Y**2))
    psnew=-1.7584761e-4*radius**2-(6.0880439e-06*radius)+1.1739247e01
    #return psnew*X/ps0,psnew*Y/ps0
    return psnew*radius/ps0

def platescale(radius):
    #X=MX
    #Y=MY
    # if we adopt plate scale from projection above 
    # which is 0.0235  --> 42.5532 mm/deg = 11.8203 um/"
    # we need to divide out the assumed plate scale
    #radius=np.sqrt((X**2)+(Y**2))
    psnew=-1.7584761e-4*radius**2-(6.0880439e-06*radius)+1.1739247e01
    #return psnew*X/ps0,psnew*Y/ps0
    return psnew*radius 


gr2 = correctFORplatescale(r/0.0235)
plt.figure()
plt.plot(r, gr2)
