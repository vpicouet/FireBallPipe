#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 08:49:12 2017

@author: dvibert
"""

from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy import wcs, coordinates, Table
from astropy import units as u
from astropy.table import Table
from focustest2 import Image

from LocalMaskProjector import LocalGuiderMaskProjector, LocalScienceMaskProjector

f3_astrometry_filename = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/perf3.0-09-01/sky_00_astrometry/stack475788.new"
f3_astrometry_filename = "/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/perf3.0-09-01/sky_00_astrometry/stack475788.new"
f1_astrometry_filename = "/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/perf3.0-09-01/sky_00_astrometry/stack500513.new"


hdu = fits.open(f1_astrometry_filename)
w  = wcs.WCS(hdu[0].header)

#field center
center = w.wcs.crval*u.deg
center = coordinates.SkyCoord(center[0], center[1])

# get field rotation from wcs header
field_center_radec = w.wcs_pix2world([w.wcs.crpix], 1)
north_point_radec = field_center_radec + [0., 1.]
north_point_pix = w.wcs_world2pix(north_point_radec, 1) - w.wcs.crpix
rot_angle = np.arctan2(north_point_pix[0,0], north_point_pix[0,1])*180/np.pi

localframe = coordinates.SkyOffsetFrame( origin=center, rotation=-rot_angle*u.deg)


# check with M13
#im = Image(filename = f3_astrometry_filename, plot = True, stack = True, stack_image = False, py = True, subDark = False, verbose = True, quick = True, Type = 'Detector', Fits = True)
m13tab = Table.read('stack475788.starpos.fits')

#fig = plt.figure()
#fig.add_subplot(111, projection=w)
#plt.imshow(im.image, origin='lower')
#plt.xlabel('RA')
#plt.ylabel('Dec')


plt.figure()
#plt.imshow(hdu[0].data, origin='lower')
plt.plot(m13tab['xcentroid'], m13tab['ycentroid'],'+r')

m13_radec = w.wcs_pix2world(m13tab['xcentroid'], m13tab['ycentroid'], 1)
m13_coords = coordinates.SkyCoord(m13_radec[0]*u.deg, m13_radec[1]*u.deg)

plt.figure()
plt.plot(m13_coords.ra, m13_coords.dec, '+')

# rotate
m13_coords_local = m13_coords.transform_to(localframe)
plt.figure()
plt.plot(m13_coords_local.lon, m13_coords_local.lat, '+')

# test the reverse rotation
m13_coords_local.transform_to(coordinates.ICRS)


