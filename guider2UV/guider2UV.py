#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 14:20:31 2017

@author: dvibert
"""
from __future__ import division, print_function

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs, coordinates
from astropy import units as u
from astropy.wcs.utils import proj_plane_pixel_scales
#from scipy.optimize import curve_fit

try:
   import cPickle as pickle # for python 2 only CPickle is faster
except:
   import pickle

from .MaskAstrometry import LocalGuiderMaskProjector, LocalScienceMaskProjector

def spherical_shift(radec, ra_c, dec_c):
    '''Given a line of sight (ra_c, dec_c), the function converts 
    right ascenssion and declination coordinates in a local frame centered
    on the line of sight'''
    ra, dec = radec
    a = coordinates.SkyCoord(ra, dec, unit="deg")
    F1_center = coordinates.SkyCoord(ra_c*u.deg, dec_c*u.deg)
    F1_Frame = coordinates.SkyOffsetFrame(origin=F1_center)
    a2 = a.transform_to(F1_Frame)
    return np.array([a2.lon.deg, a2.lat.deg]).ravel()



def register_UV_star_simple(pattern, Stars_UV_Visibe, Stars_id, UVstar_id=3):
    delta_lon = Stars_UV_Visibe.lon.deg - pattern.lon.deg[Stars_id] 
    delta_lat = Stars_UV_Visibe.lat.deg - pattern.lat.deg[Stars_id] 
    delta_lon = delta_lon.mean()
    delta_lat = delta_lat.mean()
    UVstar_guider_coord_lon = pattern[UVstar_id].lon.deg + delta_lon
    UVstar_guider_coord_lat = pattern[UVstar_id].lat.deg + delta_lat
    return np.array([UVstar_guider_coord_lon, UVstar_guider_coord_lat])

 
def diff_skycoord(s1, s2):
    return np.array((s2.lon.deg - s1.lon.deg , s2.lat.deg - s1.lat.deg))


def coord_list_to_array(coord):

    n = len(coord)
    coord_arr = np.zeros((n, 2))
    
    for i in range(n):
        c = coord[i]
        coord_arr[i] = np.array([ c.lon.deg, c.lat.deg]).T

    return coord_arr

   
def fit_model(coord, coord_obs, gamma=True, ytilt=False, weight=None):
    
    coord_arr = coord_list_to_array(coord)
    coord_obs_arr = coord_list_to_array(coord_obs)
        
    #delta = coord_obs_arr - coord_arr
    delta = coord_arr - coord_obs_arr
    data = delta.T.reshape(-1)
    
    n = len(coord)
    if gamma:
        if not ytilt:
            # with gamma
            print('Fitting rotation, translation and magnification')
            row_x = np.hstack((coord_arr*[1.,-1], np.ones((n,1)), np.zeros((n,1)) )) # xn - x =  x dgamma - y theta + dx  
            row_y = np.hstack((coord_arr[:,::-1], np.zeros((n,1)), np.ones((n,1)) )) # yn - y =  y dgamma + x theta      + dy 
            mat = np.vstack((row_x, row_y))
            if weight is not None:
                mat = np.diag(np.sqrt(weight)).dot(mat)
            matinv =  np.linalg.pinv(mat)
            sol = matinv.dot(data)
            gamma = 1 + sol[0]
            theta_rad = sol[1]
            thetay_rad = 0.
            deltax = sol[2]
            deltay = sol[3]
            theta = theta_rad*180/np.pi*60 #arcmin
            print("gamma: {}\ntheta: {} arcmin\ndx: {} arcsec\ndy: {} arcsec".format(gamma, theta, deltax*3600, deltay*3600))
            covar = matinv.dot(matinv.T)
            # accuracy, assuming 1 arcsec measurement error
            print("variances: {}\n".format(np.sqrt(np.diag(covar))/3600*[1, 180/np.pi*60, 3600, 3600])) #
        else:
            # with gamma, ytilt
            print('Fitting rotation, translation and magnification')
            row_x = np.hstack((coord_arr[:,[0,1,1]]*[1.,-1,-1], np.ones((n,1)), np.zeros((n,1)) )) # xn - x =  x dgamma - y theta - y tilt + dx  
            row_y = np.hstack((coord_arr[:,::-1], np.zeros((n,2)), np.ones((n,1)) ))               # yn - y =  y dgamma + x theta                              + dy 
            mat = np.vstack((row_x, row_y))
            if weight is not None:
                mat = np.diag(np.sqrt(weight)).dot(mat)
            matinv =  np.linalg.pinv(mat)
            sol = matinv.dot(data)
            gamma = 1 + sol[0]
            theta_rad = sol[1]
            thetay_rad = sol[2]                  
            deltax = sol[3]
            deltay = sol[4]
            theta = theta_rad*180/np.pi*60 #arcmin
            thetay = thetay_rad*180/np.pi*60 #arcmin
            print("gamma: {}\ntheta: {} arcmin\ntheta_y: {} arcmin\ndx: {} arcsec\ndy: {} arcsec".format(gamma, theta, thetay, deltax*3600, deltay*3600))
            covar = matinv.dot(matinv.T)
            # accuracy, assuming 1 arcsec measurement error
            print("variances: {}\n".format(np.sqrt(np.diag(covar))/3600*[1, 180/np.pi*60, 180/np.pi*60, 3600, 3600])) #
            
    else:
        if not ytilt:
            # without gama
            print('Fitting rotation and translation')
            row_x = np.hstack((-coord_arr[:,[1]], np.ones((n,1)), np.zeros((n,1)) )) # xn - x =  - y theta + dx  
            row_y = np.hstack((coord_arr[:,[0]], np.zeros((n,1)), np.ones((n,1)) )) # yn - y =    x theta      + dy 
            mat = np.vstack((row_x, row_y))
            if weight is not None:
                mat = np.diag(np.sqrt(weight)).dot(mat)
            matinv =  np.linalg.pinv(mat)
            sol = matinv.dot(data)
            gamma = 1.
            theta_rad = sol[0]
            thetay_rad = 0.
            deltax = sol[1]
            deltay = sol[2]
            theta = theta_rad*180/np.pi*60 #arcmin
            print("theta: {} arcmin\ndx: {} arcsec\ndy: {} arcsec".format(theta, deltax*3600, deltay*3600))
            covar = matinv.dot(matinv.T)
            # accuracy, assuming 1 arcsec measurement error
            print("variances: {}\n".format(np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 3600, 3600]))
        else:
            # without gama, with ytilt
            print('Fitting rotation and translation')
            row_x = np.hstack((-coord_arr[:,[1,1]], np.ones((n,1)), np.zeros((n,1)) )) # xn - x =  - y theta - y tilt + dx  
            row_y = np.hstack((coord_arr[:,[0]], np.zeros((n,2)), np.ones((n,1)) )) # yn - y =    x theta                 + dy 
            mat = np.vstack((row_x, row_y))
            if weight is not None:
                mat = np.diag(np.sqrt(weight)).dot(mat)
            matinv =  np.linalg.pinv(mat)
            sol = matinv.dot(data)
            gamma = 1.
            theta_rad = sol[0]
            thetay_rad = sol[1]
            deltax = sol[2]
            deltay = sol[3]
            theta = theta_rad*180/np.pi*60 #arcmin
            thetay = thetay_rad*180/np.pi*60 #arcmin
            print("theta: {} arcmin\ntheta_y: {} arcmin\ndx: {} arcsec\ndy: {} arcsec".format(theta, thetay, deltax*3600, deltay*3600))
            covar = matinv.dot(matinv.T)
            # accuracy, assuming 1 arcsec measurement error
            print("variances: {}\n".format(np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 180/np.pi*60, 3600, 3600]))
            
        
    #residuals
    data_new = mat.dot(sol)
    delta_new = data_new.reshape((2,n)).T
    residuals = delta_new - delta
    print("residuals in arcsec:", residuals*3600)
    residual_max = np.abs(residuals).max(axis=0)*3600.
    residual_rms = np.sqrt(np.square(residuals).mean(axis=0))*3600
    print("max residual in EL,CE {:.1f}, {:.1f} arcsec".format( *list(residual_max)) )     
    print("mean residual in EL,CE {:.1f}, {:.1f} arcsec".format(*list(residual_rms) ))     

    fullsol = np.array((gamma-1., theta_rad, thetay_rad, deltax, deltay))
    
    #fullsol = -fullsol

    return fullsol, residuals

     
def plot_fit(coord, coord_obs, residuals, labels=None, sol=None):
    
    n = len(coord)
    
    coord_arr = coord_list_to_array(coord)
    coord_obs_arr = coord_list_to_array(coord_obs)
        
    delta = coord_obs_arr - coord_arr
    
    plt.figure(figsize=(8,4))
    plt.subplot(121)
    plt.axis('equal')
    plt.title("model versus mesure")
    plt.plot(coord_arr[:,0]*3600, coord_arr[:,1]*3600,'or')
    if labels is not None:
        for i in range(n):
            plt.text(coord_arr[i,0]*3600, coord_arr[i,1]*3600, labels[i], color='k')
    qv = plt.quiver(coord_arr[:,0]*3600, coord_arr[:,1]*3600, delta[:,0]*3600, delta[:,1]*3600)
    plt.quiverkey(qv, .8,.9,10, '10 arcsec', color='k')
    plt.ylim(plt.ylim()[::-1])
    plt.xlabel('EL arcsec')
    plt.ylabel(' - CE arcsec')
    delta_mean = delta.mean(axis=0)
    delta_rms = np.sqrt(np.square(delta - delta_mean).mean(axis=0))*3600
    delta_mean *= 3600
    legend =  "error mean in EL, CE {:.1f}, {:.1f} arcsec\n".format(*list(delta_mean))
    legend += "error rms in EL,CE {:.1f}, {:.1f} arcsec".format(*list(delta_rms))
    plt.text(-10,7, legend)
    
    plt.subplot(122)
    plt.axis('equal')
    plt.title("residual after model fit")
    plt.plot(coord_obs_arr[:,0]*3600, coord_obs_arr[:,1]*3600,'or')
    if labels is not None:
        for i in range(n):
            plt.text(coord_obs_arr[i,0]*3600, coord_obs_arr[i,1]*3600, labels[i], color='k')
    qv = plt.quiver(coord_obs_arr[:,0]*3600, coord_obs_arr[:,1]*3600, residuals[:,0]*3600, residuals[:,1]*3600)
    plt.quiverkey(qv, .8,.9,10, '10 arcsec', color='k')
    plt.ylim(plt.ylim()[::-1])
    plt.xlabel('EL arcsec')
    plt.ylabel(' - CE arcsec')
    
    legend = ''
    if sol is not None:
        dgamma, theta_rad, thetay_rad, deltax, deltay = sol
        gamma = 1. + dgamma
        theta = theta_rad*180/np.pi*60 #arcmin
        thetay = thetay_rad*180/np.pi*60 #arcmin
        legend = "rotation: {:.2f} arcmin\n".format(theta)
        legend += "perpendicularity: {:.2f} arcmin\n".format(thetay)
        legend += "magnification {:.4f}\n".format(gamma)
        legend += "deltax: {:.4f} arcsec\n".format(deltax*3600)
        legend += "deltay: {:.4f} arcsec\n".format(-deltay*3600)
    residual_rms = np.sqrt(np.square(residuals).mean(axis=0))*3600
    legend += "mean residual in EL,CE {:.1f}, {:.1f} arcsec".format(*list(residual_rms))
    plt.text(-10,7, legend)

    plt.show()


def add_hystcomp(moves):
    
    nmoves = moves.shape[0]

    # add reverse moves
    moves = np.vstack([moves, -moves[::-1]])
    
    
    # add 1st CE compensation:
    delta = 30.
    precomp1 = np.array([0., -delta])*np.sign(moves[0,1])
    
    hystcomp_moves = [precomp1, -precomp1]
    flags = [False, True]
    
    for i in range(nmoves):
        hystcomp_moves.append(moves[i])
        flags.append(True)
        
    # add reverse move precomp
    precomp_rev = np.array([0., -delta])*np.sign(moves[nmoves,1])
    hystcomp_moves.extend([precomp_rev, -precomp_rev])
    flags.extend([False,True])
    
    # add eventual precomp in reverse moves 
    for i in np.arange(nmoves-1)+nmoves:
        if np.sign(moves[i,1]) != np.sign(moves[i+1,1]):
            precomp = np.array([[0., -delta], [0., delta]])*np.sign(moves[i+1,1])
            precomp[0] += moves[i]
            hystcomp_moves.extend([precomp[0], precomp[1]])
            flags.extend([False,True])
        else:
            hystcomp_moves.append(moves[i])
            flags.append(True)
    
    hystcomp_moves.append(moves[-1])
    flags.append(True)
            
    return np.array(hystcomp_moves), np.array(flags)

    
def compute_autocoll_moves(targets_coord, hystcomp = False, CEg = 1.02928, Elg = 1.00379 ):

    ntargets = len(targets_coord)
    theta_EL_CE = np.zeros((ntargets,2))
    moves = np.zeros((ntargets,2))
    
    print("\ntargets local coords in siderostat local frame:")
    for i, c in enumerate(targets_coord):    
        #frame move
        theta_EL_CE[i,0] =  c.lon.deg/2./Elg
        theta_EL_CE[i,1] = -c.lat.deg/2./CEg
        print("EL: {:.1f} arcsec ; CE: {:.1f} arcsec".format(theta_EL_CE[i,0]*3600, theta_EL_CE[i,1]*3600))
        
    # compute relative moves
    moves[0] = theta_EL_CE[0] # 1st move from guider center
    moves[1:] = theta_EL_CE[1:] - theta_EL_CE[:-1]
    moves *= 3600
    
    if hystcomp:
        moves, flags = add_hystcomp(moves)
    else:
        flags = np.ones(ntargets, dtype=bool)
    
    print("\nsiderostat moves sequence: ")
    # print(moves.shape, flags.shape)
    # for m,f in zip(moves, flags):
    #     print("EL: {:.1f} arcsec ; CE {:.1f} arcsec ; image: {}".format(m[0], m[1], f))
    
    return moves, flags
    

    
class Guider2UV(object):
    '''
    Determine FOV center position on Guider, with data from perf 3 test 
    assuming:
        - Science Mask distortion & magnification used to build the mask is valid
        - The Guider to sky mapping (including distortion) measured with skytest 
        globular cluster is valid, this mapping makes the transformation from 
            + local sky coordinates (centered on guider center 640 ,540 with 
            angles on horizontal & vertical axis)
            + position x,y (horizontal,vertical) in pixels
        
        - For a Given mask we have both (from perf 3 test data):
            * a guider image of the 4 RC6 stars (to know the exact pattern)
            * a guider image with 3 RC6 stars visible & corresponding to the 4th star
            (middle right on guider image) exctly out on the slit and visible in UV
            
        The method is the following:
            - from the 4 stars measured positions on guider compute the 4 stars
            local sky coord
            - from the shifted 3 stars measured positions on guider compute the
            local sky coord
            - find the translation (spherical rotation),  between the 2 local frames
            that minimise the angular distance between each of the 3 stars 
            - knowing the above translation we can determine the coordinates 
            of the 4th star in the local frame where the uv is visible. 
            - from the local coord of the slit, we can compute:
                - the slit position on guider
                - any other related position on science mask to guider pixels.
    '''
    
 
    def __init__(self, 
                 guider_astrometry_filename=None,
                 Field_center=coordinates.SkyCoord(0*u.deg, 0*u.deg), 
                 Field_rotation=0*u.deg,
                 mask_rotation=0.*u.deg,
                 FOVcenter_guider_coord=None,
                 Field_gamma=1.,  # platescale correction
                 guider_wcs=None,
                 filename=None):
        '''
        set up mask astrometry and Field astrometry
        '''
        if filename is None:
        
            # define Field Local Frame
            self.FieldP = LocalScienceMaskProjector(Field_center, Field_rotation, Field_gamma)
            self.mask_rotation = mask_rotation
            
            # define Guider Local Frame
            if guider_astrometry_filename is not None:
                hdu = fits.open(guider_astrometry_filename)
                h = hdu[0].header
                w  = wcs.WCS(h)
                #re move the cd (rotation), but keep mean magnification
                scale = proj_plane_pixel_scales(w)
                del h['CD1_1']
                del h['CD1_2']
                del h['CD2_2']
                del h['CD2_1']
                h['CDELT1'] = scale[0]
                h['CDELT2'] = scale[1]
                guider_wcs = wcs.WCS(h)        

            self.GuiderP = LocalGuiderMaskProjector(guider_wcs)
            
            if FOVcenter_guider_coord is not None:
                self.FOV_center_guider_coord = FOVcenter_guider_coord
                
        else:
            # restore 
            self.restore(filename)


    def __str__(self):
        str_ = '''
Guider2UV object:
    Local Field Projector: {}
    Guider Field Projector: {}
    mask_rotation: {}
    FOV center in guider: {}x{} pix
'''.format(self.FieldP, self.GuiderP, self.mask_rotation, 
        self.FOV_center_guider_pix[0], self.FOV_center_guider_pix[1])

        return str_

        
    def copy(self):
        return Guider2UV(guider_wcs=self.GuiderP.w.copy(), 
                     Field_center=self.FieldP.center.copy(), 
                     Field_rotation=self.FieldP.rotation.copy(),
                     mask_rotation=self.mask_rotation.copy(),
                     Field_gamma = self.FieldP.gamma,
                     FOVcenter_guider_coord = self.FOV_center_guider_coord.copy())

        
    def set_FOV_center(self, FourStarsGuider_pos, StarsGuider_pos, Slit_pos, Stars_id=[0,1,2], 
                UVStar_id=3, world=False):
        '''
        Determine the FOV center in guider from the 4 star RC6 pattern and a 3 or less star 
        positions corresponding to the 4th star visible in UV, knowing te slit
        '''
        if world:
            slit_local_coord = self.FieldP.world2local(Slit_pos)
        else:
            slit_local_coord = self.FieldP.pix2local(Slit_pos)  
          
        self.pattern_coord = self.GuiderP.pix2local(FourStarsGuider_pos)
  
        UVstar_coord  = self.UVstar_guider_coord(StarsGuider_pos, Stars_id, UVStar_id)
        
        FOV_center_guider_lon = UVstar_coord.lon - slit_local_coord.lon
        FOV_center_guider_lat = UVstar_coord.lat - slit_local_coord.lat
        self.FOV_center_guider_coord = coordinates.SkyCoord(FOV_center_guider_lon, FOV_center_guider_lat, frame=self.GuiderP.localframe)
        print("FOV center angular position in guider", self.FOV_center_guider_coord)   
        print("FOV center pixel position in guider", self.FOV_center_guider_pix)
    

    @property
    def FOV_center_guider_pix(self):
        return self.GuiderP.local2pix(self.FOV_center_guider_coord)
    

    def UVstar_guider_coord(self, StarsGuider_pos, Stars_id=[0,1,2], 
                            UVStar_id=3, FourStarsGuider_pos=None):

        if FourStarsGuider_pos is not None:
            this_pattern_coord = self.GuiderP.pix2local(FourStarsGuider_pos)
        else:
            this_pattern_coord = self.pattern_coord

        Stars_coord = self.GuiderP.pix2local(StarsGuider_pos)

        UVstar_coord = register_UV_star_simple(this_pattern_coord, Stars_coord,
                                               Stars_id, UVStar_id)
        
        return coordinates.SkyCoord(UVstar_coord[0]*u.deg, 
                                    UVstar_coord[1]*u.deg, 
                                    frame=self.GuiderP.localframe)
    
    
    def FieldLocal_to_guider(self, Fcoord_local, angle=False):
        '''
        Determine objects position (local angles or pixels) on guider from their 
        local coord in the Field 
        '''
        
        # rotate
        rot = self.mask_rotation .to('rad').value  
        mrot = np.array([[1, -rot],
                         [rot, 1]])
        lon = np.array(Fcoord_local.lon.deg)
        lat = np.array(Fcoord_local.lat.deg)

        F_local = mrot.dot(np.vstack((lon, lat)))
        
        # recenter to guider
        local_lon = self.FOV_center_guider_coord.lon + F_local[0]*u.deg
        local_lat = self.FOV_center_guider_coord.lat + F_local[1]*u.deg
        Fcoord_local = coordinates.SkyCoord(local_lon, local_lat, frame=self.GuiderP.localframe)
        if angle:
            return Fcoord_local
        else:
            return self.GuiderP.local2pix(Fcoord_local)


    def guider_to_FieldLocal(self, Gcoord_local, angle=False):
        '''
        Determine object local angles on Science Mask from their 
        coordinates (either local angles or pixdels) on guider
        '''

        if not angle:
            Gcoord_local = self.GuiderP.pix2local(Gcoord_local)
            
        # recenter to FOV center    
        F_local_lon = Gcoord_local.lon - self.FOV_center_guider_coord.lon
        F_local_lat = Gcoord_local.lat - self.FOV_center_guider_coord.lat
        
        # unrotate
        rot = self.mask_rotation .to('rad').value  
        mrot = np.array([[1, rot],
                         [-rot, 1]])
        F_local = mrot.dot(np.array([F_local_lon.deg, F_local_lat.deg]))
        return coordinates.SkyCoord(F_local[0]*u.deg, F_local[1]*u.deg, frame=self.FieldP.localframe)
    

    def SienceMask2guider(self, coords, world=False, angle=False):
        '''
        Determine objects position (local angles or pixels) on guider from their 
        coordinates (either ra/dec or mask position in mm) in the Field
        '''
        if world:
            local_coords = self.FieldP.world2local(coords)
        else:
            local_coords = self.FieldP.pix2local(coords)    
            
        return self.FieldLocal_to_guider(local_coords, angle=angle)


    def guider2ScienceMask(self, coords, angle=False, world=False):
        '''
        Determine object position (ra/dec or mm) on Science Mask from their 
        coordinates (either local angles or pixels) on guider
        '''
        Fcoord = self.guider_to_FieldLocal(coords, angle=angle)
        if world:
            return self.FieldP.local2world(Fcoord)
        else:
            return self.FieldP.local2pix(Fcoord)


    def set_detector_mapping(self, mapping, offsets = [1., 1.]):

        self.mask_det_map = mapping
        #dxm_dyd = -0.013
        #dym_dxd = -0.015
        #m = np.array([[0., dxm_dyd],[dym_dxd ,0.]])

        self.direct_map = lambda w, x, y: self.mask_det_map.map(w, x, y).T + offsets
        self.inv_map = lambda w, x, y: self.mask_det_map.inv_map(w, x - offsets[0], y - offsets[1])
    
    
    def detector2guider(self, coord_det, wave= 0.20619, angle=True ):
                        
        pos = self.inv_map(wave, coord_det[:, 0], coord_det[:, 1]).T
        
        local_coords = [self.SienceMask2guider(p, angle=angle) for p in pos]
    
        return local_coords


    def compute_autocoll_moves_slits(self, slits, slit_table, hystcomp = False, CEg = 1.02928, Elg = 1.00379):
               
        slits = np.array(slits)
        #nslits = slits.size
        slit_coords = []
        st = slit_table
        if "xmm" in st.colnames:
            x,y = "xmm", "ymm"
        else:
            x,y = "x_mm", "y_mm"


        for s in slits:
            slit_pos = np.array([st[st['Internal-count']==s][x][0], st[st['Internal-count']==s][y][0]])
            print("slit position in mm on mask:", slit_pos)
            slit_coords.append(self.SienceMask2guider(slit_pos, angle=True))
        
    
        moves, flags = compute_autocoll_moves(slit_coords, hystcomp, CEg = CEg, Elg = Elg)
        slit_coords = slit_coords + slit_coords[::-1] # revert 
        
        return moves, flags, slit_coords
    
        
    def compute_autocoll_move_stars(self, stars, star_table, hystcomp = False, CEg = 1.02928, Elg = 1.00379):
        
        stars = np.array(stars)
        #nstars = stars.size
        star_coords = []
        st = star_table
        all_star_coords = coordinates.SkyCoord(st['RA']*u.deg, st['DEC']*u.deg)

        for s in stars:
            star_pos_radec = all_star_coords[st['Internal count']==s]
            # print("star position Ra/Dec: ", star_pos_radec)
            star_coords.append(self.SienceMask2guider(star_pos_radec, world=True, angle=True))
    
        moves,flags = compute_autocoll_moves(star_coords, hystcomp, CEg =CEg, Elg = Elg)
        
        return moves, flags, star_coords


    def pattern_in_slit(self, starid , Slit_pos, world=False, FourStarsGuider_pos=None):
        '''
        Show the 4 stars RC6 pattern on guider for one of them set on a slit
        
        note: the 4 stars RC6 pattern should be the same (rotation included) 
        than the one used to determine FOV 
        '''
        if world:
            slit_local_coord = self.FieldP.world2local(Slit_pos)
        else:
            slit_local_coord = self.FieldP.pix2local(Slit_pos)    
        
        if FourStarsGuider_pos is not None:
            this_pattern_coord = self.GuiderP.pix2local(FourStarsGuider_pos)
        else:
            this_pattern_coord = self.pattern_coord
            

        pattern_origin_Field_coord_lon = slit_local_coord.lon - this_pattern_coord[starid].lon
        pattern_origin_Field_coord_lat = slit_local_coord.lat - this_pattern_coord[starid].lat
        
        pattern_Field_coord_lon = pattern_origin_Field_coord_lon + this_pattern_coord.lon
        pattern_Field_coord_lat = pattern_origin_Field_coord_lat + this_pattern_coord.lat
        pattern_Field_coord = coordinates.SkyCoord(pattern_Field_coord_lon, 
                                                   pattern_Field_coord_lat,
                                                   frame=self.FieldP.localframe)
 
        return self.FieldLocal_to_guider(pattern_Field_coord)


    def save(self, filename='Guider2UV.pkl'):
        '''
        Save the FOV registration, and astrometry configuration in a pickle file
        '''
        print("Dumping to " + filename)
        
        data = {'Fcenter': self.FieldP.center,
                'Frot': self.FieldP.rotation,
                'Gwcs': self.GuiderP.w,
                }
        
        if self.FieldP.gamma != 1.:
            data['Fgamma'] = self.FieldP.gamma
        
        if self.mask_rotation != 0*u.deg:
            data['Fmaskrot'] = self.mask_rotation 

        try:
            data['FOVcenter'] = (self.FOV_center_guider_coord.lon, self.FOV_center_guider_coord.lat)
        except AttributeError:
            pass
        
        try:
            # save pattern_coord
            data['pattern'] = (self.pattern_coord.lon, self.pattern_coord.lat)
        except AttributeError:
            pass
        
        with  open(filename,'wb') as file_handle:            
            pickle.dump(data, file_handle)
        
        
    def restore(self, filename='Guider2UV.pkl'):
        '''
        Restore a saved configuration (in pickle file)
        '''    
        with open(filename,'rb') as file_handle:
            data = pickle.load(file_handle)
        
        self.filename = filename
        
        try:
            Fgamma = data['Fgamma']
        except  KeyError:
            Fgamma = 1.
            
        try:
            self.mask_rotation = data['Fmaskrot']
        except KeyError:
            self.mask_rotation = 0*u.deg
                
        self.FieldP = LocalScienceMaskProjector(data['Fcenter'], data['Frot'], Fgamma)
        self.GuiderP = LocalGuiderMaskProjector(data['Gwcs'])
        try:
            FOV_center_guider_coord = data['FOVcenter']
        except KeyError:
            print("No FOV center")
        else:
            self.FOV_center_guider_coord = coordinates.SkyCoord(FOV_center_guider_coord[0],
                                                         FOV_center_guider_coord[1], 
                                                         frame=self.GuiderP.localframe)
            print("FOV center angular position in guider", self.FOV_center_guider_coord)

            print("FOV center pixel position in guider", self.FOV_center_guider_pix)
            
        try:
            #restore pattern
            self.pattern_coord = coordinates.SkyCoord(data['pattern'][0], 
                                                      data['pattern'][1],
                                                      frame = self.GuiderP.localframe)
        except KeyError:
                pass


    def  update_model(self, coord, coord_obs, gamma=True, ytilt=False, weight=None,
                      inplace=False, plot=False, labels=None):
        
        sol, residuals = fit_model(coord, coord_obs, gamma, ytilt, weight)
        
        if plot:
            plot_fit(coord, coord_obs, residuals, labels=labels, sol=sol)
            
        if inplace: 
            G2UVcor = self
        else:
            G2UVcor = self.copy()

        dgamma, theta_rad, thetay_rad, deltax, deltay = sol
        gamma = 1. + dgamma
        deltax = deltax * u.deg
        deltay = deltay * u.deg
        theta = theta_rad*180/np.pi*60*u.arcmin

        new_offset = [-theta_rad*G2UVcor.FOV_center_guider_coord.lat + deltax,
                       theta_rad*G2UVcor.FOV_center_guider_coord.lon + deltay]
                        
        newFOVcenter = coordinates.SkyCoord(G2UVcor.FOV_center_guider_coord.lon*gamma + new_offset[0],
                                            G2UVcor.FOV_center_guider_coord.lat*gamma + new_offset[1],
                                           frame = G2UVcor.GuiderP.localframe)
                
        G2UVcor.mask_rotation = G2UVcor.mask_rotation + theta
        G2UVcor.FOV_center_guider_coord = newFOVcenter
        
        G2UVcor.FieldP.gamma *= gamma
        
        return G2UVcor, residuals


    def correct_guider_coord(self, coord, gamma, theta, dx, dy):
        
        if not isinstance(coord, list):
            coord = [coord]
        
        n = len(coord)
        coord_arr = np.zeros((n, 2))
        for i,c in enumerate(coord):
            coord_arr[i] = np.array([ c.lon.deg, c.lat.deg]).T
        
        row_x = np.hstack((coord_arr*[1.,-1], np.ones((n,1)), np.zeros((n,1)) )) # xn - x =  x dgamma - y theta + dx  
        row_y = np.hstack((coord_arr[:,::-1], np.zeros((n,1)), np.ones((n,1)) )) # yn - y =  y dgamma + x theta      + dy 
        mat = np.vstack((row_x, row_y))
          
        sol = np.array([gamma, theta, dx, dy])
        new = mat.dot(sol).reshape((2,n)).T
        new += coord_arr
        
        new_coord = []
        for i in range(n):
            new_coord.append(coordinates.SkyCoord(new[i,0]*u.deg ,new[i,1]*u.deg, frame=self.GuiderP.localframe))
        
        if n==1:
            new_coord = new_coord[0]
        
        return new_coord, new
 
        
if __name__ == "__main__":
    
    from astropy.table import Table
    from matplotlib import pyplot as plt
    
    path = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/perf3.0-09-01/sky_00_astrometry/"
#    path = "/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/perf3.0-09-01/sky_00_astrometry/"

#    F3_astrometry_filename = path + 'stack475788.new'
#    F2_astrometry_filename = path + 'stack486607.new'
    F1_astrometry_filename = path + 'stack500513.new'
#    F4_astrometry_filename = path + 'stack502282.new'

#    Field_center = coordinates.SkyCoord(352.3424*u.deg, 0.21245*u.deg)
    Field_center = coordinates.SkyCoord(32.19*u.deg, -5.688*u.deg)

#    G2UV = Guider2UV(F3_astrometry_filename, Field_center )
    G2UV = Guider2UV(F1_astrometry_filename, Field_center )
    
    
    #################
    ################# SC-GUI test
    newpath = "/home/dvibert/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/"
#    newpath = "/Users/vincent/ownCloud/FIREBALL/Tests-at-FortSumner/SC_GUI01/"
#    new_F1_star_pattern_filename = newpath + 'F1_stack7924678.fits_table.fits' #F1
#    new_F1_star_pattern_tab = Table.read(new_F1_star_pattern_filename)  
#    new_F1_star_pattern_tab.sort('xcentroid')
#    new_F1_star_pattern_pos = np.array([new_F1_star_pattern_tab['xcentroid'], new_F1_star_pattern_tab['ycentroid']]).T
     
    new_F1_star_pattern_filename = newpath + 'F1_stack7971979.fits_table.fits'
    new_F1_star_pattern_tab = Table.read(new_F1_star_pattern_filename)  
    new_F1_star_pattern_tab.sort('xcentroid')
    new_F1_star_pattern_pos = np.array([new_F1_star_pattern_tab['xcentroid'], new_F1_star_pattern_tab['ycentroid']]).T
    

    starid = 3
    slit_pos = np.array([0.5525186, -4.7601449]) #30
    slit_pos_ang = coordinates.SkyCoord(32.2022*u.deg, -5.8006*u.deg)
    slit_pos = np.array([0, 0]) #30
    slit_pos_ang = Field_center
    
    
    
    
    #print(G2UV.pattern_in_slit(starid, slit_pos, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos))
#    prediction = G2UV.pattern_in_slit(starid, slit_pos, world=False, FourStarsGuider_pos=new_F1_star_pattern_pos)
    
    plt.figure()
    plt.plot(new_F1_star_pattern_pos[:,0], -new_F1_star_pattern_pos[:,1], '+')
    plt.xlim([0,1400])
    plt.ylim([0,-1100])
    for i in  range(4):
        plt.text(new_F1_star_pattern_pos[i,0], -new_F1_star_pattern_pos[i,1], str(i))
    
#    plt.plot(prediction[0], prediction[1], 'or')
    
    
    ### register with new data
    #star_UV_filename = path + 'starsinguiderwhenuvsource/stack1399407_table.fits' #F?
    star_UV_filename = newpath + '/Guider/img/stack8074052.fits_table.fits' #F1
    star_UV_tab = Table.read(star_UV_filename)  
    star_UV_tab.sort('xcentroid')
    star_UV_pos = np.array([star_UV_tab['xcentroid'], star_UV_tab['ycentroid']]).T
    plt.plot(star_UV_pos[:,0], -star_UV_pos[:,1], 'x')
    for i in  range(3):
        plt.text(star_UV_pos[i,0], -star_UV_pos[i,1], str(i))
    
    #slit_pos = coordinates.SkyCoord(352.445678710937*u.deg, 0.117206886410713*u.deg)
    #G2UV.set_FOV_center(new_F1_star_pattern_pos, star_UV_pos, slit_pos, UVStar_id=3, world=False)
    G2UV.set_FOV_center(new_F1_star_pattern_pos, star_UV_pos, slit_pos_ang, UVStar_id=3, world=True)
    plt.plot(G2UV.FOV_center_guider_pix[0], -G2UV.FOV_center_guider_pix[1],'o')
    plt.grid()
#
#    # check prediction with neighboring slit #25
    new_slit_pos = np.array([-0.905, -3.2085])
    new_slit_pos_ang = coordinates.SkyCoord(32.1685*u.deg, -5.7639*u.deg)
    


        

    starid = 3
    print(G2UV.pattern_in_slit(starid, new_slit_pos, world=False))
    pred = G2UV.pattern_in_slit(starid, new_slit_pos_ang, world=True)
    print(pred) 
    #[array([  350.27358904,   686.81014691,   742.7547067 ,  1068.75873879]),
    #array([  699.59574877,   311.1393818 ,  1024.16645599,   640.88014618])]
    # should have been LFS 329.3  750.5


    #new slit (3rd)
    #guider pos 430.2 283.1 #36
    new_slit_pos_ang = coordinates.SkyCoord(32.2641*u.deg, -5.7229*u.deg)
    starid = 3    
    print(G2UV.pattern_in_slit(starid, new_slit_pos, world=False))
    pred = G2UV.pattern_in_slit(starid, new_slit_pos_ang, world=True)
    print(pred) 
        
    
    pred = G2UV.pattern_in_slit(starid, new_slit_pos, world=False)
    plt.figure()
    plt.plot(star_UV_pos[:,0],  star_UV_pos[:,1], '+')
    plt.xlim([0,1300])
    plt.ylim([0,1100])

    plt.plot(pred[0], pred[1], 'or')
    ### prediction is goind wrong direction in y (vertical guider)
    
 
    #----
#    UVstar_guider_coord = test.register_UV_star_lsqr(star_pattern_local_coord, star_UV_local_coord)
#
#    plt.figure()
#    plt.plot(star_pattern_local_coord.lon.deg, star_pattern_local_coord.lat.deg,'+')
#    plt.plot(star_UV_local_coord.lon.deg, star_UV_local_coord.lat.deg, 'ob')
#    plt.plot(UVstar_guider_coord[0], UVstar_guider_coord[1],'or')
#
#    UVstar_guider_coord2 = test.register_UV_star_simple(star_pattern_local_coord, star_UV_local_coord)
#    
#    plt.figure()
#    plt.plot(star_pattern_local_coord.lon.deg, star_pattern_local_coord.lat.deg,'+')
#    plt.plot(star_UV_local_coord.lon.deg, star_UV_local_coord.lat.deg, 'ob')
#    plt.plot(UVstar_guider_coord2[0], UVstar_guider_coord2[1],'or')
#
#    #---
#    UVstar_radec = coordinates.SkyCoord()
#    UVstar_local_coord = FieldP.world2local(UVstar_radec)
    
