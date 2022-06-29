#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:15:55 2019

@author: Vincent
"""    
from astropy import units as u    
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
from scipy import special as sp
from decimal import Decimal
from astropy.convolution import Gaussian2DKernel


def z_to_proper(z, H0=67.77, OmegaM=0.315):
   """
   compute proper distance in Mpc from redshift in a flat universe
   using Adachi & Kasai (2011) Luminosity distance approximation
   """
   return z_to_comoving(z, H0, OmegaM)/(1.+z)




def z_to_comoving(z, H0=67.77, OmegaM=0.315): # H0 in km/s/MPc
   """
   compute comoving distance in Mpc from redshift in a flat universe
   using exact analytic formula
   """
   #Define useful cosmological parameters
   c = 299792.458 #Speed of light in km/s
   Dh = c/H0
   m = 1./6.
   lambdazero=1.-OmegaM
   
   #Calculate distance
   x = lambdazero*((1.+z)**(-3))/(OmegaM+(lambdazero*((1.+z)**(-3))))
   Dc = Dh*2.*(m/np.sqrt(OmegaM))*((OmegaM/lambdazero)**m)*sp.beta(m,0.5-m)*(sp.betainc(m,0.5-m,lambdazero)-sp.betainc(m,0.5-m,x))
#    """
#    compute comoving distance in Mpc from redshift in a flat universe
#    using Adachi & Kasai (2011) Luminosity distance approximation
#    """
#    c = 299792.458 # speed of light in km/s
#    x = ((1-OmegaM)/OmegaM)/(1+z)**3
#    x0 = ((1-OmegaM)/OmegaM)
#    phi_x = 1 + 1.320*x + 0.4415*x**2 + 0.02656*x**3
#    phi_x /= 1 + 1.392*x + 0.5121*x**2 + 0.03944*x**3
#    phi_x0 = 1 + 1.320*x0 + 0.4415*x0**2 + 0.02656*x0**3
#    phi_x0 /= 1 + 1.392*x0 + 0.5121*x0**2 + 0.03944*x0**3
#    Dc = 2*c/H0/np.sqrt(OmegaM)*(phi_x0 - 1/np.sqrt(1+z)*phi_x)
   return Dc



#table = Table.read('/Users/Vincent/Nextcloud/Work/Target_selection/2020Palestine/eBOOS/DR12Q_FBFOV.fits')
#table['z_Mpc'] = z_to_comoving(table['Z_VI'])
#table.write('/Users/Vincent/Nextcloud/Work/Target_selection/2020Palestine/eBOOS/DR12Q_FBFOV.fits',overwrite=True)

def z_to_luminosity(z, H0=67.3, OmegaM=0.315):  # H0 in km/s/MPc
    """
    compute luminosity distance in Mpc from redshift in a flat universe
    using Adachi & Kasai (2011) Luminosity distance approximation
    """    
    return z_to_comoving(z, H0, OmegaM) * (1 + z)



def Euclidian_distance_qso(ra1,dec1,z1,ra2,dec2,z2):
    """Nous interesse que jusqua 5Mps a utiliser avec des proper mpc!!!"""
    c1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree, distance=z_to_proper(z1)*u.mpc)
    c2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree, distance=z_to_proper(z2)*u.mpc)
    sep = c1.separation_3d(c2)
    return sep ##(in Mpc)!!!

def Impact_parameter(ra1,dec1,z1,ra2,dec2,z2):
    """Nous interesse que jusqua 60kpc a utiliser avec des proper mpc!!!"""
    gal = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree, distance=z_to_proper(z1)*u.mpc)
    qso = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree, distance=z_to_proper(z1)*u.mpc)
    sep = gal.separation(qso).value
    return sep ##(in Mpc)!!!



def Angular_distance_qso(ra1,dec1,z1,ra2,dec2,z2):
    c1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree, distance=z_to_proper(z1)*u.mpc)
    c2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree, distance=z_to_proper(z2)*u.mpc)
    sep = c1.separation(c2)
    return sep ##degreee



def Lum2Flux(lum, z):
    from astropy.cosmology import WMAP9 as cosmo
    rcm = cosmo.luminosity_distance(z).to(u.cm)
    flux = lum / (4*np.pi*rcm**2)
    return np.log10(flux.data)

def Lum2MagAbs(lum):
    rcm = 10*u.pc.to(u.cm)
    flux = lum / (4*np.pi*rcm**2)
    return np.log10(flux)


def Lum2MagAbs2(lum):
    """in erg/s
    """
    ergs2watts = 1e-7#(1*u.erg/u.s).to(u.watt)
    lum *= ergs2watts
    Ms = 4.85
    Ls = 3.0128e28#watts
    M = Ms - np.log10(lum/Ls)/0.4
    return M

def Lum2MagAbs3(lum):
    """in erg/s
    """
    lum /= 3.0128e28
    Ms = 51.605
    M = Ms - np.log10(lum)/0.4
    return M


def Lum2MagAbs4(lum):
    """in erg/s
    """
    return 4.85 - 2.5 * np.log10(lum/3.0128e28)


def MagAbs2Lum(Mag):
    return 10**(-0.4*(Mag - 51.605))

#def Flux2MagAbs(Flux):
    

def mag_rel_to_abs(mag_relative,z=None, pc=None):
    if pc is None:
        mag_abs = mag_relative -5 *  np.log10(z_to_luminosity(z) * 1e6) +5 +2.5*np.log10(1+z)
    if z is None:
        mag_abs = mag_relative -5 *  np.log10(pc) +5 
    return mag_abs


def mag_abs_to_app(mag_abs,z=None, pc=None):
    if pc is None:
        mag_apparente = mag_abs + 5 * np.log10(z_to_luminosity(z) * 1e6) -5 - 2.5 * np.log10(1+z)
    if z is None:
        mag_apparente = mag_abs + 5 * np.log10(pc) -5      
    return mag_apparente

def Erg_per_sec_2_power(Flux=6.2e-17, redshift=0.7):
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    #cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    a = cosmo.luminosity_distance(redshift)  
    dl_cm = a.to(u.cm)
    return Flux*(4*np.pi*dl_cm**2).value

def Photons2ergs(_lambda=2050):
    # c=3e8, h=6.62e-34
    eV = 1.2398/(_lambda*1e-4)#In microns
    ergs = eV*1.6022e-12
    # Verification : 1.2398/(12400*1e-4)=1
    return ergs
#Photons2ergs()

def DetectionLimit(ADULimit=50, exposure=50, ConvGain=0.53, emGain=1370, PhotonsPerCm2=208 ):
    #Needs to be beforehand convolved by the size of the slit to not be by Angstrom but by slit = element resolution
    #208 ce sont des (phel/s/A )  /   (photons/s/cm2/A)   =  phel .  cm2 / photon(hors atmosphere) 
    lim = ADULimit/exposure/ConvGain/emGain/PhotonsPerCm2
    return lim
    #Erg_per_sec_2_power()
    
    

#sur une fente sur un objet
ADULimit = 5*500/Gaussian2DKernel(x_stddev=7/2.35).array.max()#150*sqrt(9)*sqrt(8)*sqrt(3)#5*std
print('FIREBall detection limit is %0.2E Photons/cmˆ2/sec'%(Decimal(DetectionLimit(ADULimit))))
print('FIREBall detection limit is %0.2E ergs/sec/cm^2'%(Decimal(Photons2ergs(2050)*DetectionLimit(ADULimit))))
print('FIREBall detection limit is %0.2E erg/sec'%(Decimal(Erg_per_sec_2_power(Photons2ergs(2050)*DetectionLimit(ADULimit),redshift=0.7))))
####a diviser par les elements de resolution
print('FIREBall detection limit is %0.2E erg/sec for a resol'%(Decimal(Erg_per_sec_2_power(Photons2ergs(2050)*DetectionLimit(ADULimit),redshift=0.7))))




#Bruit sur le flux
ADULimit = 5*8000#150*sqrt(9)*sqrt(8)*sqrt(3)#5*std
print('FIREBall detection limit is %0.2E Photons/cmˆ2/sec'%(Decimal(DetectionLimit(ADULimit))))
print('FIREBall detection limit is %0.2E ergs/sec/cm^2'%(Decimal(Photons2ergs(2050)*DetectionLimit(ADULimit))))
print('FIREBall detection limit is %0.2E erg/sec'%(Decimal(Erg_per_sec_2_power(Photons2ergs(2050)*DetectionLimit(ADULimit),redshift=0.7))))
####a diviser par les elements de resolution
print('FIREBall detection limit is %0.2E erg/sec for a resol'%(Decimal(Erg_per_sec_2_power(Photons2ergs(2050)*DetectionLimit(ADULimit),redshift=0.7))))

###


#sur un stack de 8 lignes, 16 objets, 100 images 

ADULimit = 5*500 / (np.sqrt(16)*np.sqrt(77))/Gaussian2DKernel(x_stddev=6/2.35).array.max()#150*sqrt(9)*sqrt(8)*sqrt(3)#5*std
print('FIREBall detection limit is %0.2E Photons/cmˆ2/sec'%(Decimal(DetectionLimit(ADULimit))))
print('FIREBall detection limit is %0.2E ergs/sec/cm^2'%(Decimal(Photons2ergs(2050)*DetectionLimit(ADULimit))))
print('FIREBall detection limit is %0.2E erg/sec'%(Decimal(Erg_per_sec_2_power(Photons2ergs(2050)*DetectionLimit(ADULimit),redshift=0.7))))

#
#fit une gaussienne sur une image de std=500
#le coefficient de la gaussienne donne le flux total moyen (cf vstack)
#Johnson papier 
#Je prend mon image, soit je fit par une gaussienne partout dans le bruit et le prend la std de la mesure du flux
#siot je smooth par une gaussienne et je regarde direct le bruit





#Or on a un background de 
#(1.34e-41*u.photometric.Bol).to(u.photometric.AB)

if __name__ == '__main__':
    lum = 10e42
    rpc = 10*u.pc.to(u.cm)
    Lum2MagAbs4(37*3.0128e28)
    Lum2MagAbs4(3.0128e28)
    Lum2MagAbs4(110*3.0128e28)
    
    Lum2MagAbs4(1e42)
    
    mag_abs_to_app(-28,z=0.7)




    Lum2MagAbs2(37*3.0128e28*1e6)
    Lum2MagAbs3(3.0128e28)
    
    
    Lum2MagAbs2(10e42)
    
    MagAbs2Lum(Lum2MagAbs3(3.0128e28*1e7))
    
    Lum2MagAbs2(3.826e33)
    
    Lum2MagAbs3(3.826e33)
    
    mag_abs_to_app(-13.95,z=0.7)
    mag_abs_to_app(4.85,pc=0.7)
    
    mag_abs_to_app(4.83,pc=4.83e-6)
    
    mag_abs_to_app(1.4,pc=2.68)
    mag_abs_to_app(-0.2,pc=10.66)
    
    
    mag_abs_to_app(-0.2,pc=10.66)
    
    Lum2Flux(10e42, 4)
    Lum2MagAbs(10e42)
    
    
    
#ok: lum 2 flux, mnagabs2 mag app  , mag app 2 mag abs
# not ok: Lum 2 mag abs (only for sun?, does not seem ok) or flux to mag abs!!!