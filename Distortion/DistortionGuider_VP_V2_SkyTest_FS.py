from __future__ import division, print_function

get_ipython().magic(u'matplotlib notebook')

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
#from polynomialD import polyfit2d
from astropy.io import fits
from PolyFit2d import PolyFit2d

def sky_Coord(ra,dec,ra_c,dec_c):
    """Given a line sight, the function converts right ascenssion and declination
    into alpha_x alpha_y the small angles to center line of sight (ra_c,dec_c)"""
    numx = -np.cos(dec*np.pi/180)*np.sin((ra - ra_c)*np.pi/180)
    numy = -np.cos((ra - ra_c)*np.pi/180) * np.cos(dec*np.pi/180) * np.sin(dec_c*np.pi/180)  + np.sin(dec*np.pi/180) * np.cos(dec_c*np.pi/180)
    denom = np.cos((ra - ra_c)*np.pi/180) * np.cos(dec*np.pi/180) * np.cos(dec_c*np.pi/180)  + np.sin(dec*np.pi/180) * np.sin(dec_c*np.pi/180)
    alpha_x = - np.arctan(numx/denom)*180/np.pi
    alpha_y = np.arctan(numy/denom)*180/np.pi
    return np.array([alpha_x,alpha_y])


def find_center(stars_alt, stars_az,psf_x, psf_y,alt_c,az_c,delta = 0.01, pas=100, plot=True):
    FCx = 1200.
    FCy = 1080/2.
    l1 = []
    l2 = []
    l3 = []
    l4 = []   
    for FCalt in np.linspace(alt_c - delta,alt_c + delta,pas):
        for FCaz in np.linspace(az_c - delta,az_c + delta,pas):
            y = sky_Coord(stars_az.ravel(), stars_alt.ravel(), FCaz, FCalt)[0] * 3600 
            x = sky_Coord(stars_az.ravel(), stars_alt.ravel(), FCaz, FCalt)[1] * 3600
            valx = psf_x.ravel() - FCx
            valy = psf_y.ravel() - FCy
            cx1 = PolyFit2d(x,y,valx,[1,1])
            cy1 = PolyFit2d(x,y,valy,[1,1])
            l1.append(FCalt)
            l2.append(FCaz)
            l3.append(abs(cx1.jacobian(0,0)[0] - cy1.jacobian(0,0)[1]))
            l4.append(cx1.coeffs[0,0]**2 + cy1.coeffs[0,0]**2)
        

    plt.figure()
    plt.imshow(np.array(l4).reshape(100,100),cmap = 'hot') 
    plt.colorbar()
    #print('le minimum de decentrament x est atteint en FCx,FCy = ',l1[abs(np.array(l3)).argmin()],l2[abs(np.array(l3)).argmin()])
    #print('Il vaut alors',l3[abs(np.array(l3)).argmin()],l4[abs(np.array(l3)).argmin()])
    
    #print('le minimum de decentrament y est atteint en FCx,FCy = ',l1[abs(np.array(l4)).argmin()],l2[abs(np.array(l4)).argmin()])
    #print('Il vaut alors ',l3[abs(np.array(l4)).argmin()],l4[abs(np.array(l4)).argmin()])
    
    #FCx = l1[(abs(np.array(l3)) + abs(np.array(l4))).argmin()]
    #FCy = l2[(abs(np.array(l3)) + abs(np.array(l4))).argmin()]
    
    
    FCalt2 = l1[(abs(np.array(l4))).argmin()]
    FCaz2 = l2[(abs(np.array(l4))).argmin()]
    y = sky_Coord(stars_az.ravel(), stars_alt.ravel(), FCaz2, FCalt2)[0] * 3600 
    x = sky_Coord(stars_az.ravel(), stars_alt.ravel(), FCaz2, FCalt2)[1] * 3600
    valx = psf_x.ravel() - FCx
    valy = psf_y.ravel() - FCy
    cx1 = PolyFit2d(x,y,valx,[1,1])
    cy1 = PolyFit2d(x,y,valy,[1,1])
    cx2 = PolyFit2d(x,y,valx,[2,2])
    cy2 = PolyFit2d(x,y,valy,[2,2])
    cx3 = PolyFit2d(x,y,valx,[3,3])
    cy3 = PolyFit2d(x,y,valy,[3,3])
    cx4 = PolyFit2d(x,y,valx,[4,4])
    cy4 = PolyFit2d(x,y,valy,[4,4])
    cx5 = PolyFit2d(x,y,valx,[5,5])
    cy5 = PolyFit2d(x,y,valy,[5,5])
    print('With this center, cx1[0,0] = ',cx1.coeffs[0,0])
    print('With this center, cy1[0,0] = ',cy1.coeffs[0,0])
    print('With this center, cx2[0,0] = ',cx2.coeffs[0,0])
    print('With this center, cy2[0,0] = ',cy2.coeffs[0,0])
    print('With this center, cx3[0,0] = ',cx3.coeffs[0,0])
    print('With this center, cy3[0,0] = ',cy3.coeffs[0,0])
    print('With this center, cx4[0,0] = ',cx4.coeffs[0,0])
    print('With this center, cy4[0,0] = ',cy4.coeffs[0,0])
    print('With this center, cx5[0,0] = ',cx5.coeffs[0,0])
    print('With this center, cy5[0,0] = ',cy5.coeffs[0,0])
    return  [FCalt2, FCaz2]




def compute_magnification(a1, minutes=45, secondes = 0, duration=95, deg = [1,1],alt_c = 53.721, az_c = 316.868 ):
    FCx = 1200.
    FCy = 1080/2.
    star1 = SkyCoord("15h19m34.0s +60d14m47.8s")
    star2 = SkyCoord("15h20m4.99s +60d22m49.9s")
    star3 = SkyCoord("15h20m8.73s +60d21m23.7s")
    star4 = SkyCoord("15h19m59.67s +60d20m23.2s")
    Toulouse = EarthLocation(lat=43.561048*u.deg, lon=1.483655*u.deg, height=200*u.m)
    delta_test = np.linspace(0, duration/60,14)*u.minute  #7.173062
    time = Time('2017-6-23 0:{}:{}'.format(minutes,secondes),'UTC') #- utcoffset
    stars11 = np.array([star1.transform_to(AltAz(obstime=time+delta_test,location=Toulouse)).az.value,star1.transform_to(AltAz(obstime=time+delta_test,location=Toulouse)).alt.value])
    stars21 = np.array([star2.transform_to(AltAz(obstime=time+delta_test,location=Toulouse)).az.value,star2.transform_to(AltAz(obstime=time+delta_test,location=Toulouse)).alt.value])
    stars31 = np.array([star3.transform_to(AltAz(obstime=time+delta_test,location=Toulouse)).az.value,star3.transform_to(AltAz(obstime=time+delta_test,location=Toulouse)).alt.value])
    stars41 = np.array([star4.transform_to(AltAz(obstime=time+delta_test,location=Toulouse)).az.value,star4.transform_to(AltAz(obstime=time+delta_test,location=Toulouse)).alt.value])
    starsaz1 = np.vstack((stars11[0],stars21[0],stars31[0],stars41[0]))
    starsalt1 = np.vstack((stars11[1],stars21[1],stars31[1],stars41[1]))
    stars_x1 = np.vstack((a1[a1['NUMBER']==1]['X_IMAGE'],a1[a1['NUMBER']==2]['X_IMAGE'],a1[a1['NUMBER']==3]['X_IMAGE'],a1[a1['NUMBER']==4]['X_IMAGE']))
    stars_y1 = np.vstack((a1[a1['NUMBER']==1]['Y_IMAGE'],a1[a1['NUMBER']==2]['Y_IMAGE'],a1[a1['NUMBER']==3]['Y_IMAGE'],a1[a1['NUMBER']==4]['Y_IMAGE']))
    FCalt1, FCaz1 = find_center(stars_alt = starsalt1, stars_az=starsaz1, psf_x = stars_x1, psf_y = stars_y1, alt_c = alt_c, az_c = az_c,delta=0.005 )
    y1 = sky_Coord(starsaz1.ravel(), starsalt1.ravel(), FCaz1, FCalt1)[0] * 3600 
    x1 = sky_Coord(starsaz1.ravel(), starsalt1.ravel(), FCaz1, FCalt1)[1] * 3600    
    valx1 = stars_x1.ravel() - FCx
    valy1 = stars_y1.ravel() - FCy    
    cx1 = PolyFit2d(x1,y1,valx1,deg)
    cy1 = PolyFit2d(x1,y1,valy1,deg)
    print (cx1.jacobian(0,0)[0],cy1.jacobian(0,0)[1])
    return cx1,cy1


def AnalyseDistortion(a,name = 'masque', scale = 0.2):

    table = Table(a.data)
    valx = table['xcentroid']
    valy = table['ycentroid']
    x = table['X']
    y = table['Y']
    w = np.square(table['flux']/table['FWHMx']) / table['FWHMx']

    cx1  = PolyFit2d(x,y,valx, [1,1])                                                                                                              
    cy1 = PolyFit2d(x,y,valy,[1,1])                                                                                                                    
    cx1w  = PolyFit2d(x,y,valx, [1,1],w=w)                                                                                                              
    cy1w = PolyFit2d(x,y,valy,[1,1],w=w) 
    cx2  = PolyFit2d(x,y,valx,[2,2])                                                                                                                    
    cy2 = PolyFit2d(x,y,valy,[2,2])                                                                                                               
    cx3  = PolyFit2d(x,y,valx,[3,3])                                                                                                              
    cy3 = PolyFit2d(x,y,valy,[3,3])                                                                                                                    
    cx4  = PolyFit2d(x,y,valx,[4,4])                                                                                                              
    cy4 = PolyFit2d(x,y,valy,[4,4])                                                                                                               

    cx5  = PolyFit2d(x,y,valx,[5,5])                                                                                                              
    cy5 = PolyFit2d(x,y,valy,[5,5]) 
    
#    plt.figure()
#    plt.plot(a['X'],a['Y'],'*', label = 'Theoric position')
#    plt.plot(a['xcentroid'],a['ycentroid'],'.',label='Star in guider camera')
#    plt.plot(cx2(x,y),cy2(x,y),'.',label='After fitting the distortion')    
#    plt.legend(loc='lower left', borderaxespad=0.)
#    plt.xlabel('Azimuth')
#    plt.ylabel('Altitude')
#    plt.title('Position of the stars in the guider camera')
#    plt.show()
    dist = np.sqrt(np.square( y-valy)+np.square(x-valx))
    plt.figure()
    qv1 = plt.quiver(valy, valx, y-valy, x-valx, scale_units='xy', angles='xy', color='b', label = 'Fit initial distortion',scale=0.02 )                                         
    plt.quiver(valy, valx, cy1(x,y)-valy, cx1(x,y)-valx, scale_units='xy', angles='xy', color='r', label = 'Fit final distortion',scale=0.02   )                                       
    plt.quiverkey(qv1, .2,-0.14, 1, 'scale: 13 microns', coordinates='axes', color='r')                                                                  
    plt.legend()
    plt.figtext(0.13,0.13,'mean = %.2f - %.2f - %.2f \nvar = %.2f - %.2f - %.2f \nmax = %.2f - %.2f - %.2f'%(np.mean(y-valy),np.mean(x-valx),np.mean(dist),np.var(y-valy),np.var(x-valx),np.var(dist),np.max(y-valy),np.max(x-valx),np.max(dist)))
    plt.title(name)    
    plt.xlabel('y')
    plt.ylabel('x')
    plt.grid(10)
    plt.show()

    dist = np.sqrt(np.square( y-valy)+np.square(x-valx))
    plt.figure()
    qv1 = plt.quiver(valy, valx, y-valy, x-valx, scale_units='xy', angles='xy', color='b', label = 'Fit initial distortion',scale=0.02 )                                         
    plt.quiver(valy, valx, cy1w(x,y)-valy, cx1w(x,y)-valx, scale_units='xy', angles='xy', color='r', label = 'Fit final distortion',scale=0.02   )                                       
    plt.quiverkey(qv1, .2,-0.14, 1, 'scale: 13 microns', coordinates='axes', color='r')                                                                  
    plt.legend()
    plt.figtext(0.13,0.13,'mean = %.2f - %.2f - %.2f \nvar = %.2f - %.2f - %.2f \nmax = %.2f - %.2f - %.2f'%(np.mean(y-valy),np.mean(x-valx),np.mean(dist),np.var(y-valy),np.var(x-valx),np.var(dist),np.max(y-valy),np.max(x-valx),np.max(dist)))
    plt.title(name + ' with wieght')    
    plt.xlabel('y')
    plt.ylabel('x')
    plt.grid(10)
    plt.show()

    
    
#    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')                                                                  
#    qv1 = ax1.quiver(valy, valx, valy-cy1(x,y), valx-cx1(x,y), scale_units='xy', angles='xy', scale=0.01, color='b')                                          
#    ax2.quiver(valy, valx, valy-cy2(x,y), valx-cx2(x,y), scale_units='xy', angles='xy', scale=0.01, color='b')                                                
#    ax3.quiver(valy, valx, valy-cy3(x,y), valx-cx3(x,y), scale_units='xy', angles='xy', scale=0.01, color='b')                                          
#    ax4.quiver(valy, valx, valy-cy4(x,y), valx-cx4(x,y), scale_units='xy', angles='xy', scale=0.01, color='b')                                          
#    ax1.quiverkey(qv1, 1,-1.12, 1, 'scale: 1 pixels', coordinates='axes', color='r')                                                                  
#    ax1.set_title('Deg 1')                                                                                                                        
#    ax2.set_title('Deg 2')                                                                                                                        
#    ax3.set_title('Deg 3')                                                                                                                        
#    ax4.set_title('Deg 4')                                                                                                                        
#    ax1.grid()                                                                                                                                    
#    ax2.grid()                                                                                                                                    
#    ax3.grid()                                                                                                                                    
#    ax4.grid()                                                                                                                                    
#    f.show()
    
    xx = np.arange(0,1180,20)                                                                                                                           
    yy = np.arange(0,1280,20)                                                                                                                     
    xxx,yyy = np.meshgrid(xx,yy) 
    
#    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#    scale = 0.2
#    qv2 = ax1.quiver(yyy, xxx,  cy1(xxx,yyy)-yyy, cx1(xxx,yyy)-xxx,scale_units='xy', angles='xy', scale=scale)
#    ax2.quiver(yyy, xxx, cy2(xxx,yyy)-yyy, cx2(xxx,yyy)-xxx,scale_units='xy', angles='xy', scale=scale)
#    ax3.quiver(yyy, xxx,  cy3(xxx,yyy)-yyy, cx3(xxx,yyy)-xxx,scale_units='xy', angles='xy', scale=scale)
#    ax4.quiver(yyy, xxx,  cy4(xxx,yyy)-yyy, cx4(xxx,yyy)-xxx,scale_units='xy', angles='xy', scale=scale)
#    ax1.quiverkey(qv2, 1.1,-0.2, 3, 'scale: 3 pixels', coordinates='axes', color='r')
#    ax1.quiverkey(qv2, 1.1,-1.5, 10, 'scale: 10 pixels', coordinates='axes', color='r')
#    ax1.plot(x,y,'*', c='b')
#    ax2.plot(x,y,'*', c='b')
#    ax3.plot(x,y,'*', c='b')
#    ax4.plot(x,y,'*', c='b')
#    ax1.set_title('Deg 1')
#    ax2.set_title('Deg 2')
#    ax3.set_title('Deg 3')
#    ax4.set_title('Deg 4')
#    f.suptitle(name)
#    f.show()
    
    return [cx1,cy1]

m159 = fits.open('502352_matchV.fits')[1]
m121 = fits.open('475788_matchV.fits')[1]
m161 = fits.open('486607_matchV.fits')[1]
m119 = fits.open('500513_matchV.fits')[1]
m119b = fits.open('500633_matchV.fits')[1]
m119c = fits.open('500731_matchV.fits')[1]
m159b = fits.open('502179_matchV.fits')[1]
m159c = fits.open('502282_matchV.fits')[1]


#AnalyseDistortion(m119,scale = 0.2,name='mask 1')
#AnalyseDistortion(m119b,scale = 0.2,name='mask 1b')
#AnalyseDistortion(m119c,scale = 0.2,name='mask 1c')
AnalyseDistortion(m159,scale = 0.2,name='mask 4')
#AnalyseDistortion(m159b,scale = 0.2,name='mask 4b')
#AnalyseDistortion(m159c,scale = 0.2,name='mask4c')
#AnalyseDistortion(m121,scale = 0.2,name='mask 3')
#AnalyseDistortion(m161,scale = 0.2,name='mask 2')
#




#
#
#plt.figure()
#plt.plot(a['X'],a['Y'],'*', label = 'Theoric position')
#plt.plot(a['xcentroid'],a['ycentroid'],'.',label='Star in guider camera')
##plt.plot(poly.polyval2d(x,y,cy1),poly.polyval2d(x,y,cx1),'.',c='black',label = 'polyfit')
#plt.plot(poly.polyval2d(x,y,cx1),poly.polyval2d(x,y,cy1),'.',c='black',label = 'polyfit')
#plt.legend(loc='lower left', borderaxespad=0.)
#plt.xlabel('Pixel Y on GOBC')
#plt.ylabel('Pixel Xon GOBC')
#plt.show()



############################                      
#cx1 = PolyFit2d(np.hstack((x1,x2)),np.hstack((y1,y2)),np.hstack((valx1,valx2)),[1,1])
#cy1 = PolyFit2d(np.hstack((x1,x2)),np.hstack((y1,y2)),np.hstack((valy1,valy2)),[1,1])
#cx2 = PolyFit2d(np.hstack((x1,x2)),np.hstack((y1,y2)),np.hstack((valx1,valx2)),[2,2])
#cy2 = PolyFit2d(np.hstack((x1,x2)),np.hstack((y1,y2)),np.hstack((valy1,valy2)),[2,2])
#cx3 = PolyFit2d(np.hstack((x1,x2)),np.hstack((y1,y2)),np.hstack((valx1,valx2)),[3,3])
#cy3 = PolyFit2d(np.hstack((x1,x2)),np.hstack((y1,y2)),np.hstack((valy1,valy2)),[3,3])
#cx4 = PolyFit2d(np.hstack((x1,x2)),np.hstack((y1,y2)),np.hstack((valx1,valx2)),[4,4])
#cy4 = PolyFit2d(np.hstack((x1,x2)),np.hstack((y1,y2)),np.hstack((valy1,valy2)),[4,4])
#############################




#plt.figure()
#plt.plot(x1,y1,'x',c='r',label = 'test 1', color='r')
#plt.plot(x2,y2,'+',c='r',label = 'test 2',color='b')
#plt.show()
#
#plt.figure()
#plt.plot(valx1,valy1,'x',c='r',label = 'test 1',color='r')
#plt.plot(valx2,valy2,'+',c='r',label = 'Ftest 2',color='b')
#plt.show()
#
#
#
#cx1 = PolyFit2d(valx2,valy2,x2,[1,1])
#cy1 = PolyFit2d(valx2,valy2,y2,[1,1])
#
#cx2 = PolyFit2d(valx2,valy2,x2,[2,2])
#cy2 = PolyFit2d(valx2,valy2,y2,[2,2])
#
#cx3 = PolyFit2d(valx2,valy2,x2,[3,3])
#cy3 = PolyFit2d(valx2,valy2,y2,[3,3])
#
#plt.figure()
#plt.plot(x2,y2, '*',label='coord arcsec')
#plt.plot(valx2,valy2, '*',label='pixels')
#plt.plot(-723.5,-19.7,'P',label='pixel etoile didier')
#plt.plot(cx2(valx2,valy2),cy2(valx2,valy2),'.',label='pixel etoile didier')
##plt.plot(cx1(-723.5,-19.7),cy1(-723.5,-19.7),'X',label='coord fit1 etoile')
##plt.plot(cx2(-723.5,-19.7),cy2(-723.5,-19.7),'X',label='coord fit2 etoile')
##plt.plot(cx3(-723.5,-19.7),cy3(-723.5,-19.7),'X',label='coord fit3 etoile')
#plt.axis('equal')
#plt.legend(loc='lower right', borderaxespad=0.)
#plt.xlabel('Azimuth')
#plt.ylabel('Altitude')
#plt.title('Position of the stars in the sky during the test')
#plt.legend(bbox_to_anchor=(0.7, 0.2), loc=2, borderaxespad=0.)
#plt.show() 
#




##########Mesures d eprecision on prend les polynomes de degre 2

#
#
#
#from mpl_toolkits.mplot3d import axes3d
#import matplotlib.pyplot as plt
#from matplotlib import cm
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#xx = np.arange(-0.3*3600,0*3600,20)
#yy = np.arange(-0.15*3600,0.15*3600,20)
#xxx,yyy = np.meshgrid(xx,yy)
#ax.plot_surface(xxx, yyy, xxx*cx1.jacobian(0,0)[0] - cx2(xxx,yyy), rstride=8, cstride=8, alpha=0.3, color='r')
##ax.plot_surface(xxx, yyy, xxx*cx1.jacobian(0,0)[0] - cx3(xxx,yyy), rstride=8, cstride=8, alpha=0.3)
##ax.plot_surface(xxx, yyy, xxx*cx1.jacobian(0,0)[0] - cx1(xxx,yyy), rstride=8, cstride=8, alpha=0.3)
#ax.plot_surface(xxx, yyy, yyy*cy1.jacobian(0,0)[1]-cy2(xxx,yyy), rstride=8, cstride=8, alpha=0.3,color='b')
##cset = ax.contour(xxx, yyy, xxx*cx1.jacobian(0,0)[0] - cx2(xxx,yyy), zdir='z', offset=-100, cmap=cm.coolwarm)
##cset = ax.contour(xxx, yyy, xxx*cx1.jacobian(0,0)[0] - cx2(xxx,yyy), zdir='x', offset=-40, cmap=cm.coolwarm)
##cset = ax.contour(xxx, yyy, xxx*cx1.jacobian(0,0)[0] - cx2(xxx,yyy), zdir='y', offset=40, cmap=cm.coolwarm)
#ax.scatter(x1,y1,x1*cx1.jacobian(0,0)[0] - cx2(x1,y1),'*', c='b')
#ax.scatter(x1,y1,y1*cy1.jacobian(0,0)[1] - cy2(x1,y1),'*', c='b')
#ax.set_xlabel('X - Guider Camera')
#ax.set_ylabel('Y - Guider Camera')
#ax.set_zlabel('Distortion difference (pix)')
#plt.title('Distortion difference on X (red) and Y (blue)  axis to center \n magnigication  with a deg 2 polynomial fit in (x,y)')
##ax.set_zlim(-100, 100)
#plt.show()
#