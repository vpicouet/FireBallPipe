#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 19:29:40 2018

@author: Vincent
"""

from DS9Utils import *
def parse_data(data, nan=np.nan, map=map, float=float):
    vals = []
    xy = []
    for s in data.split('\n'):
        coord,val = s.split('=')
        val = val.strip() or nan
        xy.append(list(map(float, coord.split(','))))
        vals.append(float(val))
    vals = np.array(vals)
    xy = np.floor(np.array(xy)).astype(int)
    x,y = xy[:,0], xy[:,1]
    w = x.ptp() + 1
    h = y.ptp() + 1
    arr = np.empty((w, h))
    X = np.empty((w, h))
    Y = np.empty((w, h))
    indices = x-x.min(), y-y.min()
    arr[indices] = vals
    X[indices] = x
    Y[indices] = y
    return X.T, Y.T, arr.T

def process_region(region, win):
    name, info = region.split('(')
    coords = [float(c) for c in info.split(')')[0].split(',')]
    if name == 'box':
        xc,yc,w,h,angle = coords
        dat = win.get("data physical %s %s %s %s no" % (xc-old_div(w,2.),yc-old_div(h,2.),w,h))
        X,Y,arr = parse_data(dat)
        box = namedtuple('Box', 'data xc yc w h angle')
        return box(arr, xc, yc, w, h, angle)
    elif name == 'circle':
        xc,yc,r = coords
        dat = win.get("data physical %s %s %s %s no" % (xc-r,yc-r,2*r,2*r))
        X,Y,arr = parse_data(dat)
        Xc,Yc = np.floor(xc), np.floor(yc)
        inside = (X - Xc)**2 + (Y - Yc)**2 <= r**2
        circle = namedtuple('Circle', 'data databox inside xc yc r')
        return circle(arr[inside], arr, inside, xc, yc, r)
#    elif name == 'annulus':
#        xc,yc,r1,r2 = coords
#        w = 2*r2
#        h = 2*r2
#        dat = win.get("data physical %s %s %s %s no" % (xc-r2,yc-r2,w,h))
#        X,Y,arr = parse_data(dat)
#        Xc,Yc = np.floor(xc), np.floor(yc)
#        inside = between((X - Xc)**2 + (Y - Yc)**2, r1**2, r2**2)
#        annulus = namedtuple('Annulus', 'data databox inside xc yc r1 r2')
#        return annulus(arr[inside], arr, inside, xc, yc, r1, r2)
    elif name == 'ellipse':
        if len(coords) == 5:
            xc,yc,a2,b2,angle = coords
        else:
            xc,yc,a1,b1,a2,b2,angle = coords
        w = 2*a2
        h = 2*b2
        dat = win.get("data physical %s %s %s %s no" % (xc-a2,yc-b2,w,h))
        X,Y,arr = parse_data(dat)
        Xc,Yc = np.floor(xc), np.floor(yc)
        inside = (old_div((X - Xc),a2))**2 + (old_div((Y - Yc),b2))**2 <= 1
        if len(coords) == 5:
            ellipse = namedtuple('Ellipse',
                                 'data databox inside xc yc a b angle')
            return ellipse(arr[inside], arr, inside, xc, yc, a2, b2, angle)

        inside &= (old_div((X - Xc),a1))**2 + (old_div((Y - Yc),b1))**2 >= 1
        annulus = namedtuple('EllipticalAnnulus',
                             'data databox inside xc yc a1 b1 a2 b2 angle')
        return annulus(arr[inside], arr, inside, xc, yc, a1, b1, a2, b2, angle)
    else:
        raise ValueError("Can't process region %s" % name)

def getregion(win, debug=False):
    """ Read a region from a ds9 instance.

    Returns a tuple with the data in the region.
    """
    rows = win.get("regions")
    rows = [row for row in rows.split('\n') if row]
    if len(rows) < 3:
        print( "No regions found")
        return
    units = rows[2]
    assert units == 'physical'
    if debug:
        print (rows[4])
        if rows[5:]:
            print('discarding %i regions' % len(rows[5:]) )
    return process_region(rows[-1], win)




def throughfocus(center, files,x = np.linspace(11.95,14.45,11)[::-1][3:8],fibersize=100,center_type=None, SigmaMax = 4):
    fwhm = []
    EE50 = []
    EE80 = []    
    for file in files:
        image = fits.open(file)[0].data
        d = plot_rp2_convolved_wo_latex(image,center=center,fibersize=fibersize,center_type=center_type, SigmaMax = SigmaMax)
        fwhm.append(d['FWHM'])
        EE50.append(d['EE50'])
        EE80.append(d['EE80'])
    f = lambda x,a,b,c: a * (x-b)**2 + c#a * np.square(x) + b * x + c
    opt1,cov1 = curve_fit(f,x,fwhm)
    opt2,cov2 = curve_fit(f,x,EE50)
    opt3,cov3 = curve_fit(f,x,EE80)
    xtot = np.linspace(x.min(),x.max(),200)
    fig, axes = plt.subplots(3, 1, figsize=(8,8))
    axes[0].plot(x,fwhm,'-o')
    axes[0].plot(xtot,f(xtot,*opt1),linestyle='dotted')
    axes[0].plot(np.ones(2)*xtot[np.argmin(f(xtot,*opt1))],[min(fwhm),max(fwhm)])
    axes[0].set_ylabel('FWHM')
    axes[0].set_xlabel('Best actuator = %0.3f'%(xtot[np.argmin(f(xtot,*opt1))]))
    axes[1].plot(x,EE50,'-o')
    axes[1].plot(xtot,f(xtot,*opt2),linestyle='dotted')
    axes[1].set_ylabel('EE50')
    axes[1].set_xlabel('Best actuator = %0.3f'%(xtot[np.argmin(f(xtot,*opt2))]))
    axes[1].plot(np.ones(2)*xtot[np.argmin(f(xtot,*opt2))],[min(EE50),max(EE50)])
    axes[2].plot(x,EE80,'-o')
    axes[2].plot(xtot,f(xtot,*opt3),linestyle='dotted')
    axes[2].set_ylabel('EE80')
    axes[2].set_xlabel('Best actuator = %0.3f'%(xtot[np.argmin(f(xtot,*opt3))]))
    axes[2].plot(np.ones(2)*xtot[np.argmin(f(xtot,*opt3))],[min(EE80),max(EE80)])
    fig.tight_layout()    
    plt.show()
    return fwhm, EE50, EE80


if __name__ == '__main__':

    print('''\n\n\n\n      START IMAGE ANALYSIS \n\n\n\n''')
    
    #if sys.version_info.major == 3:
    n1 = input("Number of first image[ex:10 or press ENTER take all images from folder]: ")
    n2 = input("Number of first image[ex:21 or press ENTER take all images from folder]: ")
    #if sys.version_info.major == 2:
    #    n1 = input("Number of first image[ex:10 or press ENTER take all images from folder]: ")
    #    n2 = input("Number of first image[ex:21 or press ENTER take all images from folder]: ")
    
    d = DS9()
    filename = d.get("file")
    files = []
    path = []
    if (type(n1)==int) & (type(n2)==int):
        for number in np.arange(n1,n2+1):
            #path = os.path.dirname(filename) + '/image%06d.fits'%(number)
            path.append(os.path.dirname(filename) + '/image%06d.fits'%(number))
            x = np.arange(n1,n2+1)
    else:
        path = glob.glob(os.path.dirname(filename) + '/*.fits')
        x = np.arange(len(path))
    #    with fits.open(path) as f:
    #        files.append(f[0].data)        
    #    os.path.dirname(filename)
    path = np.sort(path)
    print(path)
    
    
    a = getregion(d)
    rp = plot_rp2_convolved_wo_latex(fits.open(filename)[0].data,center = [np.int(a.xc),np.int(a.yc)],fibersize=100)
    print('''\n\n\n\n     Centring on barycentre of the DS9 image (need to be close to best focus) : %0.1f, %0.1f --> %0.1f, %0.1f \n\n\n\n'''%(a.xc,a.yc,rp['Center'][0],rp['Center'][1]))
    throughfocus(center = rp['Center'], files=path,x = x,fibersize=100,center_type=None,SigmaMax=6)