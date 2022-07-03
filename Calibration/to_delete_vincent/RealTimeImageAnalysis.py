#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 14:53:08 2018

@author: Vincent
"""

#path = '/Users/Vincent/Nextcloud/Work/FTS2018_FLIGHT/test/InstrumentCentering_180819/MaksTF/'
#for passt in glob.glob(path+'*'):
#    for repertory in glob.glob(passt+'/*'):
#       # for repertory in glob.glob(repertories+'*'):
#        files = glob.glob(repertory + '/*.fits')
#        print(files[5])
#        
#        F = Focus(filename = files[5], HumanSupervision=False, source='Zn', shape='gaussian', windowing=False, peak_threshold=50,plot=False)
#        for i in range(len(F.table)):
#           # i=0
#            path = glob.glob(os.path.dirname(files[5])+'/*fits')
#            x = np.arange(len(path))
#            throughfocus(center = [F.table['X_IMAGE'][i],F.table['Y_IMAGE'][i]], 
#                         files=path,x = x,fibersize=0,
#                         center_type=None,SigmaMax=6,Plot=True)
#csvfiles = []
#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180821/thrufocus/'
#for passt in glob.glob(path+'*'):
#    for repertory in glob.glob(passt+'/*'):
#        name = repertory + '/Throughfocus.csv'
#        csvfiles.append(name)
#t = Table.read(csvfiles[0])
#t.remove_rows([0, 1, 2])
#for file in csvfiles:
#    t0 = Table.read(files)
#    t = vstack((t,t0))
#t.writeto(path + 'TotalThroughfocus.csv')
#    

path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180823/xycalib/Catalogs'
def mergeCat(path):
    cats = glob.glob(path + '/Analysis_*.csv')
    t = Table.read(cats[0])
    for file in cats[1:]:
        print(file)
        t0 = Table.read(file)
        t = vstack((t,t0))
    t.write(path + '/AnalysisTotal.csv',overwrite=True)
    return
#mergeCat(path)


def create_PA(A=15.45,B=13.75,C=14.95,pas=0.15,nombre=11):
    a = np.linspace(A-int(nombre/2)*pas, A+int(nombre/2)*pas, nombre)
    b = np.linspace(B-int(nombre/2)*pas, B+int(nombre/2)*pas, nombre)
    c = np.linspace(C-int(nombre/2)*pas, C+int(nombre/2)*pas, nombre)
    print(ENCa)
    print(ENCb)
    print(ENCc)
    return a[::-1],b[::-1],c[::-1]
#    ENCa, ENCb, ENCc = create_PA()

def ENCA(x):
    a = (ENCa[10]-ENCa[0])/(10) * x + ENCa[0]
    b = (ENCb[10]-ENCb[0])/(10) * x + ENCb[0]
    c = (ENCc[10]-ENCc[0])/(10) * x + ENCc[0]
    return a,b,c

path1 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180825/xycalib/image000297.fits'
path2 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180825/xycalib/image000380.fits'

def CompareImages(path1,path2, mask='F3'):
    if path1[-3:]  == 'csv':
        from astropy.table import Table
        F1table = Table.read(path1)
        F2table = Table.read(path2)
    else:        
        from focustest import Focus
        F1 = Focus(filename = path1,  quick=False, windowing=True, mask=mask, HumanSupervision=False, source='Zn', threshold = [7], fwhm = [9,12.5],shape='slits', peak_threshold=50,plot=False)
        F2 = Focus(filename = path2,  quick=False, windowing=True, mask=mask, HumanSupervision=False, source='Zn', threshold = [7], fwhm = [9,12.5],shape='slits', peak_threshold=50,plot=False)
        F1table, F2table = F1.table, F2.table
    if len(F1table) < len(F2table):
        F1table, F2table = F2table, F1table
    F1n, F2n, = F1table.copy(), F2table.copy()
    F1n.remove_rows(np.arange(len(F1n)))
    F2n.remove_rows(np.arange(len(F2n)))
    for slit in F1table[F1table['id_slit']>0]:
        if len( F2table[(F2table['id_slit']==slit['id_slit']) & (F2table['wavelength']==slit['wavelength'])])==1:
            print('match')
            F1n.add_row(slit)
            F2n.add_row(F2table[(F2table['id_slit']==slit['id_slit']) & (F2table['wavelength']==slit['wavelength'])][0])
    
    x, y = 'X_IMAGE', 'Y_IMAGE'
    #x, y = 'xcentroid', 'ycentroid'
    Deltax = F1n[x] - F2n[x]
    Deltay = F1n[y] - F2n[y]
    id = (abs(Deltax)<np.percentile(abs(Deltax),97)) & (abs(Deltay)<np.percentile(abs(Deltay),97))
    Deltax, Deltay = Deltax[id], Deltay[id]
    delta = np.array([Deltax,Deltay]).T    
    position = np.array([F2n[x][id],F2n[y][id]]).T 
    theta_rad, deltax, deltay = ComputeSmallRotationOffset(delta,position)

    
    plt.figure()
    Q = plt.quiver(position[:,0], position[:,1], delta[:,0],  delta[:,1])
    qk = plt.quiverkey(Q, 0.5, 0.2, 10, '10 pixels', color='r', labelpos='E',coordinates='figure')
    plt.show()
    return theta_rad, deltax, deltay

path1 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180825/xycalib/image000296_table.csv'
path2 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180825/xycalib/image000381_table.csv'
path2 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/image-000025-000034-Zinc-with_dark-161-stack_table.csv'
CompareImages(path1,path2,mask='F2')

path1 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180825/xycalib/image000297_table.csv'
path2 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180825/xycalib/image000380_table.csv'
path2 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018/AIT-Optical-FTS-201805/180612/image-000075-000084-Zinc-with_dark-121-stack_table.csv'
CompareImages(path1,path2,mask='F3')


path1 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180826/xycalib/image000000_table.csv'
path2 = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/AIT-Optical-FTS-2018-Flight/XYCalibration/DiffuseMasksIllumination180823/StackedImage_24-43-NoDark_table.csv'
CompareImages(path1,path2,mask='F1')
   
def ComputeSmallRotationOffset(delta,position):
    data = np.concatenate((delta[:,0], delta[:,1]))
    row_x = np.hstack((-position[:,[1]], np.ones((len(position),1)), np.zeros((len(position),1)) )) # xn - x =  x dgamma - y theta + dx
    row_y = np.hstack((position[:,[0]], np.zeros((len(position),1)), np.ones((len(position),1)) )) # yn - y =  y dgamma + x theta      + dy
    mat = np.vstack((row_x, row_y))
    matinv =  np.linalg.pinv(mat)
    sol = matinv.dot(data)
    theta_rad = sol[0]
    deltax = sol[1]
    deltay = sol[2]
    theta = theta_rad*180/np.pi*60 #arcmin
    print("theta: {} arcmin\ndx: {} pixel\ndy: {} pixel".format(theta, deltax, deltay))
    covar = matinv.dot(matinv.T)
    # accuracy, assuming 1 arcsec measurement error
    print("variances: {}\n".format(np.sqrt(np.diag(covar))/3600*[180/np.pi*60, 3600, 3600])) #
    #residual
    data_new = mat.dot(sol)
    #print("residuals in arcsec:", (data_new - data).reshape((len(position),2))*3600)
    return theta_rad, deltax, deltay


from astropy.table import Table
from DS9Utils import *
from focustest import Focus
#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180821/thrufocus/'
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/'
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180831/thrufocus'

for path in glob.glob(path+'*')[2:]:
    for repertory in glob.glob(path+'/*/'):
        #files = glob.glob(repertory + '/*.fits')
        #print(files[5])
        #F = Focus(filename = files[5], HumanSupervision=False, source='Zn', shape='gaussian', windowing=False, peak_threshold=50,plot=False)
        files = glob.glob(repertory + '/*table.csv')
        try:
            table = Table.read(files[0])
        except IndexError:
            files = glob.glob(repertory + '/*.fits')
            F = Focus(filename = files[8], HumanSupervision=False, source='Zn', shape='gaussian', windowing=False, peak_threshold=50,plot=False)
        try:
            files = glob.glob(repertory + '/*table.csv')
            table = Table.read(files[0])
    
            for i in range(len(table)):
               # i=0
                path = glob.glob(repertory+'/*fits')
                #x = np.arange(len(path))
                throughfocus(center = [table['X_IMAGE'][i],table['Y_IMAGE'][i]], 
                             files=path,fibersize=0,
                             center_type=None,SigmaMax=6,Plot=True, ENCa_center=16.45, pas=0.15)
        except IndexError:
            pass



def mergeCatThroughfocus(path,stage):
    from astropy.table import Table, vstack
    csvfiles = []
    for files in glob.glob(path+'*'):
        if stage == 2:
            for files in glob.glob(files+'/*'):
                name = files + '/Throughfocus.csv'
                csvfiles.append(name)
        if stage == 1:
            name = files + '/Throughfocus.csv'
            csvfiles.append(name)            
    t = Table.read(csvfiles[0])

    t.remove_rows(np.arange(len(t)))
    for file in csvfiles:
        t0 = Table.read(file)
        t = vstack((t,t0))
    t.write(path + 'TotalThroughfocus.csv')
    return t
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180831/thrufocus/'
cat = mergeCatThroughfocus(path,stage=1)


path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180831/FOCUS180831/F1/'
cat = mergeCatThroughfocus(path,stage=1)
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180831/FOCUS180831/F2/'
cat = mergeCatThroughfocus(path,stage=1)
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180831/FOCUS180831/F3/'
cat = mergeCatThroughfocus(path,stage=1)
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180831/FOCUS180831/F4/'
cat = mergeCatThroughfocus(path,stage=1)


path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180823/xycalib/Catalogs'
def mergeCat(path):
    """
    Veirfier que ya pas de doublons
    """
    cats = glob.glob(path + '/Analysis_*.csv')
    t = Table.read(cats[0])
    for file in cats[1:]:
        print(file)
        t0 = Table.read(file)
        t = vstack((t,t0))
    t.write(path + '/AnalysisTotal.csv',overwrite=True)
    return
mergeCat(path)


#
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/snape/180821/thrufocus/'
path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/data/guider/180821/MaksTF/'
def delete(extension = '/*.png', path=path):
    for passt in glob.glob(path+'*'):
        for repertory in glob.glob(passt+'/*'):
           # for repertory in glob.glob(repertories+'*'):
            files = glob.glob(repertory + extension)
            print(files)
            #files = glob.glob(repertory + '/Throughfocus*.csv')
            for file in files:
                print(file)
                print ('Removed')
                os.remove(file)
#                answer = raw_input('Delete these files?[y/n]')
#                if answer == '':
#                    os.remove(file)
#                    print ('Removed')
    return
                
delete(extension = '/*.png', path=path)               

