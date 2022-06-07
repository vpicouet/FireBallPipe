#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:59:42 2018

@author: Vincent
"""
import os, sys, glob
from astropy.io import fits
import time   
import numpy as np
import timeit
#Give in entry fist image number and last
#Will stack every image in this range 1-2, 1-2-3... while they are descending, if they are not here we wait
#each time a stack is done, we add it's right and left part and save it with stack number-number name
#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/ReductionPipeline/TestStacking/BighImages/'

def CropImage(path,names='image??????.fits'):
    for file in glob.glob(path+names):
        print(file)
        fitsImage = fits.open(file)[0]
        fitsImage.data = fitsImage.data[:,1053:2133]
        if 'NAXIS3' in fitsImage.header:
            fitsImage.header.remove('NAXIS3') 
        fitsImage.writeto(os.path.dirname(os.path.dirname(path)) + '/Croped/' + os.path.basename(file)[:-5] + '.crop.thresh.fits',overwrite=True)
        return 
#CropImage(path=path,names='image??????.fits')

def ReturnFiles(paht, n0, n):
    files2stack = []
    for file in glob.glob(path + names):
        if (int(os.path.basename(file)[5:11])>=n0) & (int(os.path.basename(file)[5:11])<=n):
            files2stack.append(file)
    #10print(files2stack)
    return files2stack, n0, n


def StackFlight(files2stack, n0, n):
    start = timeit.default_timer()
    Image0 = fits.open(files2stack[0])[0]
    lx,ly = Image0.data.shape
    stack = np.zeros((lx,ly,len(files2stack)))
    for i,file in enumerate(files2stack):
        with fits.open(file) as f:
            stack[:,:,i] = f[0].data
    stack = np.mean(stack,axis=2)
    stack = AddParts2Image(stack)
    Image0.data = stack
    name = os.path.dirname(files2stack[0])+'/stack/' + 'image%06d-%06d.fits'%(n0,n)
    Image0.writeto(name)
    stop = timeit.default_timer()
    print("""
    ********************************************************************
                        Stack saved: {}\n duration= {}s      
    ********************************************************************
    """.format(name, stop - start)) 
    return stack
        
def AddParts2Image(array,size=[1053,3216-2133]):
    lx,ly = array.shape
    left = np.ones((lx,size[0]))*array.mean()
    right = np.ones((lx,size[1]))*array.mean()
    BigStack = np.hstack((left,array,right))
    return BigStack
    
if __name__ == '__main__':
    print(__file__)
    entry = sys.argv[1]#'/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/ReductionPipeline/TestStacking/Croped/image000000.crop.thresh.fits'#sys.argv[1]
    path = os.path.dirname(entry) + '/'
    name = os.path.basename(entry)
    extension = name[11:]
    
    n1 = raw_input("Stacking: First Image : ")
    n2 = raw_input("Stacking: Last Image : ")
    if n2=='*':
        n2=1000
    n1,n2 = int(n1), int(n2)
    
    names = 'image??????' + extension
    n = n1+1
    while n<n2:
        print(n,n2)
        if os.path.isfile(path+'stack/' + 'image%06d-%06d.fits'%(n1,n)):
            print('Stack already existing, stacking next one')
            n += 1
        else:
            print(n,n2)
            for number in np.arange(n,n2+1):
                #print('blable')
                if os.path.isfile(path + 'image%06d'%(number)+extension):
                    print('Stack not existing, applying stacking...')
                    files2stack, n1 ,number  = ReturnFiles(path, n1, number)
                    stackedImage = StackFlight(files2stack, n1, number)
                    print('Stack created.')
                    n = number + 1
                    break
                else:
                    print('File not dowloaded yet, waiting for data...') 
            time.sleep(1)
    print('All stack created.')
