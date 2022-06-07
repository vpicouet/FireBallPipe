#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:59:42 2018

@author: Vincent


This stacking function aims at stacking images coming from 1Mbit trasnmitter for FB-2 flight
Give the first and last image you want to stack ('*'=1000 for inf) even if this image is not taken yet.
The code deals with:
    - Delays if the image has not been taken yet
    - If the telemetry could not downlinked an image because it's always the last one which is downloaded
    - The stacked saved are resized to nominal size to make easier the real time DS9 analysis as only 1/3 of the images are downlinked 
Just give the path of an image as argumebt of the function and you will be asked to enter firt and last number to stack
"""

import os, glob
from astropy.io import fits
import time   
import numpy as np
import timeit
from shutil import copyfile

def CropImage(path,names='image??????.fits'):
    for file in glob.glob(path+names):
        print(file)
        fitsImage = fits.open(file)[0]
        fitsImage.data = fitsImage.data[:,1053:2133]
        if 'NAXIS3' in fitsImage.header:
            fitsImage.header.remove('NAXIS3') 
        fitsImage.writeto(os.path.dirname(os.path.dirname(path)) + '/Croped/' + os.path.basename(file)[:-5] + '.crop.thresh.fits',overwrite=True)
        #CropImage(path=path,names='image??????.fits')
        return 

def ReturnFiles(paht, n0, n):
    files2stack = []
    for file in glob.glob(path + names):
        if (int(os.path.basename(file)[5:11])>=n0) & (int(os.path.basename(file)[5:11])<=n):
            files2stack.append(file)
    #10print(files2stack)
    return files2stack, n0, n


def StackFlight(files2stack, n0, n):
    """First version, take all the images to stack them.
    The function takes several seconds for tens of images which is too long... -> v1
    """
    start = timeit.default_timer()
    Image0 = fits.open(files2stack[0])[0]
    lx,ly = Image0.data.shape
    stack = np.zeros((lx,ly,len(files2stack)))
    for i,file in enumerate(files2stack):
        with fits.open(file) as f:
            stack[:,:,i] = f[0].data
    stack = np.sum(stack,axis=2)
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
  
def StackFlight_v1(stack1, stack2, n):
    """New version to not stack all the images. Just take the last image stacked
    and add the last image downloaded
    """
    start = timeit.default_timer()
    Image0 = fits.open(stack1)[0]
    lx,ly = Image0.data.shape
    try:
        stack = Image0.data + fits.open(stack2)[0].data
        stack = AddParts2Image(stack)
    except ValueError:
        stack = Image0.data + AddParts2Image(fits.open(stack2)[0].data)
    Image0.data = stack
    name = stack1[:-11] +'%06d.fits'%(n)
    Image0.writeto(name)
    stop = timeit.default_timer()
    print("""
    ********************************************************************
                        Stack saved: {}\n duration= {}s      
    ********************************************************************
    """.format(name, stop - start)) 
    return    
      
def AddParts2Image(array,size=[1053,3216-2133]):
    """Add the right and left parts of the images which have been discarded because of the 
    1Mbit download which takes too long (longer than the 30s exposure time)
    """
    lx,ly = array.shape
    left = np.zeros((lx,size[0]))#+array.mean()
    right = np.zeros((lx,size[1]))#+array.mean()
    BigStack = np.hstack((left,array,right))
    return BigStack
    
if __name__ == '__main__':
    print(__file__)
    entry = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/ReductionPipeline/TestStacking/Croped/image000000.crop.thresh.fits'#sys.argv[1]
    path = os.path.dirname(entry) + '/'
    name = os.path.basename(entry)
    extension = name[11:]
    
    n1 = 100#raw_input("Stacking: First Image : ")
    n2 = 120#raw_input("Stacking: Last Image : ")
    if n2 == '*':
        n2 = 1000
    n1,n2 = int(n1), int(n2)
    names = 'image??????' + extension
    try:
        os.makedirs(path+'stack')
    except FileExistsError:
        pass
    copyfile(path + 'image%06d' % (n1) + extension, path+'stack/' + 'image%06d-%06d.fits'%(n1,n1))


    n = n1 + 1
    while n < n2:
        print(n,n2)
        if os.path.isfile(path + 'stack/' + 'image%06d-%06d.fits' % (n1, n)):
            print('Stack already existing, stacking next one')
            n += 1
            
            
        else:
            print(n,n2)
            for number in np.arange(n,n2+1):
                lastfile = path + 'image%06d' % (number) + extension
                
                
                if os.path.isfile(lastfile):
                    print('Stack not existing, Searching for last stack...')
                    stacked_files = glob.glob(path + 'stack/'  + 'image%06d-??????.fits'%(n1))
                    numbers = [int(os.path.basename(file)[-11:-5]) for file in stacked_files]
                    print('last_skacked_file = ', np.argmax(numbers))
                    if os.path.isfile(path + 'stack/' + 'image%06d-%06d.fits' % (n1, number)) :
                        n = number + 1
                    else:
                        last_skacked_file = stacked_files[np.argmax(numbers)]
                        stackedImage = StackFlight_v1(stack1 = last_skacked_file, stack2=lastfile, n=number)
                        print('Stack created.')
                        n = number + 1
                    break
                
                
                else:
                    print('File not dowloaded yet, waiting for data...') 
            time.sleep(1)
    print('All stack created.')
