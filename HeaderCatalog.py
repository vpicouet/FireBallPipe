#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 22:18:58 2018

@author: Vincent
"""
from __future__ import print_function
import os, sys
import datetime
import time
import numpy as np
from astropy.table import Table
import subprocess

def HeaderCatalog(path):
    with open(path, 'r') as myfile:
        header=myfile.read().replace('\n', 'newline')
    a = header.split('newline')
    fields = []
    for line in a[:91]:
        line = line.split('=')
        line = line[0].split(' ')
        fields.append(line[0])
    print (fields)
    
    values = np.zeros(((len(a)+1)/len(fields),len(fields)))
    values_flat = np.zeros(len(a)-1)
    for i,line in enumerate(a[:-1]):
        line = line.replace('/ ','=').split('=')
        #print(line[1])
        try:
            values_flat[i] = float(line[1]) 
        except ValueError:
            #hour = line[1].replace("'",'T').split('T')[2]
            #values_flat[i] = datetime.timedelta(hours=sec.tm_hour,minutes=sec.tm_min,seconds=sec.tm_sec).total_seconds()
            DT = line[1].split("'")[1]
            sec = datetime.datetime.strptime(DT, "%Y-%m-%dT%H:%M:%S")
            time.mktime(sec.timetuple())
            values_flat[i] = time.mktime(sec.timetuple())
    values = values_flat.reshape(values.shape)
    t = Table(values, names=fields)
    t.write(os.path.dirname(path) + '/HeaderCatalog.csv',overwrite=True) 
    print('Table saved :', os.path.dirname(path) + '/HeaderCatalog.csv')
    return t    

#path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/SkyTest110918/Guider/stack22364375.fits'
if __name__ == '__main__':
    
    path = sys.argv[1]
#    path = '/Users/Vincent/Nextcloud/FIREBALL/TestsFTS2018-Flight/E2E-AIT-Flight/SkyTest110918/Guider/stack22364375.fits'

    function="""files=$(l\\s %s/stack*fits);for file in $files;do fold $file | head -105 | grep '\|IMGTAG\|EXPOSURE\|FRAMESTA\|FRAMEEND\|IMGCOUNT\|ROTENC\|LINAENC\|LINBENC\|LINCENC\|RA_DTU\|DEC_DTU\|ROLL_DTU\|AZ\|EL\|MROT\|CCDTEMP\|CAMTEMP\|PWSTEMP\|TRTD\|PRESSURE\|VALID0\|USE0\|TX0\|TY0\|CX0\|CY0\|FLUX0\|SIGMAX0\|SIGMAY0\|VALID1\|USE1\|TY1\|CX1\|CY1\|FLUX1\|SIGMAX1\|SIGMAY1\|VALID2\|USE2\|TX2\|TY2\|CX2\|CY2\|FLUX2\|SIGMAX2\|SIGMAY2\|VALID3\|USE3\|TX3\|TY3\|CX3\|CY3\|FLUX3\|SIGMAX3\|SIGMAY3\|VALID4\|USE4\|TX4\|TY4\|CX4\|CY4\|FLUX4\|SIGMAX4\|SIGMAY4\|VALID5\|USE5\|TX5\|TY5\|CX5\|CY5\|FLUX5\|SIGMAX5\|SIGMAY5\|VALID6\|USE6\|TX6\|TY6\|CX6\|CY6\|FLUX6\|SIGMAX6\|SIGMAY6\|VALID7\|USE7\|TX7\|TY7\|CX7\|CY7\|FLUX7\|SIGMAX7\|SIGMAY7\|DATE';done >> %s/fitsheader.txt;"""%(os.path.dirname(path),os.path.dirname(path))
    #function=r"""files=$(\ls stack*fits);for file in $files;do fold $file | head -105 | grep '\|IMGTAG\|EXPOSURE\|FRAMESTA\|FRAMEEND\|IMGCOUNT\|ROTENC\|LINAENC\|LINBENC\|LINCENC\|RA_DTU\|DEC_DTU\|ROLL_DTU\|AZ\|EL\|MROT\|CCDTEMP\|CAMTEMP\|PWSTEMP\|TRTD\|PRESSURE\|VALID0\|USE0\|TX0\|TY0\|CX0\|CY0\|FLUX0\|SIGMAX0\|SIGMAY0\|VALID1\|USE1\|TY1\|CX1\|CY1\|FLUX1\|SIGMAX1\|SIGMAY1\|VALID2\|USE2\|TX2\|TY2\|CX2\|CY2\|FLUX2\|SIGMAX2\|SIGMAY2\|VALID3\|USE3\|TX3\|TY3\|CX3\|CY3\|FLUX3\|SIGMAX3\|SIGMAY3\|VALID4\|USE4\|TX4\|TY4\|CX4\|CY4\|FLUX4\|SIGMAX4\|SIGMAY4\|VALID5\|USE5\|TX5\|TY5\|CX5\|CY5\|FLUX5\|SIGMAX5\|SIGMAY5\|VALID6\|USE6\|TX6\|TY6\|CX6\|CY6\|FLUX6\|SIGMAX6\|SIGMAY6\|VALID7\|USE7\|TX7\|TY7\|CX7\|CY7\|FLUX7\|SIGMAX7\|SIGMAY7\|DATE';done >> fitsheader.txt;"""
    print (function.replace('\\\\','\\'))
    subprocess.call(function.replace('\\\\','\\'), shell=True)
    
    #function="""files=$(ls %s/stack*fits);for file in $files;do echo $file;done"""%(os.path.dirname(path))
    #subprocess.call(function, shell=True
    HeaderCatalog(path=os.path.dirname(path) + '/fitsheader.txt')
