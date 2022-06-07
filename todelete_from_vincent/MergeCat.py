#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 22:18:58 2018

@author: Vincent
"""
from __future__ import print_function
from astropy.table import Table, vstack
import os,sys

def mergeCatThroughfocus(paths):
    t = Table.read(paths[0])
    for file in paths[1:]:
        print(file)
        t0 = Table.read(file)
        t = vstack((t,t0))
    t.write(os.path.dirname(file) + '/MergedCatalog.csv',overwrite=True)
    return t    

mergeCatThroughfocus(paths=sys.argv[1:])
 
