#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 19:11:41 2018

@author: dvibert
"""
import sys
from os import path

filename = '/home/dvibert/Téléchargements/waterpress.csv'

def bytes_from_file_skip_zeros(filename, chunksize=8192):
    with open(filename, "rb") as f:
        while True:
            chunk = f.read(chunksize)
            if chunk:
                for b in chunk:
                    b = bytes([b])
                    if b != b'\x00':
                        yield b
            else:
                break

def convert(filename):
    byte_str = b''
    for b in bytes_from_file_skip_zeros(filename):
        byte_str += b


    filenameout, ext = path.splitext(filename)
    filenameout += '_ok'+ ext


    with open(filenameout,'w') as f:
        f.write(byte_str.decode())
    
if __name__ == '__main__':
    if sys.version_info.major == 3:
        print('Python 3, running converter pressure script')
        convert(sys.argv[1])    
    else:
        print('Not in python 3, I can not run this script')
        pass
    