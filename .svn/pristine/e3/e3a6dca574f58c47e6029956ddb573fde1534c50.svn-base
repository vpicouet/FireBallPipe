#! /usr/bin/env python2
from __future__ import print_function
from __future__ import division


import os, sys
import glob
from astropy.io import fits
from astropy.table import Table

def csv2fits(path):
    
    print("converting cvs table from {} to fits".format(path))
    files = glob.glob(path+'*.csv')

    # filter out non csv files
    #files = [f  for f in files if os.path.splitext(f)[1] == '.csv']
    print (files)
    for fname in files:
        # read csv
        t = Table.read(fname, format='csv')
        
        # read originial fits header
        fits_name = fname.replace('_table.csv','.fits')
        header = fits.getheader(fits_name)
        del header[:6]
        
        # convert 2 fits table
        newname = fname.replace('.csv', '.fits')        
        hdu = fits.PrimaryHDU()
        #hdu.header.extend(header)        
        thdu  = fits.table_to_hdu(t)
        thdu.header.extend(header)
        thdulist = fits.HDUList([hdu, thdu])
        thdulist.writeto(newname, overwrite=True )
        print("Table written to {}".format(newname))
 
      
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(prog="csv2fits")
    
    parser.add_argument('-p','--path', default=os.getcwd(), 
                        help='''The path of the files to convert''')
        
    args = parser.parse_args()
    csv2fits(args.path)    


