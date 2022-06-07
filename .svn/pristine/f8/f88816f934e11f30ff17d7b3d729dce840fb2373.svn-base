#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
Created on Sun Aug 12 09:21:47 2018

@author: Vincent
'''

import os,sys
import glob
import datetime
import pandas as pd
import numpy as np
from astropy.table import Table

date = datetime.datetime.now()

if __name__ == '__main__':
    print(__file__)
    path = sys.argv[1]
    print (path)
    days = glob.glob(path+'/1809*')
    days.sort() 
    print(days)
    total_P = pd.read_csv(days[0] +  '/pressure.csv')[0:0]
    for day in days:
	try:
            df = pd.read_csv(day +  '/pressure.csv')
	except IOError:
            pass
        else:    
            total_P = pd.concat([total_P,df])
        #print(total.size)
    totallog = total_P.copy()
    totallog['pressure[mbar]'] =  np.log10(totallog['pressure[mbar]'])
    total_P.to_csv(path + '/TotalPressure.csv', index=False)
    totallog.to_csv(path + '/TotalPressureLog.csv', index=False)

    total_T = pd.read_csv(days[0] +  '/alltemps.csv')[0:0]
    for day in days:
	try:
            df = pd.read_csv(day +  '/alltemps.csv')
	except IOError:
            pass
        else:    
            total_T = pd.concat([total_T,df])
        #print(total.size)
    total_T.to_csv(path + '/Totalalltemps.csv', index=False)
    a = Table.read(path + '/Totalalltemps.csv')
    new_order = ['time','CuClamp1[C]','Getter[C]','CuClamp2[C]','bot-tank-lks[C]','top-tank-lks[C]','EMCCDBack[C]','Reject[C]','Coldhead[C]','DOBC-PV[C]','Water-Temp-bot[C]','MUX-Temp[C]','Water-Temp-bot-pt100[C]','Water-Temp-top[C]']
    t_new = a[new_order]
    
#    threshold1 = 200
#    threshold2 = -220
#    for header in new_order[1:]:
#        t_new[header][t_new[header]>threshold1] = np.nan
#        t_new[header][t_new[header]<threshold2] = np.nan
    
    
    t_new.write(path + '/Totalalltemps.csv',overwrite=True)
    total_Pow = pd.read_csv(days[0] +  '/power.csv')[0:0]
    for day in days:
	try:
            df = pd.read_csv(day +  '/power.csv')
	except IOError:
            pass
        else:    
            total_Pow = pd.concat([total_Pow,df])
        #print(total.size)
    total_Pow.to_csv(path + '/TotalPower.csv', index=False)
    
    date = datetime.datetime.now()
    print('Succesfully merge: date = [{}]'.format(date))  
    #date = '180811'#sys.argv[1]
    #pathdate = path + date
#    with open(pathdate + '/temps.csv', 'rb') as f:
#        reader = csv.reader(f)
#        for row in reader:
#            print row
#    
#fname = cbook.get_sample_data(pathdate + '/pressure.csv', asfileobj=False)
#plt.plotfile(fname, (0, 5, 6))
        
    #df = pd.read_csv(pathdate + '/pressure.csv')
    #df[1:-1:60]#
#banck.

    total_wp = pd.read_csv(days[0] +  '/waterpress_ok.csv')[0:0]
    for day in days:
	try:
            df = pd.read_csv(day +  '/waterpress_ok.csv')
	except IOError:
            pass
        else:    
            total_wp = pd.concat([total_wp,df])
        #print(total.size)
    total_wp.to_csv(path + '/Totalallwaterpress.csv', index=False)
