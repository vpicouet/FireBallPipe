# uncompyle6 version 3.2.3
# Python bytecode 3.6 (3379)
# Decompiled from: Python 3.6.4 |Anaconda, Inc.| (default, Jan 16 2018, 18:10:19) 
# [GCC 7.2.0]
# Embedded file name: /home/dvibert/work/FireBallPipe/read_ds9_ctr.py
# Compiled at: 2018-09-13 21:38:47
# Size of source mod 2**32: 1102 bytes
"""
Created on Thu Sep 13 19:29:40 2018

@author: dvibert
"""
import numpy as np

def readline_skip_comment(f, comment='#'):
    for line in f:
        if line.find(comment) != 0:
            yield line


def read_one_contour(f):
    for line in readline_skip_comment(f):
        if line.find('(') == 0:
            ctr = []
            for line in readline_skip_comment(f):
                if line.find(')') != 0:
                    ctr.append(np.asfarray(line.split()))
                else:
                    break

            yield np.array(ctr)


def read_ds9_ctr(filename):
    with open(filename, 'r') as (f):
        ctrs = [ctr for ctr in read_one_contour(f)]
    return ctrs


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    ctrs = read_ds9_ctr('Calibration/Slits/F1.ctr')
    plt.figure()
    for c in ctrs:
        if c.shape[0] > 10:
            plt.plot(c[:, 0], c[:, 1], ':k')