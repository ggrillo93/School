#!/usr/bin/env python
import pyfits
import numpy as np
import glob
from pylab import *

def box(matrix,n):
    rows = matrix[-n:]
    result = []
    for row in rows:
        result.append(row[-n:])
    return result

def noise(m):
    m = m[0]
    b = box(m,20)
    noise = np.var(b)-1.4697
    return noise

def signal(m):
    m = m[0]
    b = box(m,20)
    signal = np.mean(b)
    return signal

def noiselist(files):
    noiselist = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        noiselist.append(noise(m))
    return noiselist

def signallist(files):
    signallist = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        signallist.append(signal(m))
    return signallist