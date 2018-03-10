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

files=glob.glob('*.fits')
signal=sorted(signallist(files))
noise=sorted(noiselist(files))

(m,b)=polyfit(signal, noise,1)
print(b)
print(m)
yp=polyval([m,b],signal)
plot(signal,yp)
scatter(signal,noise)
grid(True)
xlabel('Signal S(DN)')
ylabel('Noise (sigma_S^2-sigma_R^2)')
title('Determination of G')
show()
