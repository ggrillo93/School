#!/usr/bin/env python
import pyfits
import math
import numpy as np
import glob
from pylab import *
import scipy.stats

def box(matrix,n):
    data = np.array(matrix[0,:,:])
    small = data[-n-1:-1,-n-1:-1]
    return small

def noise(m):
    b = box(m,20)
    noise = np.var(b)-1.4697
    return noise

def signal(m):
    b = box(m,20)
    signal = np.mean(b)
    return signal

def noiselist(files):
    lista2 = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista2.append((noise(m)))
    print 'noise before log:' + str(lista2)
    return lista2

def signallist(files):
    lista = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista.append((signal(m)))
    print 'signal before log:' + str(lista)
    return lista

files=glob.glob('*.fits')
print 'files: ' + str(files)
signal2=sorted(signallist(files))
print 'sorted signal after log: ' + str(signal2)
noise=sorted(noiselist(files))
print 'sorted noise after log: ' + str(noise)

(m,b)=polyfit(noise, signal2, 1)
h=scipy.stats.linregress(signal2,noise)
print(h)
print(b)
print(m)
yp=polyval([m,b],noise)
plot(noise,yp)
scatter(noise,signal2)
grid(True)
ylabel('Signal = S (DN)')
xlabel('Noise = (sigma_S^2 - sigma_R^2) (DN)')
title('G Estimate')
show()
