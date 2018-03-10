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
    noise = np.std(b)
    return noise

def signal(m):
    b = box(m,20)
    signal = np.mean(b)-188.23
    return signal

def noiselist(files):
    lista2 = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista2.append((noise(m)))
    print 'noise before log:' + str(lista2)
    noiselistlog = np.log10(lista2)
    return noiselistlog

def signallist(files):
    lista = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista.append((signal(m)))
    print 'signal before log:' + str(lista)
    signallistlog = np.log10(lista)
    return signallistlog

def idealnoise(files):
    lista = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista.append(math.sqrt(signal(m)))
    print 'ideal noise:' + str(lista)
    signallistlog = np.log10(lista)
    return signallistlog

files=glob.glob('*.fits')
print 'files: ' + str(files)
signal2=sorted(signallist(files))
print 'sorted signal after log: ' + str(signal2)
noise=sorted(noiselist(files))
print 'sorted noise after log: ' + str(noise)
ideal=sorted(idealnoise(files))

(m,b)=polyfit(signal2, noise, 1)
h=scipy.stats.linregress(signal2,noise)
print(h)
print(b)
print(m)
yp=polyval([m,b],signal2)
plot(signal2,yp)
scatter(signal2,noise)
grid(True)
xlabel('log(Signal) = log(S) (DN)')
ylabel('log(Noise) = log(Sigma_S) (DN)')
title('Photon Transfer Curve: Region 2')
show()
