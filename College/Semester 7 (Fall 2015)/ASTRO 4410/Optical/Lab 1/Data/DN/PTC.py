#!/usr/bin/env python
import pyfits
import numpy as np
import glob
from pylab import *

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
    signal = np.mean(b)-188.83
    return signal

def noiselist(files):
    lista2 = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista2.append((noise(m)))
    noiselistlog = np.log10(lista2)
    return noiselistlog

def signallist(files):
    lista = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista.append((signal(m)))
    signallistlog = np.log10(lista)
    return signallistlog

files=sorted(glob.glob('*.fits'))
signal=signallist(files)
noise=noiselist(files)

scatter(signal,noise)
grid(True)
xlabel('log(Signal) = log(S) (DN)')
ylabel('log(Noise) = log(Sigma_S) (DN)')
title('Photon Transfer Curve')
show()
