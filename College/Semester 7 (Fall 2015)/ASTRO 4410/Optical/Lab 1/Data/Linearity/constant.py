#!/usr/bin/env python
import pyfits
import numpy as np
import glob
from pylab import *

def box(matrix,n):
    data = np.array(matrix[0,:,:])
    small = data[-n-1:-1,-n-1:-1]
    return small

def box2(matrix,n):
    data = np.array(matrix[0,:,:])
    small = data[0:n,0:n]
    return small

def signal2(m):
    b = box2(m,20)
    signal = np.mean(b)
    return signal

def signal(m):
    b = box(m,20)
    signal = np.mean(b)
    return signal

def signallist(files):
    lista = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista.append((signal(m)))
    signallistlog = lista
    return signallistlog

def signallist2(files):
    lista = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista.append((signal2(m)))
    signallistlog = lista
    return signallistlog

files=glob.glob('*.fits')
signal=sorted(signallist(files))
signal2=sorted(signallist2(files))
exptime1=[5,10,20,30,40,45,50,55,60,70,80,100,150,300]

def exppowerlinear(listactm):
    lista2ctm = []
    for i in listactm:
        a = i**0.8972
        lista2ctm.append(a)
    return lista2ctm

exptime = exppowerlinear(exptime1)

    
(m1,b1)=polyfit(signal2, signal, 1)
print(b1)
print(m1)
yp=polyval([m1,b1],signal2)
plot(signal2,yp)
scatter(signal2,signal)

grid(True)
ylabel('log(Signal) = log(S) (DN)')
xlabel('log(Exposure time) (s)')
title('Linearity curve')
show()