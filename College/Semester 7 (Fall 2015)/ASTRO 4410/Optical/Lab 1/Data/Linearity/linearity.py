#!/usr/bin/env python
import pyfits
import numpy as np
import glob
from pylab import *
import scipy.stats
import math

def box(matrix,n):
    data = np.array(matrix[0,:,:])
    small = data[-n-1:-1,-n-1:-1]
    return small

def box2(matrix,n):
    data = np.array(matrix[0,:,:])
    small = data[0:20,0:20]
    return small

def signal(m):
    b = box(m,20)
    signal = np.mean(b)
    return signal

def signal2(m):
    b = box2(m,20)
    signal = np.mean(b)
    return signal

def signallist(files):
    lista = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista.append((signal(m)))
    return lista

def signallist2(files):
    lista = []
    for fits in files:
        m = pyfits.getdata(fits,0)
        lista.append((signal2(m)))
    return lista

files=glob.glob('*.fits')
signal=sorted(signallist(files))
signal2=sorted(signallist2(files))
exptime=[55,60,70,80,100,150,300]

def listofvalues(m,b,lista):
    values = []
    for x in lista:
        y = m*x + b
        values.append(y)
    return values

def maximumdeviation(ideal,actual):
    l1 = np.array(ideal)
    l2 = np.array(actual)
    a = l2 - l1
    return [max(a),math.fabs(min(a)),max(l2)]

def nonlinearity(list):
    nl = 100*((list[0]+list[1])/list[2])
    return nl

#print 'stats for line 1: ' + str(scipy.stats.linregress(exptime,signal))
#print 'stats for line 2: ' + str(scipy.stats.linregress(exptime,signal2))

(m1,b1)=polyfit(exptime, signal2, 1)
#print(b1)
#print(m1)
yp=polyval([m1,b1],exptime)
plot(exptime,yp)
scatter(exptime,signal2)
#print 'Actual data: ' +str(signal2)
#print 'Best fit data: ' +str(listofvalues(m1,b1,exptime))
a = maximumdeviation(listofvalues(m1,b1,exptime), signal2)
print(nonlinearity(a))

(m,b)=polyfit(exptime, signal, 1)
# print(b)
# print(m)
yp=polyval([m,b],exptime)
plot(exptime,yp)
scatter(exptime,signal)
# print 'Actual data: ' +str(signal)
# print 'Best fit data: ' +str(listofvalues(m,b,exptime))
a2 = maximumdeviation(listofvalues(m,b,exptime), signal)
print(nonlinearity(a2))
print((nonlinearity(a2)+nonlinearity(a))/2)
grid(True)
ylabel('Signal (DN)')
xlabel('Exposure time (s)')
title('Linearity curve')
show()