#!/usr/bin/env python
import pyfits
import glob
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt

def images(transitfiles):
    l=[]
    for f in transitfiles:
        m=pyfits.getdata(f,0)
        l.append(m)
    return l

def bias(biasfiles):
    l=[]
    for f in biasfiles:
        m=pyfits.getdata(f,0)
        l.append(m)
    a=np.array(l)
    return np.mean(a,axis=0)
    
def adjustments(images,finalflat,bias):
    biasadjust=[]
    finallist=[]
    for f in images:
        m=np.subtract(f,bias)
        biasadjust.append(m)
    for b in biasadjust:
        final=np.divide(b,finalflat)
        finallist.append(final)
    return finallist

def smallmatrix(matrix,yi,yf,xi,xf):
    small=matrix[yi:yf,xi:xf]
    return small

def maximum(matrix):
    a=np.max(matrix)
    for row_i, row in enumerate(matrix):
        for column_i, column in enumerate(row):
            if column == a:
                y = row_i
                x = column_i
                break
    return (x,y,a)

def star(r,m,yi,xi):
    (x, y, a) = r
    top = yi + y - 15
    bottom = yi + y + 15
    left = xi + x - 15
    right = xi + x + 15
    # print("x, y, a", x, y, a)
    # print("mini_m", top, bottom, left, right)
    return m[top:bottom,left:right]

def evolutiontarget(adjustedtransit):
    l=[]
    for f in adjustedtransit:
        small=smallmatrix(f,250,350,1230,1260)
        r=maximum(small)
        s=star(r,f,250,1230)
        brightness=np.mean(s)
        l.append(brightness)
    return l

def evolutionbenchmark(adjustedtransit):
    l=[]
    for f in adjustedtransit:
        small=smallmatrix(f,150,230,980,1015)
        r=maximum(small)
        s=star(r,f,150,980)
        brightness=np.mean(s)
        l.append(brightness)
    return l

def scale(evolutionbenchmark,evolutiontarget):
    b=evolutionbenchmark[175]
    t=evolutiontarget[175]
    quot=b/t
    adjusted=np.divide(evolutionbenchmark,quot)
    return adjusted

def domeflat(flats):
    l=[]
    for f in flats:
        m=pyfits.getdata(f,0)
        l.append(m)
    a=np.array(l)
    median=np.median(a,axis=0)
    mean=np.mean(a,axis=0)
    return np.divide(median,mean)

# open flat files, transit files, bias files
flats=sorted(glob.glob('DomeFlat/*.fit'))
transitfiles=sorted(glob.glob('Transit/*.fit'))
biasfiles=sorted(glob.glob('Bias/*.fit'))


# calculate bias
bias=bias(biasfiles)

# calculate flat
finalflat=domeflat(flats)

# transform transit files into list of matrices
images=images(transitfiles)

# adjust transit image files using bias and flats
adjustedtransit=adjustments(images,finalflat,bias)

# calculate target star brightness through time
targetevolution=evolutiontarget(adjustedtransit)

# calculate benchmark star brightness through time
benchmarkevolution=evolutionbenchmark(adjustedtransit)

# scale benchmark star brightness through time
# substract benchmark star brightness from target star brightness

scaled1=scale(benchmarkevolution,targetevolution)
compensatedevolution=np.subtract(targetevolution,scaled1)
scaled=np.divide(compensatedevolution,targetevolution[175])
scaledplus1=np.add(1,scaled)

# plot
time=range(0,230)
scatter(time,scaledplus1)
plt.grid(True, which='both')
plt.minorticks_on
show()