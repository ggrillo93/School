#!/usr/bin/env python
import pyfits
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt
import glob
import os

def columnToList(text,colnumber):
    with open(text) as f:
        lista=[]
        for line in f:
            spl=line.split()
            lista.append(spl[colnumber])
    return lista

data1=sorted(glob.glob("*0.5.dat"))
data2=sorted(glob.glob("*0.05.dat"))

colors=['r','g','orange']
legends=['Backward Diff','Central Diff', 'Forward Diff']

f1=columnToList(data1[0],1)
f2=columnToList(data2[0],1)
x1=columnToList(data1[0],0)
x2=columnToList(data2[0],0)

plt.figure()
plot(x1,f1,color='blue',label='Original Function')
for n in range(len(data1)):
	fp=data1[n]
	der=columnToList(fp,2)
	plot(x1,der,color=colors[n],label=legends[n])
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Numerical differentiation, h=0.5')

plt.figure()
plot(x2,f2,color='blue',label='Original Function')
for n in range(len(data2)):
	fp=data2[n]
	der=columnToList(fp,2)
	plot(x2,der,color=colors[n],label=legends[n])
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Numerical differentiation, h=0.05')

plt.figure()
plot(x2,f2,color='blue',label='Original Function')
sder=glob.glob('secder.dat')[0]
plot(x2,columnToList(sder,2),color='red',label='Second derivative')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Numerical differentiation #2, h=0.05')

show()