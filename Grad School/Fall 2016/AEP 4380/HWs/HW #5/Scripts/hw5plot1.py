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

text = "oscillator.dat"
y = columnToList(text, 1)
yprime = columnToList(text, 2)
t = columnToList(text, 0)
plot(t,y)
plt.xlabel(r'$t$')
plt.ylabel(r'$x(t)$')
plt.minorticks_on
plt.grid(b=True, which='major')
plt.figure()
plot(y,yprime)
plt.xlabel(r'$x(t)$')
plt.ylabel(r'$\frac{dx}{dt}$',size=16)
plt.xlim(-1.05, 1.05)
plt.ylim(-1.05, 1.05)
plt.minorticks_on
plt.grid(b=True, which='major')
show()