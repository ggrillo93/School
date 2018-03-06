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

f="rootplot.dat"
x=columnToList(f,0)
j0y0=columnToList(f,1)
j2y2=columnToList(f,2)
g=[]
for n in range(len(x)):
    gi=float(j0y0[n])-float(j2y2[n])
    g.append(gi)
plot(x,g,color='red')
y=[0,0]
x2=[0,20]
plot(x2,y,ls='dashed',color='black')
plt.xlabel(r'$x$')
plt.ylabel(r'$J_{0}(x)Y_{0}(x) - J_{2}(x)Y_{2}(x)$')
plt.ylim(-0.2,0.2)
plt.minorticks_on
plt.xticks(np.arange(0, 20, 1.0))
plt.grid(b=True, which='major')
plt.show()