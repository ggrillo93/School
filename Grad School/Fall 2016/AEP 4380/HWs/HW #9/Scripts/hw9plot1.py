#!/usr/bin/env python
import pyfits
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt

def columnToList(text,colnumber):
    with open(text) as f:
        lista=[]
        for line in f:
            spl=line.split()
            lista.append(spl[colnumber])
    return lista

y1 = columnToList("histogram.dat", 0)
x1 = np.linspace(0,99, 100)
plt.bar(x1,y1,1,color='green')
plt.minorticks_on
plt.ylim(0,1100)
plt.xlim(0,100)
plt.xlabel("Range")
plt.ylabel("N")
plt.savefig("histogram.png")

x2 = columnToList("pairs.dat", 0)
y2 = columnToList("pairs.dat", 1)
plt.figure()
plt.scatter(x2,y2)
plt.xlabel(r'$x_i$')
plt.ylabel(r'$x_{i+1}$')
plt.xlim(-0.1,1.1)
plt.ylim(-0.1,1.1)
plt.savefig("pairs.png")
plt.show()
