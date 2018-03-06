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

text = "wvsiterations.dat"
w = columnToList(text, 0)
n = columnToList(text, 1)
plot(w,n)
plt.xlabel(r'$\omega$')
plt.ylabel(r'$N$')
plt.xlim(1,2)
plt.grid(b=True, which='major')
plt.minorticks_on
show()