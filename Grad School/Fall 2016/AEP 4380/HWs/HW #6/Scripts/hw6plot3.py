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

data = "task4.dat"
z = columnToList(data, 0)
v = columnToList(data, 1)
Er = columnToList(data, 2)
Ez = columnToList(data, 3)
plot(z,v, label=r'$V$'+' '+r'$(volts)$')
plot(z,Er, label= r'$E_r$'+' '+r'$(V/mm)$')
plot(z, Ez, label= r'$E_z$'+' '+r'$(V/mm)$')
plt.xlim(-10,10)
plt.xlabel(r'$z$' + ' ' + r'$(mm)$')
plt.grid(b=True, which='major')
plt.minorticks_on
plt.legend()
show()