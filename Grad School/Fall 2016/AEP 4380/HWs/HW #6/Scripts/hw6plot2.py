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

data = "task3.dat"
z = columnToList(data, 0)
v = columnToList(data, 1)
E = columnToList(data, 2)

plot(z,E,label=r'$E_z$' + ' ' +r'$(V/mm)$')
plot(z,v,label=r'$V$'+' '+r'$(volts)$')
plt.xlabel(r'$z$'+' '+r'$(mm)$')
plt.grid(b=True, which='major')
plt.minorticks_on
plt.legend()
show()