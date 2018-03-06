#!/usr/bin/env python
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

data = 'lanemden.dat'
xi = columnToList(data,2)
theta = columnToList(data, 0)
plot(xi,theta)
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\theta$')
plt.savefig('fig4.png')
show()
