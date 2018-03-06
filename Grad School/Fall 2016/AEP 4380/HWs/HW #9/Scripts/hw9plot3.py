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

def strToFloat(lista):
	nlista=[]
	for val in lista:
		try:
			newval=float(val)
		except ValueError:
			newval=0
		nlista.append(newval)
	return nlista

x = strToFloat(columnToList("data2.dat", 1))
y = strToFloat(columnToList("data2.dat", 2))
x1 = []
y1 = []
for n in range(136,180):
    x1.append(x[n])
    y1.append(y[n])
plot(x1, y1, marker = 'o', markersize=10)
plt.minorticks_on
plt.xlabel("x")
plt.ylabel("y")
plt.ylim(min(y1)-1, max(y1)+1)
plt.xlim(min(x1)-2, max(x1)+3)
plt.grid()
plt.savefig("t=1e7b.png")
