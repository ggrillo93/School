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

t = columnToList("data1.dat", 0)
d = columnToList("data1.dat", 1)
plot(t, d)
plt.minorticks_on
plt.xlabel("Time step")
plt.ylabel("End-to-end distance")
plt.savefig("tvsdb.png")

E = columnToList("data1.dat", 2)
plt.figure()
plt.plot(t,E)
plt.xlabel("Time steps")
plt.ylabel("Energy")
plt.savefig("tvsEb.png")
plt.show()
