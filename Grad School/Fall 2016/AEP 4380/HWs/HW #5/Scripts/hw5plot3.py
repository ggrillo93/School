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

text = "solarsystem2.dat"
xearth = columnToList(text, 1)
yearth = columnToList(text, 2)
xmars = columnToList(text, 3)
ymars = columnToList(text, 4)
xjupiter = columnToList(text, 5)
yjupiter = columnToList(text, 6)
xsun = columnToList(text, 7)
ysun = columnToList(text, 8)
plot(xsun[1:], ysun[1:], label = "Sun")
# plot(xearth[1:], yearth[1:], label = "Earth")
# plot(xmars[1:], ymars[1:], label = "Mars")
# plot(xjupiter[1:], yjupiter[1:], label = "Jupiter")
plt.legend()
plt.ylabel("y-coordinate")
plt.xlabel("x-coordinate")
show()