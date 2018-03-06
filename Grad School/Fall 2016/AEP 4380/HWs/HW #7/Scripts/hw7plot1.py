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

text = "solution.dat"
m = columnToList(text, 0)
P = columnToList(text, 1)
r = columnToList(text, 2)
plot(m, P)
plt.figure()
plot(m, r)
show()
