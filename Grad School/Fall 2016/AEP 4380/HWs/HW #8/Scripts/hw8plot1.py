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

data = "output1.dat"
data2 = "output3.dat"
data3 = "output7.dat"
t = columnToList(data, 0)
co2 = columnToList(data, 1)
fit = columnToList(data3, 2)
residuals = columnToList(data3, 3)
quadfit = columnToList(data2, 2)
plot(t, co2, label = r'$Original$' + ' ' + r'$data$')
plot(t, fit, label = r'$Least \, squares$' + ' ' + r'$fit$')
plot(t, quadfit, label = r'$Least \, squares$' + ' ' + r'$fit$' + ' ' + r'$(quadratic$' + ' ' + r'$only)$')
plt.ylabel(r'$CO_2$' + ' ' + r'$(ppm)$')
plt.xlabel(r'$Time$' + ' ' + r'$(months)$')
plt.minorticks_on
plt.legend(loc=2, fontsize=10)
plt.xlim(0,650)
plt.figure()
scatter(t, residuals)
plot([0,700],[0,0], ls="dashed", color = 'black')
plt.ylabel(r'$Residual$' + ' ' + r'$(ppm)$')
plt.xlabel(r'$Time$' + ' ' + r'$(months)$')
plt.xlim(0,650)
plt.minorticks_on
plt.figure()
plot(t, co2)
plt.ylabel(r'$CO_2$' + ' ' + r'$(ppm)$')
plt.xlabel(r'$Time$' + ' ' + r'$(months)$')
plt.xlim(0,650)
plt.minorticks_on
show()
