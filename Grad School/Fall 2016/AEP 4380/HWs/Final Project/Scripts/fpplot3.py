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

T = [1e7,2e7,4e7,8e7]
nu = [4.60, 3.54, 2.72, 2.08]
x = log(np.linspace(1e7, 8e7, 1000))
fit  = np.polyfit(log(T), nu, 1)
plot(log(T), nu, label = r'$Data$')
plot(x, fit[0]*x+fit[1], label = r'$Least-squares\,fit$')
plt.legend(loc=0)
plt.xlabel(r'$\ln{T}\,(K)$')
plt.ylabel(r'$\nu$')
plt.savefig('fig3.png')
show()
