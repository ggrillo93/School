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

A = [0.5,1,2]
x = linspace(0,2,1000)
maxsk = [log(70),log(175),log(1000)]
fit = np.polyfit(A,maxsk,1)
plot(A,maxsk, label = r'$Data$')
plot(x,log(linspace(0,1000,1000)), label = r'$Least-squares\,fit$')
plt.legend(loc=0)
show()
