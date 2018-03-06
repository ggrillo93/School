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

data = "solution.dat"
m = columnToList(data, 0)
r = columnToList(data, 1)
R = []
for i in r:
    R.append(float(i))
rf = []
ma = R[-1]
for i in R:
    rf.append(i/ma)
P = columnToList(data, 2)
rho = columnToList(data, 3)
data2 = "model.dat"
r2 = columnToList(data2, 0)
P2 = columnToList(data2, 3)
rho2 = columnToList(data2, 2)
plot(rf,P, label = 'Grillo Model')
plot(r2, P2, label = 'Model S')
plt.xlabel(r'$r/R$')
plt.ylabel(r'$P \, (dyne/cm^2)$')
plt.xlim(0,1)
plt.legend()
savefig('final1.png')
plt.figure()
plot(rf,rho, label = 'Grillo Model')
plot(r2,rho2, label = 'Model S')
plt.ylabel(r'$\rho \, (g/cm^3)$')
plt.xlabel(r'$r/R$')
plt.xlim(0,1)
plt.legend()
savefig('final2.png')
