#!/usr/bin/env python
import pyfits
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

text = "chenleestable.dat"
t = strToFloat(columnToList(text, 0))
x = strToFloat(columnToList(text, 1))
y = strToFloat(columnToList(text, 2))
z = strToFloat(columnToList(text, 3))
# var = [x,y,z]
# # ylabels = [r'$x(t)$', r'$y(t)$', r'$z(t)$']
# for i in range(len(var)):
# 	plt.figure()
# 	plot(t,var[i])
# 	plt.minorticks_on
# 	plt.xlabel(r'$t$')
# 	plt.ylabel(ylabels[i])
# plot(y,z)
# plt.xlabel(r'$y(t)$')
# plt.ylabel(r'$z(t)$')
# plt.figure()
# plot(x,z)
# plt.xlabel(r'$x(t)$')
# plt.ylabel(r'$z(t)$')
fig = figure(5)
ax = fig.gca(projection='3d')
plt.plot(x,y,z,'k')
plt.xlabel(r'$x(t)$')
plt.ylabel(r'$y(t)$')
ax.set_zlabel(r'$z(t)$')
show()