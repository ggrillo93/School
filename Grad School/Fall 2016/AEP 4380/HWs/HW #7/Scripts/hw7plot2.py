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

data = "task2.dat"
t = strToFloat(columnToList(data, 0))
x = strToFloat(columnToList(data,1))
psi = strToFloat(columnToList(data, 2))

psi1 = []
psi2 = []
psi3 = []
psi4 = []
psi5 = []
xlist = []
for n in range(len(t)):
    if t[n] == 1e-14:
        psi1.append(psi[n])
        xlist.append(x[n])
    if t[n] == 2e-14:
        psi2.append(psi[n])
    if t[n] == 3e-14:
        psi3.append(psi[n])
    if t[n] == 4e-14:
        psi4.append(psi[n])
    if t[n] == 5e-14:
        psi5.append(psi[n])

l = [psi1,psi2,psi3,psi4,psi5]
titles = ['t = 1e-14 sec', 't = 2e-14 sec', 't = 3e-14 sec', 't = 4e-14 sec', 't = 5e-14 sec']
labels = [r'$x$' + ' ' + r'$(\dot{A})$', r'$|\psi|^2$']
for i in l:
    s = 0
    for n in i:
        s = s+n
    print(s)
#    plt.figure()
#    plot(xlist,l[n])
#    plt.title(titles[n])
#    plt.xlabel(labels[0])
#    plt.ylabel(labels[1])
#    plt.minorticks_on
#    plt.savefig("figure"+str(n+5)+".png", format="png")

