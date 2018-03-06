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

data=glob.glob("ivsx0.dat")
x=columnToList(data[0],0)
i=columnToList(data[0],1)
plot(x,i)
plt.ylabel(r'$\frac{I}{I_{0}}$',fontsize=20)
plt.xlabel(r'$x_{0}$'+ ' ' +r'$(microns)$', fontsize=16)
plt.title(r'$\frac{I}{I_{0}}$'+ ' ' + r'$vs$' + ' '+ r'$x_{0}$', fontsize=20, y=1.02)
show()
