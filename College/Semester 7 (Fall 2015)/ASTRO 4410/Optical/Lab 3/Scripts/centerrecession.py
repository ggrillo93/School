#!/usr/bin/env python
import pyfits
import glob
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt

galaxy=pyfits.getdata('Spectra/NGC7331_Spectrum_002.fit',0)
m=galaxy[117:317,0:200]
# find center
means=[]
for row in m:
    mean=np.mean(row)
    means.append(mean)
print(np.argmax(means))
center=m[184]
slope=(5461.0-4358.0)/(807.0-1504.0)
n=5461-slope*807
# group=1
# averaged=center.reshape(-1, group).mean(axis=1)
# x=np.add(np.multiply(arange(0,len(averaged)),slope),n)
x=np.add(np.multiply(arange(0,len(center)),slope),n)
print(x[np.argmax(center)])
plot(x,center)
show()