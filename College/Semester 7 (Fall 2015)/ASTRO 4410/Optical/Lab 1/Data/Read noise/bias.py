#!/usr/bin/env python
import glob
import numpy as np
import pyfits

files=glob.glob('*.fits')
result=0
for fits in files:
    data = pyfits.getdata(fits,0)
    noise = np.mean(data)
    result = result + noise
av = result/len(files)
print av
