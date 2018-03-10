#!/usr/bin/env python
import pyfits
import numpy as np

def box(matrix,n):
    rows = matrix[-n:]
    result = []
    for row in rows:
        result.append(row[-n:])
    return result

def noise(fits):
    m = pyfits.getdata("fits",0)
    m = m[0]
    box = box(m,20)
    noise = np.std(box)
    return noise

def signal(fits):
    m = pyfits.getdata("fits",0)
    m = m[0]
    box = box(m,20)
    signal = np.mean(box)-1.2123
    return signal
