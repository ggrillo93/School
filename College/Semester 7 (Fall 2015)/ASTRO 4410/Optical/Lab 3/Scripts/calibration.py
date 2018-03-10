#!/usr/bin/env python
import pyfits
import glob
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt

m=pyfits.getdata('Calibration_Hg.fit',0)
center=m[212]
x=arange(0,10)
plot(x,center[1500:1510])
show()
