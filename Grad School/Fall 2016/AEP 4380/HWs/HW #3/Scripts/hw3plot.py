#!/usr/bin/env python
import pyfits
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt
import glob
import os

t = np.linspace(-10,10,1000)
f = t*np.sin(t)/2.0
plot(t,f)
show()
