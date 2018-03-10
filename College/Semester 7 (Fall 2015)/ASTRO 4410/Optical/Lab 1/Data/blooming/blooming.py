#!/usr/bin/env python

import pyfits
import numpy as np
import glob
from pylab import *
import scipy.stats
import math

time = [1,2,4,8]
maximum = [3862,6759,8929,9606]
plot(time,maximum)
scatter(time, maximum)
grid(True)
ylabel('Maximum signal (DN)')
xlabel('Exposure time (s)')
title('Blooming')
show()