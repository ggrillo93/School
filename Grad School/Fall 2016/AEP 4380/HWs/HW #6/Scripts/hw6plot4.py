#!/usr/bin/env python
from pylab import *

clf()
data = loadtxt("cplotmain.dat", 'float')
r = loadtxt("cplotr.dat")
z = loadtxt("cplotz.dat")
cs = contour(r,z, data, 10)
clabel(cs)
ylim(-2,2)
xlim(2,11)
xlabel(r'$r$' + ' ' + r'$(mm)$')
ylabel(r'$z$' + ' ' + r'$(mm)$')
show()