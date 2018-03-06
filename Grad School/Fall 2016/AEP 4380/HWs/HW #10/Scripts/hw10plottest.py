from pylab import *

x = linspace(0,500,512)
y = loadtxt('2ndpsitest.dat', 'float')
plot(x,y)
show()
