from pylab import *
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = linspace(0,500,512)
y = linspace(0,500,512)

data = loadtxt('t=0.dat', 'float')

hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(x, y)
ha.plot_surface(X, Y, data)
plt.xlabel(r'$x\,(m)$')
plt.ylabel(r'$y\,(m)$')
ha.set_zlabel(r'$Re(\psi)\,(m)$')
plt.savefig("plot13d.png")
plt.xlim(0,500)
plt.ylim(0,500)

data = loadtxt('t=0.1.dat', 'float')

hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(x, y)
ha.plot_surface(X, Y, data)
plt.xlabel(r'$x\,(m)$')
plt.ylabel(r'$y\,(m)$')
ha.set_zlabel(r'$Re(\psi)\,(m)$')
plt.xlim(0,500)
plt.ylim(0,500)
plt.savefig("plot23d.png")

data = loadtxt('t=0.2.dat', 'float')

hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, data)
plt.xlabel(r'$x\,(m)$')
plt.ylabel(r'$y\,(m)$')
ha.set_zlabel(r'$Re(\psi)\,(m)$')
plt.xlim(0,500)
plt.ylim(0,500)
plt.savefig("plot33d.png")

data = loadtxt('t=0.4.dat', 'float')

hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, data)
plt.xlabel(r'$x\,(m)$')
plt.ylabel(r'$y\,(m)$')
ha.set_zlabel(r'$Re(\psi)\,(m)$')
plt.xlim(0,500)
plt.ylim(0,500)
plt.savefig("plot43d.png")
