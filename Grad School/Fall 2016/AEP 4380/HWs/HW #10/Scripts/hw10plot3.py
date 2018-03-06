from pylab import *

pixa = loadtxt('t=0.1.dat', 'float')
img = imshow(pixa, aspect='equal', extent=(0,500,0,500))
img.set_cmap('gray')
colorbar()
xlabel('x (in cm)')
ylabel('y (in cm)')
title('wave at t = 0.1 sec')
show()
