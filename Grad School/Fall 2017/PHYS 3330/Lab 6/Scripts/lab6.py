from mpl_toolkits.mplot3d import Axes3D
from scipy.fftpack import fft2, ifft2, fftshift
from matplotlib import pyplot as plt
from scipy.signal import fftconvolve
import numpy as np

pi = np.pi
n = 1.33182
k = 2*pi/632.8e-6

def fresnelsqr(a, d):
    return d/(k*a**2)

def hfunc(upxgrid, upygrid, a, d):
    uFsqr = fresnelsqr(a, d)
    arg = (upxgrid**2 + upygrid**2)/(2.*uFsqr)
    return np.exp(1j*arg)

def drop(a, h, ux, uy, c, p):
    alp = ux**2+uy**2
    if alp > 1:
        return 0
    else:
        A = 0.5*(1-n)*(a/h)
        beta = c*A*np.real((ux+1j*uy)**p)
        return 1 + A*alp + beta

def solvekdi(a, d, h, gsize, c, umax, p):
    xxrange = np.linspace(-umax, umax, gsize)
    yyrange = np.linspace(-umax, umax, gsize)
    uFsqr = fresnelsqr(a, d)
    grid = np.meshgrid(xxrange, yyrange)
    uxgrid, uygrid = grid
    ucomb = np.array(grid).T.reshape(-1,2)
    extent = (-umax, umax, -umax, umax)
    grid2 = np.zeros(gsize*gsize)
    for i in range(len(ucomb)):
        comb = ucomb[i]
        grid2[i] = drop(a, h, comb[0], comb[1], c, p)
    g = grid2.reshape(gsize, gsize)
    g = np.exp(1j*k*h*(1-n)*g)
    h2 = hfunc(uxgrid, uygrid, a, d)
    trans = 2*umax**2*fftconvolve(h2, g, mode = 'same')/gsize**2
    soln = np.abs((trans/(2*pi*uFsqr))**2)
    plt.figure()
    plt.imshow(soln.T, origin = 'lower', extent = extent, aspect = 'auto', cmap = 'jet')
    plt.xlabel(r"$u_x'$", fontsize = 16)
    plt.ylabel(r"$u_y'$", fontsize = 16)
    plt.title(r"$\vec{u}'$" + ' vs ' + r"$G$", fontsize = 16)
    plt.colorbar()
    params = (a, d, h, gsize, c, umax, p)
    #plt.savefig('/home/gian/Documents/School/Grad School/Fall 2017/PHYS 3330/Lab 6/test/' + str(d) + '.png')
    plt.show()
    return

# dvec = np.linspace(30, 60, 40)
# for d in dvec:
#     solvekdi(1.5, d, 1., 2500, 0.5, 3., 5)

solvekdi(1.5, 50, 1., 1500, 0.5, 1., 5)
