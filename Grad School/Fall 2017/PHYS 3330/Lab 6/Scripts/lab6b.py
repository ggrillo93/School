from scipy.fftpack import fft2, ifft2, fftshift
from matplotlib import pyplot as plt
from scipy.signal import fftconvolve
from matplotlib import gridspec as gs
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

def drop(ux, uy, c, p):
    alp = ux**2 + uy**2
    if alp > 1:
        return 0
    else:
        beta = c*np.real((ux+1j*uy)**p)
        soln = 1 - alp - beta
        if soln > 0:
            return soln
        else:
            return 0

def plotdrop(c, p, umax):
    plt.close()
    nx = np.linspace(-umax, umax, 500)
    ny = np.linspace(-umax, umax, 500)
    ucomb = np.array(np.meshgrid(nx, ny)).T.reshape(-1,2)
    grid = np.zeros(500*500)
    for i in range(len(ucomb)):
        comb = ucomb[i]
        grid[i] = drop(comb[0], comb[1], c, p)
    g = grid.reshape(500, 500)
    extent = [-umax, umax, -umax, umax]
    plt.figure(figsize = (10, 10))
    plt.imshow(g.T, cmap = 'jet', origin= 'lower', aspect = 'auto', extent = extent)
    plt.xlabel(r"$u_x'$", fontsize = 20)
    plt.ylabel(r"$u_y'$", fontsize = 20)
    plt.title(r'$f(\vec{u})$', fontsize = 20)
    plt.colorbar()
    #plt.savefig(str(p) + '.png')
    plt.show()
    return

def solvekdi(a, h, d, c, umax, p, gsize, path = ''):
    uFsqr = fresnelsqr(a, d)
    xxrange = np.linspace(-umax, umax, gsize)
    yyrange = np.linspace(-umax, umax, gsize)
    grid = np.meshgrid(xxrange, yyrange)
    uxgrid, uygrid = grid
    ucomb = np.array(grid).T.reshape(-1,2)
    extent = (-umax, umax, -umax, umax)
    grid2 = np.zeros(gsize*gsize)
    for i in range(len(ucomb)):
        comb = ucomb[i]
        grid2[i] = drop(comb[0], comb[1], c, p)
    g = grid2.reshape(gsize, gsize)
    g2 = np.exp(1j*k*h*(1-n)*g)
    h2 = hfunc(uxgrid, uygrid, a, d)
    field = 2*umax**2*fftconvolve(h2, g2, mode = 'same')/gsize**2
    soln = np.abs((field/(2*pi*uFsqr))**2)
    fig = plt.figure(figsize = (10, 10))
    plt.imshow(soln.T, origin = 'lower', extent = extent, aspect = 'auto', cmap = 'jet')
    plt.xlabel(r"$u_x'$", fontsize = 20)
    plt.ylabel(r"$u_y'$", fontsize = 20)
    plt.title(r"$\left|\epsilon(\vec{u}')\right|^2$", fontsize = 20)
    plt.colorbar()
    plt.show()
    plt.close()
    return
plotdrop(0.2, 8, 3.5)
plotdrop(0.2, 15, 3.5)
path = '/home/gian/Documents/School/Grad School/Fall 2017/PHYS 3330/Lab 6/Plots/'
solvekdi(1.5, 0.1, 50, 0.3, 3., 8., 2500, path = path)
