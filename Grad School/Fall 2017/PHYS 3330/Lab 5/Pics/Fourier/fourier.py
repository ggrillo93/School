import numpy as np
from matplotlib import pyplot as plt
from numpy.fft import *
from scipy import misc
import glob

def fourierread(name):
    image = misc.imread(name, mode = 'L')
    return image

files = glob.glob('*.JPG')
for f in files:
    fig, axarr = plt.subplots(1, 2, figsize = (14, 7))
    im = fourierread(f)
    axarr[0].imshow(im, cmap = 'jet', origin = 'lower', aspect = 'auto')
    axarr[0].set_axis_off()
    fftim = fftshift(fft2(im))
    axarr[1].imshow(np.log10(np.abs(fftim)), cmap = 'jet', origin = 'lower', aspect = 'auto')
    axarr[1].set_axis_off()
    plt.savefig(str(f[:-3]) + 'fourier.jpg')
plt.close('all')
