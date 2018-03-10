#!/usr/bin/env python
import pyfits
import glob
import numpy as np
import math
from pylab import *
from matplotlib import pyplot as plt

c=3*10**5

def bias(biasfiles):
    l=[]
    for f in biasfiles:
        m=pyfits.getdata(f,0)
        l.append(m)
    a=np.array(l)
    return np.mean(a,axis=0)
    
def adjustments(spectrum,finalflat,bias):
    biasadjust=np.subtract(spectrum,bias)
    final=np.divide(biasadjust,finalflat)
    return final

def domeflat(flats):
    l=[]
    for f in flats:
        m=pyfits.getdata(f,0)
        l.append(m)
    a=np.array(l)
    median=np.median(a,axis=0)
    mean=np.mean(a,axis=0)
    return np.divide(median,mean)

# open flat files, spectra files, bias files
flats=sorted(glob.glob('Flats/*.fit'))
spectrum1=pyfits.getdata('Spectra/NGC7331_Spectrum_001.fit',0)
spectrum2=pyfits.getdata('Spectra/NGC7331_Spectrum_002.fit',0)
biasfiles=sorted(glob.glob('Bias/*.fit'))

# scale spectrum 2 and add to spectrum 1
adjusteds2=np.divide(spectrum2,2)
adjusteds=np.add(spectrum1,adjusteds2)

# calculate bias
bias=bias(biasfiles)

# calculate flat
finalflat=domeflat(flats)

# adjust spectrum using bias and flats
adjustedspectrum=adjustments(adjusteds,finalflat,bias)

# locate center of the galaxy
means=[]
for i in arange(0,len(adjustedspectrum)):
    row=adjustedspectrum[i]
    mean=np.mean(row)
    means.append(mean)
centerloc=np.argmax(means)

# select box where the Halpha line is
m1=adjustedspectrum[centerloc-80:centerloc+80,90:120]

# bin vertically
m2=[]
for i in np.linspace(0,158,80):
    newrow=m1[i]+m1[i+1]
    m2.append(newrow)
m=[]
for n in np.linspace(0,78,40):
    newrow=m2[int(n)]+m2[int(n+1)]
    m.append(newrow)
    
# calculate conversion factors from pixel number to wavelength
slope=(5461.0-4358.0)/(807.0-1504.0)
n=5461-slope*807

# calculate observed Halpha at center
hloc=np.argmax(m[-1])+90
hcenter=hloc*slope+n

# calculate velocity at center
vcenter=c*(hcenter-6563.0)/6563.0
print(vcenter)
print(hcenter)

# find location (in terms of pixel number) of maximum in each row and create a list of these
hlist=[]
for i in arange(0,len(m)):
    row=m[i]
    halphaloc=np.argmax(row)+90
    hlist.append(halphaloc)

# transform pixel number to wavelength
hlambda=np.add(np.multiply(hlist,slope),n)

# calculate velocity
vlist=[]
for h in hlambda:
    v=c*(h-6563.0)/6563.0
    vlist.append(v)

# subtract galactocentric velocity
newvlist=[]
for v in vlist:
    newv=v-vcenter
    newvlist.append(newv)
    
# plot velocity vs arbitrary radius
radius=list(reversed(arange(0,len(newvlist))))
scatter(radius,np.absolute(newvlist))
plt.grid(True, which='both')
plt.minorticks_on
show()