import numpy as npi
from numpy import *
from pylab import *
from matplotlib import pyplot as plt


def crread(infile='chart_recorder.57331.7872.out'):
   """
   Reads in a chart_recorder_rspec file and plots the time series.
   Note data file may have to be edited to remove spurious lines.
   """
   x,y = loadtxt(infile, usecols=(0,1), skiprows=0, unpack=True, dtype='float')
   
   return x,y

#system temperature and calibration temperature
data=crread(infile='chart_recorder.57331.7872.out')
amplitude=data[1]
diode1=amplitude[114:154]
diode2=amplitude[626:661]
diode=npi.concatenate((diode1,diode2))
V2=npi.mean(diode2)
print('V2 = ' + str(V2))
absorber1=amplitude[464:502]
absorber2=amplitude[799:840]
absorber=npi.concatenate((absorber1,absorber2))
V3=npi.mean(absorber)
print("V3 = " + str(V3))
linearity1=amplitude[287:311]
V4=npi.mean(linearity1)
print("V4 (13 db) = " + str(V4))
linearity2=amplitude[195:241]
V5=npi.mean(linearity2)
print("V5 (7 db) = " + str(V5))
baseline1=amplitude[0:114]
baseline2=amplitude[159:194]
baseline3=amplitude[242:286]
baseline4=amplitude[312:459]
baseline5=amplitude[512:625]
baseline6=amplitude[662:787]
baseline7=amplitude[849:899]
baseline=npi.concatenate((baseline1,baseline2,baseline3,baseline4,baseline5,baseline6,baseline7))
V1=npi.mean(baseline)
print("V1= " + str(V1))
g=(V3-V1)/300.0
print('g = '+ str(g))
Tsys=V1/g
print('Tsys = '+ str(Tsys))
Tcal=(V2-g*Tsys)/g
print('Tcal = ' + str(Tcal))
print('Linearity 7 db = ' + str(V5/(g*Tsys)))
print('Linearity 13 db = ' + str(V4/(g*Tsys)))

#temperature of the sun
data=crread(infile='chart_recorder.57331.8031.out')
amplitude=data[1]
print('Non-adjusted peak = '+str(npi.amax(amplitude)))
adjustedamplitude=npi.subtract(amplitude,V1)
peak=npi.amax(adjustedamplitude)
hm=peak/2
print('Half-maximum = '+str(hm))
count=0
for n in adjustedamplitude:
   if n<hm:
      count=count+1
   else:
      break
# print(count)
# print(adjustedamplitude[count])
lista=adjustedamplitude[count:]
count2=0
for n in lista:
   if n>hm:
      count2=count2+1
   else:
      break
# print(count2+count)
# print(lista[count2])
bw=1012*npi.cos(0.27558749)*0.261799/3600
print('Beam width = '+ str(bw))
sabw=npi.pi*(bw/2)**2
print('Peak = '+str(peak))
sasun=2*npi.pi*(1-npi.cos(0.00872665))
Tsun=(peak/g)
print(g)
print('Tsun = '+ str(Tsun))
# plot(data[0],adjustedamplitude)
# plt.grid(True, which='both')
# plt.minorticks_on
# show()