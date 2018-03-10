# 2013 Dec 3
#
import numpy as npi
import subprocess
from numpy import *
from pylab import *
import os
from os import getpid
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy import *
from scipy import stats
from scipy import integrate
from scipy.interpolate import interp1d
from matplotlib import pyplot
#from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
     #grid, savefig, show
from matplotlib.patches import FancyArrowPatch
from datetime import *
from scipy import constants as spconstants
#from cosmolopy import *         # includes constants, distances, etc.
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import pyplot
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
     grid, savefig, show, close
import subprocess
from subprocess import Popen, call, PIPE
import string
#from astro_jmc import *

def Process(infile=''):
    """ 
    Processes the spectral data in the input file by:
	1.  averaging the spectra
	2.  fitting a polynomial baseline and subtracting/dividing by it
    3.  outputting the baselined spectrum to a file 
    """

    if not infile:
       print "Enter the file name of your spectral-line data:"
       print "E.g. spectra.AzEl.180.47.6.56261.9133.out"
       infile = raw_input('Enter and hit return: ')
       print "Input file = ", infile

    outfile = infile.split('.out')[0] + '_calibrated_spectrum.txt'
    plotfile = infile.split('.out')[0] + '_calibrated_spectrum.pdf'
    fout = open(outfile, 'w')

    # read the data file
    frequencies, spectra = read_spectra(infile=infile)
    freqselect, specselect, specresid =  process_spectra(frequencies, spectra)
    
    print >> fout, "Frequency   Spectrum   FinalSpec "
    for n in range(size(freqselect)):
        print >> fout, '%8.4f %10.0f %10.5f'%(frequencies[n], specselect[n], specresid[n])

    print "Spectrum plot is saved as a PDF file: ", plotfile
    print "Output spectrum is saved in text file: ", outfile
    savefig(plotfile)
    fout.close()
 

#def get_rtel_info():
    #RtelLatitude, RtelLongitude, RtelAltitude = a4410_rtel()
    #temperature = 10.           # Celcius
    #pressure = 1000.            # millibars
    #Rtel = make_observer_on_surface(\
        #RtelLatitude, RtelLongitude, RtelAltitude, temperature, pressure)
    #return Rtel

def vlsr(l, b):
    """
    Returns the Sun's peculiar motion in the direction given by l, b.

    l,b in degrees
    returned velocity is in km/s
    """
    vrsun_constant = 19.5		# km/s
    lsun_deg = 57. 
    bsun_deg = 23.
    vrsun = vrsun_constant * (\
            cos(deg2rad(b))*cos(deg2rad(bsun_deg)) * \
            cos(deg2rad(l-lsun_deg)) + sin(deg2rad(b))*sin(deg2rad(bsun_deg))\
            )
    return vrsun

def read_header(infile='spectra.AzEl.180.47.6.56261.9133.out'):
    """
    infile = 'spectra.lb.0.0.56246.9239.out'
    Read header lines from spectrum file written by spectrum_analyzer6.
    """
    #statinfo = os.stat(infile)
    fin = open(infile, 'r')
    text = open(infile).read()
    nheaderlines=12
    #for line in text.splitlines():
        #if line[0] == 'S' or line[0] == '0':
    count = 0
    while count < nheaderlines:
        line = fin.readline()
        if count == 3: BWch = float(line.split()[0]) 
        if count == 5: Dt = float(line.split()[0]) 
        if count == 7: fcenter = float(line.split()[0]) 
        if count == 8: BW = float(line.split()[0]) 
        if count == 9: Nch = int(line.split()[0]) 
        if count == 10: f1 = float(line.split()[0]) 
        if count == 11: f2 = float(line.split()[0]) 
        count += 1
    nspectra = string.count(text, "\n") / (Nch+1)
    fin.close()
    print BWch, Dt, fcenter, BW, Nch, f1, f2
    return BW, nspectra, Nch, BWch, Dt, fcenter, f1, f2

def read_spectra(infile='spectra.AzEl.180.47.6.56261.9133.out', nspectra_test=100000):
    """
    Reads the spectrum file written by spectrum_analyzer6 and outputs a 2D array
    of spectra and a vector of frequencies. 
    """
    fout = open('read_spectra_metadata.out', 'w')
    #Rtel = get_rtel_info()		# structure needed by NOVAS s/w
    nheaderlines=12

    BW, nspectra, Nch, BWch, Dt, fcenter, f1, f2 = read_header(infile=infile)

    if nspectra_test < nspectra:
        nspectra = nspectra_test
    print "reading ", nspectra, " spectra"

    fin = open(infile, 'r')

    # skip header lines:
    for n in range(nheaderlines): line = fin.readline()

    print >> fout, ' Julian Date    RArad    DECrad    ldeg      bdeg'
    # now get spectra
    spectra = zeros(shape=(Nch, nspectra))
    frequencies = zeros(Nch)
    vradvec=[]
    for ns in range(nspectra):
        # check that there isn't a blank spectrum by checking that
        # next line isn't another spectrum header line
        pre_pos = fin.tell()
        nblanks = 0
        nSlines = 0
        #print pre_pos
        while fin.readline().split()[0] == 'Spectrum:':
            if nSlines > 0: 
               pre_pos = last_pos
               nblanks += 1
            nSlines += 1
            last_pos = fin.tell()
            #print pre_pos, last_pos
        fin.seek(pre_pos)		# back up
        line = fin.readline()
        spectrum_number = int(line.split()[1])
        MJD = float(line.split()[3])
        RArad = float(line.split()[7])
        DECrad = float(line.split()[9])
        ldeg = float(line.split()[19])
        bdeg = float(line.split()[21])

        print >> fout, '%f %f %f %f %f'%(MJD, RArad, DECrad, ldeg, bdeg)

        if nblanks > 0: \
          print "Spectrum: ", spectrum_number, pre_pos, "Nblanks = ",   nblanks
        #print "S: ", spectrum_number, pre_pos
        #fin.seek(last_pos)
        #line = fin.readline()
        spectrum = zeros(Nch)
        for nch in range(Nch):
            line = fin.readline()
            nch_in = int(line.split()[0])
            spectrum[nch] = float(line.split()[1])
            frequencies[nch] = float(line.split()[2])
            spectra[nch, ns] = spectrum[nch]
    return frequencies, spectra

def plot_spectrum(frequencies, spectra):
    """
    Plots a single spectrum that is an average of any multiple 
    spectra in the input array, spectra.
    """
    plotfile1='spectrum_plot.pdf'

    nspectra = shape(spectra)[1]
    avespectrum = mean(spectra, axis=1) 
    plot(frequencies, avespectrum)
    xlabel(r'$\rm Frequency\,\, (MHz)$', fontsize=14)
    ylabel(r'$\rm Spectrum$', fontsize=14)
    title(r'$\rm Average\,\,of\,\,%d \,\,spectra$'%(nspectra))

    savefig(plotfile1)
    print "spectrum is saved as a PDF file: ", plotfile1
    #raw_input('hit return when done viewing spectrum on screen')
    #close()

def process_spectra(frequencies, spectra):
    """
    Plots a single spectrum that is an average of any multiple
    spectra in the input array, spectra.  Does polynomial fitting
    and outputs a baseline subtracted and divided spectrum. 
    """
    #plotfile1='processed_spectrum.pdf'


    subplots_adjust(hspace=0.4)
 
    subplot(311)			# complete spectrum
    nspectra = shape(spectra)[1]
    avespectrum = mean(spectra, axis=1)
    plot(frequencies, avespectrum, 'k-', lw=2)
    xlabel(r'$\rm Frequency\,\, (MHz)$', fontsize=14)
    ylabel(r'$\rm Raw\,\, Spectrum$', fontsize=13)
    title(r'$\rm Average\,\,of\,\,%d \,\,spectra$'%(nspectra))

    fzoom1 = 1418.0
    fzoom2 = 1422.0

    subplot(312)			# zoomed-in spectrum
    indices = where((fzoom1 <= frequencies) & (frequencies <= fzoom2))
    freqselect = frequencies[indices]
    specselect = avespectrum[indices]
    smax = specselect.max()
    plot(freqselect, specselect/smax, 'k-', lw=2)
    xticks((1418, 1419, 1420, 1421, 1422), (r'$1418$', r'$1419$', r'$1420$', r'$1421$', r'$1422$'))
    annotate(r'$\rm Zoomed\,\,spectrum$', xy=(0.1,0.75), xycoords='axes fraction', ha='left', fontsize=13)
    ylabel(r'$\rm Spectrum / Maximum$', fontsize=13)

    polydeg=5
    nranges=3
    franges, ffit, sfit, params2, poly2 = choose_frequencies(freqselect, specselect, polydeg=polydeg, nranges=nranges)
    plot(freqselect, poly2/smax, 'r-', lw=2)
    xlabel(r'$\rm Frequency\,\, (MHz)$', fontsize=14)

    subplot(313)			# residual spectrum
    specresid = (specselect - poly2)/poly2
    plot(freqselect, specresid, 'k-', lw=2)    
    xticks((1418, 1419, 1420, 1421, 1422), (r'$1418$', r'$1419$', r'$1420$', r'$1421$', r'$1422$'))
    xlabel(r'$\rm Frequency\,\, (MHz)$', fontsize=14)
    ylabel(r'$\rm T(\nu) / T_{\rm sys}$', fontsize=14)

    #print "spectrum is saved as a PDF file: ", plotfile1
    #raw_input('hit return when done viewing spectrum on screen')
    #close()
    #savefig(plotfile1)

    return freqselect, specselect, specresid 
   
def choose_time_range_for_spectrum_plot():

    happy=0
    while not(happy):
       print "  Select time range for plotting a time slice"
       print "  Use cursor to select 2 points per, earlier time first"
       print "  First: click on plot window, then click on two time points in increasing time order"
       xxx = ginput(2, timeout=0)
       t1 = xxx[0][1]    
       t2 = xxx[1][1]    
       happy = 1
       print t1, t2

    return t1, t2

def choose_frequencies(frequencies, spectrum, polydeg=3, nranges=4):
    """
    Select frequency ranges and does least-squares fit of polynomial
    Degree of polynomial = polydeg
    Number of frequency ranges to select = nranges

    Assumes spectrum is displayed on screen. 
    """

    # choose frequencies to use for fitting
    a,b,c,d = axis() 

    franges = zeros(shape=(nranges,2))		# pairs of frequencies
    happy=0
    while not(happy):
       print "  Select ", nranges,  " frequency ranges to use in polynomial fit"
       print "  Use cursor to select 2 points per range, left-hand point first"
       print "  Frequency ranges should not overlap"
       print "  Avoid obvious spectral lines"
       print "  First: click on plot window at top:"
       print "  Then select ranges from the plot with label 'Zoomed Spectrum'."
       for nr in range(nranges):
         print "select frequency range #", nr, " click at two frequencies, low frequency first"
         xxx = ginput(2, timeout=0)
         franges[nr,0] = xxx[0][0]    
         franges[nr,1] = xxx[1][0]    
         plot(franges[nr], c+0.1*(d-c)*ones(2), 'r-', lw=2)
       happy = 1
       print franges
      
    ffit=[]
    sfit=[]
    ffit=array(ffit)
    sfit=array(sfit)
    for nch, f in enumerate(frequencies):
        for nr in range(nranges):
            #print "   ", franges[nr,0],  franges[nr,1]
            if franges[nr,0] <= f and f <= franges[nr,1]:
               ffit = append(ffit, f) 
               sfit = append(sfit, spectrum[nch])

    #plot(ffit, sfit, 'go', markersize=3)
    params2 = polyfit(ffit, sfit, polydeg)
    poly2 = polyval(params2, frequencies)
    #plot(frequencies, poly2, 'k-', lw=3) 

    return franges, ffit, sfit, params2, poly2 


