from numpy import *
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
from astropy.table import Table
from astropy.stats import LombScargle
from scipy import signal
from detect_peaks import *
import uncertainties as unc
from uncertainties import unumpy as unp
import glob

path = '/home/gian/Documents/School/Grad School/Spring 2017/ASTRO 6523/Final Project/Figures/'
def gaussian(height,center,fwhm,xvec):
    e=(-4*np.log(2)**2*(xvec-center)**2)/(fwhm**2)
    g=height*np.exp(e)
    return -g

def whiten(timeflux, minpeak, quarter = '', nit = 5, p = False):
    n = 1
    if p == True:
        f, axarr = plt.subplots(nit-1, sharex = True, figsize=(15,15))
    time = timeflux[0]
    phase = timeflux[1]
    flux = timeflux[2]
    fluxerr = timeflux[3]
    while n < nit:
        flux = flux-mean(flux)
        fluxerr = abs(fluxerr-mean(fluxerr))
        if p == True:
            axarr[n-1].scatter(phase,flux,marker='+',color='b')
        freq, power = LombScargle(phase,flux).autopower()
        plt.figure()
        plt.plot(freq,power)
        plt.title('Lomb-Scargle Periodogram, Q'+quarter)
        plt.ylabel('Power (arb)')
        plt.xlabel('Frequency (1/phase)')
        plt.savefig(path+'periodogram'+quarter+'.png')
        plt.close()
        peak = power.argmax()
        best = freq[peak]
        if best < minpeak:
            print("Peak too small")
            break
        else:
            y = LombScargle(phase,flux).model(phase,best)
            if p == True:
                axarr[n-1].plot(phase,y,color='r')
                axarr[n-1].set_ylim(min(flux)-0.01,max(flux))
            flux = subtract(flux, y)
            n = n + 1
    plt.setp([a.get_xticklabels() for a in axarr[0:-1]], visible=False)
    plt.xlabel('Phase')
    f.text(0.04,0.5,'Normalized flux', va = 'center', rotation='vertical')
    axarr[0].set_title('Lomb-Scargle Whitening' + ', Q'+quarter)
    savefig(path+'white' +quarter+'.png')
    return asarray([time,phase,flux,fluxerr])

def readclean(name, fmt, hd, cadence):
    data = asarray(Table.read(name, format = fmt, hdu = hd))
    time = []
    flux = []
    fluxerr = []
    phase = []
    tinit, tfin = data[0][0], data[-1][0]
    diff = tfin-tinit
    if cadence == 'long':
        for i in data:
            if isnan(i[4]) == False:
                time.append(i[0])
                flux.append(i[4])
                phase.append((i[0]-tinit)/diff)
        return asarray([time,phase,flux,tinit,tfin])
    elif cadence == 'short':
        for n in range(len(data)):
            if isnan(data[n][3]) == False and n < len(data)-1:
                diff2 = abs(data[n+1][3]-data[n][3])
                if diff2 < 50:
                    time.append(data[n][0])
                    flux.append(data[n][3])
                    fluxerr.append(data[n][4])
                    phase.append((data[n][0]-tinit)/diff)
        return asarray([time,phase,flux,fluxerr])
    else:
        print("Cadence not valid")
        return

def scplot(data, q):
    q = str(q)
    fig = plt.figure(figsize = (10,5))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.scatter(data[0],data[2], marker = '+', c = 'b')
    ax1.set_xlabel('Time (BJD days)')
    ax2.plot(data[1], ones(len(data[0])))
    ax2.set_xlabel('Phase')
    ax1.set_ylim(min(data[2])-150, max(data[2])+150)
    ax1.set_ylabel('Flux (electrons/s)')
    plt.title('Short Cadence Light Curve for Kepler-75, (Q'+q+')', y=1.15)
    plt.tight_layout()
    plt.savefig('sclc'+q+'.png')
    return

def findwidth(transit,amp):
    l = []
    for n in range(len(transit[0])-1):
        time = transit[0][n]
        flux1 = transit[1][n]
        flux2 = transit[1][n+1]
        if abs(flux1) >= amp/2. and flux1 < 0 and abs(flux2) >= amp/2. and flux2 < 0:
            l.append(time)
    return l[-1]-l[0]

fits = sorted(glob.glob('*slc.fits'))

Rlist = []
bad = [0,1,2,3,8,11,12,13,14,16]
nitlist = [7,7,0,0,7,8,6,6,7,6,7,0,0,0,0,8,0,7]
diff = []
for n in range(len(fits)):
    if n not in bad:
        d = readclean(fits[n],'fits',1,'short') # read, clean short cadence data quarter
        m = mean(d[2]) # calculate mean of the flux
        d[2] = d[2]/mean(d[2]) # normalize flux
        d[3] = d[3]/mean(d[3]) # normalize error

        w = whiten(d, 0.01, quarter = str(n+1), nit = nitlist[n], p = True) # whiten quarter
        transit = detect_peaks(w[2], valley = True, mph = 0.015, mpd = 100) # roughly detect transit points

        amp = [] # calculate rough amplitude of each transit
        for t in transit:
            amp.append(abs(w[2][t]))
        tamp = mean(amp) # average of rough amplitude, used to create template

        wilist = [] # roughly calculate width of each transit
        for i in range(len(amp)):
            l,u = transit[i]-200,transit[i]+200
            trange = [w[1][l:u], w[2][l:u]]
            width = findwidth(trange, amp[i])
            wilist.append(width)
        wiav = mean(wilist) # average of rough transit width, used to create template

        tvec = linspace(0, 0.01, 1000) # create gaussian template
        gauss = gaussian(tamp, 0.005, wiav, tvec)

        ccf = signal.correlate(w[2], gauss, mode = 'same') # perform matched filtering
        figure()
        plot(w[1], ccf)
        title('CCF (Q'+str(n+1)+')')
        ylabel(r'$CCF$')
        xlabel('Phase')
        ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        tight_layout()
        savefig(path+'ccf'+str(n+1)+'.png') # plot CCF

        if n+1 != 6: # detect ccf peaks
            ploc = detect_peaks(ccf, mph = 1.5e-2, mpd = 1000)
        if n+1 == 6: # special case for Q6
            ploc = detect_peaks(ccf, mph = 2e-2, mpd = 1000)

        uflux = unp.uarray(w[2],w[3])

        avamp = []
        for p in ploc: # better calculation of dip amplitude
            l,u = p-10,p+10 # average among neighbors
            av = mean(uflux[l:u])
            avamp.append(abs(av))

        f, axarr = plt.subplots(len(ploc), sharey = True, figsize = (5,10)) # plot transit + template
        for a in range(len(ploc)):
            p = ploc[a]
            l,u = p-200,p+200
            trange = [w[1][l:u],w[2][l:u]]
            width = findwidth(trange, avamp[a])
            axarr[a].scatter(trange[0],trange[1])
            axarr[a].plot(trange[0],gaussian(unp.nominal_values(avamp[a]),w[1][p],width,trange[0]), c = 'r', linewidth = 3.0)
            axarr[a].set_xlim([min(trange[0]),max(trange[0])])
            axarr[a].set_ylim(-0.025,0.015)
            xlabel('Phase')
            f.text(0.01,0.5,'Normalized Flux', va = 'center', rotation='vertical')
            f.subplots_adjust(left=0.2)
            axarr[0].set_title('Transits' + ', Q'+str(n+1))
            plt.tight_layout
            savefig(path+'transits'+str(n+1)+'.png')

        dettimes = [] # find detection times
        r = 0
        for p in ploc:
            dettimes.append(w[0][p])
        while r < len(dettimes)-1: # calculate interval between transits
            t1 = dettimes[r]
            t2 = dettimes[r+1]
            d = t2-t1
            if d < 10.:
                diff.append(d)
            r = r + 2

        R = unc.ufloat(0.88,0.04)*9.951*mean(avamp)**0.5 # radius of planet according to quarter data
        Rlist.append(R)

period, perr = mean(diff), std(diff)/(len(diff))**0.5
radius = mean(Rlist)
print([period, perr, radius])
