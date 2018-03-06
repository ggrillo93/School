from numpy import *
from matplotlib import *
from matplotlib import pyplot
from matplotlib.pyplot import hist
from matplotlib.pyplot import *
from pandas import read_csv


def readdata():
   """
   Reads birthdata expressed as day numbers in year
   Data are from A6523 student and their facebook friends
   Based on code from Ross Jennings, 2017 April
   """
   infile1 = 'BirthdaysSheet1.csv'
   phases = read_csv(infile1)['Phase'].values

   #The ['Phase'] part extracts the column with that header, and the .values converts it into a numpy array.

   #Reading the second sheet (with 1156 birthdays from our Facebook friends) was a little more complicated, especially because I wanted to get it in the same format as the data from the first sheet. Here's what I did:

   infile2 = 'BirthdaysSheet2.csv'
   sheet2 = read_csv(infile2)
   phase, total = sheet2['Phase'].values[:-1], sheet2['Total'].values[:-1]

   #The [:-1] removes the last element from each array. For phase, this is a NaN, and for total, it's the total number of birthdays (1156). To change it into a format resembling that of the data from the first sheet, I ended up doing the following:

   phases_all = list()
   for i in range(365):
       for j in range(total[i]): phases_all.append(phase[i])
   phases_all = np.array(phases_all)
   return phases, phases_all

# Main
#if  __name__  == '__main__':

phases, phases_all = readdata()


# fig = figure()
# range = ((0,1.))
# hphall, binsall, patches = hist(phases_all, bins=30, range=range, normed=False, label='large sample')
# hph, bins, patches = hist(phases, bins=30, range=range, normed=False, label='small sample')
# xlabel('Phase (cycles)' )
# ylabel('Counts' )
# title('Coarse histogram')
# legend(loc=0)
# show()
# savefig('bday_data_coarse_histograms.pdf')
#
#
# fig = figure()
# range = ((0,1.))
# hphall, binsall, patches = hist(phases_all, bins=365, range=range, normed=False, label='large sample')
# hph, bins, patches = hist(phases, bins=300, range=range, normed=False, label='small sample')
# xlabel('Phase (cycles)' )
# ylabel('Counts' )
# title('Fine histogram')
# legend(loc=0)
# show()
# savefig('bday_data_fine_histograms.pdf')
#
# fig = figure()
phivec = arange(0., 1., 0.01)

def pdfmodel1(phi):
    """
    PDF model 1 = constant
    Range of phi is (0,1)
    """
    pdfmodel1 = ones(size(phi))
    return pdfmodel1

def pdfmodel2(phi, B, offset):
    """
    PDF model 2 = constant + sinusoid
    Need B < 1.
    Assume range of phi is [0,1]
    Normalization is unity because sine has one-year period
    """
    return  1 + B*sin(2.*pi*(phi-offset))

def pdfmodel3(phi, A, H, B=0.5, C=0.25):
    """
    PDF model 3 = constant + rectangular bump with fixed location and width
    A = bump size
    B = start phase of bump
    C = bump width
    need B+C < 1.
    Assume range of phi is [0,1]
    """
    norm = A + H * C
    pdfvec= zeros(size(phi))
    inds1 = where(phi<=B)
    inds2 = where((B < phi) & (phi <= B+C))
    inds3 = where(B+C < phi)
    pdfvec[inds1] = A
    pdfvec[inds2] = A+H
    pdfvec[inds3] = A
    pdfvec/= norm
    return pdfvec

# plot(phivec, pdfmodel1(phivec), drawstyle='steps', label='model 1')
# plot(phivec, pdfmodel2(phivec, 1., 0.2), drawstyle='steps', label='model 2')
# plot(phivec, pdfmodel3(phivec, 1., 0.3, B=0.4, C=0.4), drawstyle='steps', label='model 3')
# legend(loc=1)
# show()
# savefig('bday_pdfmodels.pdf')
#
#
# # Likelihood functions for phases_all
#
# # Model 1 = constant:
# LL1 = sum(log((pdfmodel1(phases_all))))
# print "Model 1: ln(likelihood) = ", LL1

# Model 2 = constant + sinusoid

Asin = 1.		# held constant since normalization couples A,B,C
Bmin = 1.e-10
Bmax = Asin/2.
Bvec = arange(Bmin, Bmax, 0.01)	# amplitude of sinusoid
Cvec = arange(0., 1., 0.01)	# offset in cycles

# prior for B is 1/B except at B=epsilon where prior is set to next value
priorB = zeros(size(Bvec))
priorB[1:] = (1./ Bvec[1:]) * log(Bvec[-1] / Bvec[0])
priorB[0] = priorB[1]
priorC = ones(size(Cvec)) / Cvec.max()

#ll2vec = zeros((size(Avec), size(Bvec), size(Cvec)))
ll2array= zeros((size(Bvec), size(Cvec)))
lposterior2_unnorm = zeros((Bvec.size, Cvec.size))

for nB, B in enumerate(Bvec):
   lpB = log(priorB[nB])
   for nC, C in enumerate(Cvec):
        lpC = log(priorC[nC])
        LL2 = sum(log((pdfmodel2(phases_all, B, C))))
        ll2array[nB, nC] = LL2
        lposterior2_unnorm[nB, nC] = lpB + lpC + LL2

# location of maximum likelihood in array
nBbest, nCbest = unravel_index(ll2array.argmax(), ll2array.shape)
Bbest = Bvec[nBbest]
Cbest = Cvec[nCbest]
print "Model 2: ",  nBbest, nCbest, "   ",  Bbest, Cbest, " ln(likelihood) = ", ll2array[nBbest, nCbest]

# location of maximum posterior
nBpost, nCpost= unravel_index(lposterior2_unnorm.argmax(), lposterior2_unnorm.shape)
Bpost= Bvec[nBpost]
Cpost= Cvec[nCpost]
print "Model 2: ",  nBpost, nCpost, "   ",  Bpost, Cpost, " ln(posterior_u) = ", lposterior2_unnorm[nBpost, nCpost]

# imshow of likelihood
fig=figure()
extent=(Bvec.min(), Bvec.max(), Cvec.min(), Cvec.max())
imshow(ll2array.T, aspect='auto', origin='lower', extent=extent)
plot(Bbest, Cbest, '+', ms=12)
#imshow(ll2array.T, aspect='auto', origin='lower')
xlabel('B' )
ylabel('C' )
title(r'$\rm Log \ likelihood \ Model \ 2 \ = \ 1 \ + \ B\ \sin(phi-C)$', fontsize=11)
print(len(ll2array.T[0]))
colorbar()
show()

# savefig('bday_imshow_likelihood_model2.pdf')

# # imshow of unnormalized posterior
# fig=figure()
# extent=(Bvec.min(), Bvec.max(), Cvec.min(), Cvec.max())
# imshow(lposterior2_unnorm.T, aspect='auto', origin='lower', extent=extent)
# plot(Bpost, Cpost, '+', ms=12)
# #imshow(ll2array.T, aspect='auto', origin='lower')
# xlabel('B' )
# ylabel('C' )
# title(r'$\rm Unnormalized \ posterior \ Model \ 2 \ = \ 1 \ + \ B\ \sin(phi-C)$', fontsize=11)
# colorbar()
# show()
# savefig('bday_imshow_posterior_unnorm_model2.pdf')
#
#
# # Model 3 = constant + tophat
# # variable parameters = tophat height and location (H, B)
# A = 1		# constant (but gets renormalized in pdfmodel3 function)
# C = 0.25	# tophat width
# dH = 0.01
# dB = 0.01
# epsilon=1.e-10				# to avoid 0.
# Hvec = arange(epsilon, 0.5+dH, dH)		# tophat height
# Bvec = arange(0., (1.-C)+dB, dB)	# tophat leading edge
#
# priorH = (1. - Hvec/Hvec.max()) * 2. / Hvec.max()
# priorH = (1. - Hvec/Hvec.max())**2 * 3. / Hvec.max()
# priorH = (1. - (Hvec/Hvec.max())**2+epsilon) * (1. / (2./3. + epsilon))  /  Hvec.max()
# priorB = ones(size(Bvec)) / Bvec.max()
#
# ll3array = zeros((Hvec.size, Bvec.size))
# lposterior3_unnorm = zeros((Hvec.size, Bvec.size))
# for nH, H in enumerate(Hvec):
#    lpH = log(priorH[nH])
#    for nB, B in enumerate(Bvec):
#       lpB = log(priorB[nB])
#       LL3 = sum(log((pdfmodel3(phases_all, A, H, B=B, C=C))))
#       ll3array[nH, nB] = LL3
#       lposterior3_unnorm[nH, nB] = lpH + lpB + LL3
# # maximum likelihood
# nHbest, nBbest = unravel_index(ll3array.argmax(), ll3array.shape)
# Hbest = Hvec[nHbest]
# Bbest = Bvec[nBbest]
# print "Model 3: ", nHbest, nBbest, "   ", Hbest, Bbest, " ln(likelihood) = ", ll3array[nHbest, nBbest]
#
# # maximum posterior:
# nHpost, nBpost= unravel_index(lposterior3_unnorm.argmax(), lposterior3_unnorm.shape)
# Hpost= Hvec[nHpost]
# Bpost= Bvec[nBpost]
# print "Model 3: ", nHpost, nBpost, "   ", Hpost, Bpost, " ln(posterior_u = ", lposterior3_unnorm[nHpost, nBpost]
#
# # imshow of log likelihood
# fig=figure()
# extent=(Hvec.min(), Hvec.max(), Bvec.min(), Bvec.max())
# imshow(ll3array.T, aspect='auto', origin='lower', extent=extent)
# plot(Hbest, Bbest, '+', ms=12)
# xlabel('H' )
# ylabel('B' )
# title(r'$\rm log \ likelihood \ Model \ 3 \ = \ constant \ + \ tophat \  location =  B, \  height = H$', fontsize=11)
# colorbar()
# show()
# savefig('bday_imshow_likelihood_model3.pdf')
#
# # imshow of unnormalized posterior
# fig=figure()
# extent=(Hvec.min(), Hvec.max(), Bvec.min(), Bvec.max())
# imshow(lposterior3_unnorm.T, aspect='auto', origin='lower', extent=extent)
# plot(Hpost, Bpost, '+', ms=12)
# xlabel('H' )
# ylabel('B' )
# title(r'$\rm Unnormalized \ posterior  \ Model \ 3 \ = \ constant \ + \ tophat \ with \ location =  B \quad\quad height = H$', fontsize=11)
# colorbar()
# show()
# savefig('bday_imshow_posterior_unnorm_model3.pdf')
# raw_input('hit return')
# close('all')
