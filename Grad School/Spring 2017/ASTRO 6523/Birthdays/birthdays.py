from numpy import *
from pandas import read_csv
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
from scipy import special

sheet2 = read_csv('Birthdays - Sheet2.csv')
phase, total = sheet2['Phase'].values[:-1], sheet2['Total'].values[:-1]

phases_all = list()
for i in range(365):
    for j in range(total[i]): phases_all.append(phase[i])
phases_all = np.array(phases_all)

def sinemodel(phi, A, offset):
    return A*sin(2.*pi*(phi-offset))+1

def betamodel(x,a,b):
    return special.gamma(a+b)/(special.gamma(a)*special.gamma(b)) * x**(a-1) * (1-x)**(b-1)

Amin = 1.e-10
Amax = 0.5
size = 50
Avec = linspace(Amin,Amax,size)
offvec = linspace(0.,1.,size)

priorA = zeros(len(Avec))
priorA[1:] = 1./(Avec[1:]*log(Amax/Amin))
priorA[0] = priorA[1]
prioros = ones(len(offvec))/max(offvec)

papbvec = []
for pa in priorA:
    for pb in prioros:
        papbvec.append(log(pa)+log(pb))
papbvec = asarray(papbvec).reshape(size,size)

Aos = []
for A in Avec:
    for os in offvec:
        Aos.append([A,os])

Avec = []
osvec = []
likelihood = []
for nAos in Aos:
    A = nAos[0]
    os = nAos[1]
    probabilities = []
    for phi in phases_all:
        prob = sinemodel(phi, A, os)
        probabilities.append(prob)
    likepoint = sum(log(asarray(probabilities)))
    likelihood.append(likepoint)
    osvec.append(os)
    Avec.append(A)

Avec = asarray(Avec)
osvec = asarray(osvec)
likelihood = asarray(likelihood)
extent = (Avec.min(), Avec.max(), osvec.min(), osvec.max())
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Avec,osvec,likelihood)

likelihood = asarray(likelihood).reshape(size,size)
plt.figure()
imshow(likelihood.T, aspect='auto', origin='lower', extent=extent)
colorbar()

plt.figure()
imshow(likelihood.T + papbvec.T, aspect='auto', origin='lower', extent=extent)
colorbar()
show()
