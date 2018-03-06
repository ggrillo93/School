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

amin = 1
amax = 10
size = 70
avec = linspace(amin,amax,size)
bvec = avec
priora = zeros(len(avec))
priora[1:] = 1./(avec[1:]*log(amax/amin))
priora[0] = priora[1]
priorb = ones(len(bvec))/max(bvec)

papbvec = []
for pa in priora:
    for pb in priorb:
        papbvec.append(log(pa)+log(pb))
papbvec = asarray(papbvec).reshape(size,size)

abvec = []
for a in avec:
    for b in bvec:
        abvec.append([a,b])

papbvec = asarray(papbvec)
avec = []
bvec = []
likelihood = []
for n in abvec:
    a = n[0]
    b = n[1]
    probabilities = []
    for phi in phases_all:
        prob = betamodel(phi, a, b)+1e-10
        probabilities.append(prob)
    likepoint = sum(log(asarray(probabilities)))
    likelihood.append(likepoint)
    avec.append(a)
    bvec.append(b)

avec = asarray(avec)
bvec = asarray(bvec)
likelihood = asarray(likelihood)
extent = (avec.min(), avec.max(), bvec.min(), bvec.max())
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(avec,bvec,likelihood)

m = argmax(likelihood)
maxlike = [avec[m],bvec[m]]
likelihood = asarray(likelihood).reshape(size,size)
plt.figure()
imshow(likelihood.T, aspect='auto', origin='lower', extent=extent)
colorbar()

plt.figure()
imshow(likelihood.T+papbvec.T, aspect='auto', origin='lower', extent=extent)
colorbar()
print(max((likelihood+papbvec).reshape(size*size,1)))
show()
