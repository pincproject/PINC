import h5py
import pylab as plt
import numpy as np
from itertools import izip as zip, count

import sys
sys.path.append('../framework')
from pincClass import *

ini = "../../local.ini"
pinc = PINC(iniPath = ini)

# ANALYTICAL TRAJECTORY

timeSteps = np.array([2**(-x) for x in range(3,16)])
trueSizes = np.array([2**x    for x in range(3,16)])

trueSizes[:] = 10
# timeSteps[:] = 0.0001

stepSizes = 2.0/trueSizes


pinc.nTimeSteps = 100
omega = 1

errorMax = np.zeros(timeSteps.shape)
error2 = np.zeros(timeSteps.shape)

for (i,timeStep,trueSize,stepSize) in zip(count(),timeSteps,trueSizes,stepSizes):

	pinc.timeStep = timeStep
	pinc.trueSize = np.array([trueSize])

	n = np.linspace(0,pinc.timeStep*pinc.nTimeSteps,pinc.nTimeSteps+1)
	# posAna = pinc.trueSize[0]*(0.5-0.25*np.cos(omega*n))
	posAna = pinc.trueSize[0]*(0.5+0.25*np.sin(omega*n))

	pinc.runCommand("rm *.h5")
	pinc.puErrorScaling()

	# EXTRACTING ENERGY
	hist = h5py.File('test_history.xy.h5','r')
	pot = hist['/energy/potential/total']
	kin = hist['/energy/kinetic/total']

	kin = kin[:,1];
	pot = pot[:,1];
	tot = pot+kin;

	# EXTRACTING TRAJECTORY OF SINGLE PARTICLE
	pop = h5py.File('test_pop.pop.h5','r')
	N = len(pop['/pos/specie 0/'])
	pos = np.zeros(N)
	for j in xrange(N):
		pos[j]=pop['/pos/specie 0/n=%.1f'%j][:]

	# PLOTTING TRAJECTORIES

	if(i==-1):
		plt.plot(pos,label='Numerical')
		plt.plot(posAna,label='Analytical')
		plt.legend()
		plt.show()
		plt.legend(loc='lower left')

	diff = pos-posAna
	errorMax[i] = np.max(abs(diff))
	error2[i] = np.sqrt(np.sum(diff**2)/(pinc.nTimeSteps+1))

print errorMax

order = np.log(errorMax[-1]/errorMax[-2])/np.log(timeSteps[-1]/timeSteps[-2])
print order

plt.figure()
plt.loglog(timeSteps,errorMax,'o-')
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='black', linestyle='--', alpha=0.5)
plt.minorticks_on()
plt.show()

# order = np.log(errorMax[-1]/errorMax[-2])/np.log(stepSizes[-1]/stepSizes[-2])
# print order
#
# plt.figure()
# plt.loglog(stepSizes,errorMax,'o-')
# plt.grid(b=True, which='major', color='k', linestyle='-')
# plt.grid(b=True, which='minor', color='black', linestyle='--', alpha=0.5)
# plt.minorticks_on()
# plt.show()

# PLOTTING ENERGY

avgEn = np.average(tot)
maxEn = np.max(tot)
minEn = np.min(tot)
absError = max(maxEn-avgEn,avgEn-minEn)
relError = absError/avgEn;
print "Relative error: %.2f%%\n"%(relError*100)

# plt.figure()
# plt.plot(pot,label='potential')
# plt.plot(kin,label='kinetic')
# plt.plot(tot,label='total')
# plt.legend(loc='lower left')
# plt.show()
