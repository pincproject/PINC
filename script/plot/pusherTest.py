import h5py
import pylab as plt
import numpy as np

import sys
sys.path.append('../framework')
from pincClass import *

ini = "../../local.ini"
pinc = PINC(iniPath = ini)

# ANALYTICAL TRAJECTORY

timeSteps = 0.314*np.array([1,0.5,0.25,0.125,0.0625,0.03125,0.015625])
# timeSteps = 0.314*np.array([0.125,0.0625,0.03125,0.015625])
pinc.nTimeSteps = 100
pinc.trueSize = np.array([10])
omega = 1

errorMax = np.zeros(timeSteps.shape)
error2 = np.zeros(timeSteps.shape)

for (it,timeStep) in enumerate(timeSteps):

	pinc.timeStep = timeStep

	n = np.linspace(0,pinc.timeStep*pinc.nTimeSteps,pinc.nTimeSteps+1)
	posAna = pinc.trueSize[0]*(0.5-0.25*np.cos(omega*n))

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
	for i in xrange(N):
		pos[i]=pop['/pos/specie 0/n=%.1f'%i][:]

	# PLOTTING TRAJECTORIES

	# plt.plot(pos,label='Numerical')
	# plt.plot(posAna,label='Analytical')
	# plt.legend()
	# plt.show()
	# plt.legend(loc='lower left')

	diff = pos-posAna
	errorMax[it] = np.max(abs(diff))
	error2[it] = np.sqrt(np.sum(diff**2)/(pinc.nTimeSteps+1))

	pinc.timeStep /= 2

print errorMax

plt.figure()
plt.loglog(timeSteps,errorMax,'o-')
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='black', linestyle='--', alpha=0.5)
plt.minorticks_on()
plt.show()

order = np.log(errorMax[-1]/errorMax[-2])/np.log(timeSteps[-1]/timeSteps[-2])
print order

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
