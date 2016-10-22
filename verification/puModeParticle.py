import h5py
import pylab as plt
import numpy as np
from itertools import izip as zip, count
from math import pi

import sys
sys.path.append('../framework')
from pinc import *

pinc = Pinc(ini = "local.ini")

test = 'dt'
assert test in ['dt','dx']

# ANALYTICAL TRAJECTORY
if test=='dt':
	start = 2
	stop = 16
if test=='dx':
	start = 2
	stop = 20


timeSteps = np.array([2**(-x) for x in range(start,stop)])
trueSizes = np.array([2**x    for x in range(start,stop)])

if(test=='dt'): trueSizes[:] = 8
stepSizes = 2.0/trueSizes
if(test=='dx'): timeSteps[:] = 0.00001#0.5*stepSizes

pinc["time:nTimeSteps"] = 100
omega = 1

errorMax = np.zeros((2,len(timeSteps)))
error2 = np.zeros((2,len(timeSteps)))

accs = ['puAccND0KE','puAccND1KE']

for m,acc in enumerate(accs):
	pinc["methods:acc"] = acc

	for (i,timeStep,trueSize,stepSize) in zip(count(),timeSteps,trueSizes,stepSizes):

		print("Run %d of %d (dt,dx)=(%f,%f)"%(i+1,stop-start,timeStep,stepSize))

		pinc["time:timeStep"] = timeStep
		pinc["grid:trueSize"] = np.array([trueSize])

		n = np.linspace(0,pinc["time:timeStep"]*pinc["time:nTimeSteps"],pinc["time:nTimeSteps"]+1)
		# posAna = pinc.trueSize[0]*(0.5-0.25*np.cos(omega*n))
		# posAna = pinc.trueSize[0]*(0.5+0.25*np.sin(omega*n))
		# posAna = (0.5-0.25*np.cos(omega*n))
		# posAna = (0.5+0.25*np.sin(omega*n))
		equilibrium = 0
		posAna = (equilibrium+0.25*np.sin(omega*n+pi/4))%1

		pinc.clean()
		pinc.run()

		# EXTRACTING ENERGY
		hist = h5py.File('../../data/history.xy.h5','r')
		pot = hist['/energy/potential/total']
		kin = hist['/energy/kinetic/total']

		kin = kin[:,1];
		pot = pot[:,1];
		tot = pot+kin;

		# EXTRACTING TRAJECTORY OF SINGLE PARTICLE
		pop = h5py.File('../../data/pop.pop.h5','r')
		N = len(pop['/pos/specie 0/'])
		pos = np.zeros(N)
		for j in xrange(N):
			pos[j]=pop['/pos/specie 0/n=%.1f'%j][:]

		pos /= pinc["grid:trueSize"][0]
		pos%=1

		# PLOTTING TRAJECTORIES

		if i==0:
			plt.plot(pos,label='Numerical')
			plt.plot(posAna,label='Analytical')
			plt.legend()
			plt.show()
			plt.legend(loc='lower left')

		diff = pos-posAna
		errorMax[m][i] = np.max(abs(diff))
		error2[m][i] = np.sqrt(np.sum(diff**2)/(pinc["time:nTimeSteps"]+1))

xaxis = timeSteps if test=='dt' else stepSizes
error = errorMax

order = np.log(error[-1]/error[-2])/np.log(xaxis[-1]/xaxis[-2])
print order

plt.figure()
plt.loglog(xaxis,error[0],'o-b',label="puAccND0")
plt.loglog(xaxis,error[1],'o-g',label="puAccND1")
if test=='dt':
	plt.loglog(xaxis,(error[0][-1]/xaxis[-1]**2)*xaxis**2,'--b',label="O(dt^2)")
	plt.loglog(xaxis,(error[1][-1]/xaxis[-1]**3)*xaxis**3,'--g',label="O(dt^3)")
if test=='dx':
	plt.loglog(xaxis,(error[0][7]/xaxis[7])*xaxis,'--',label="O(dx)")
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='black', linestyle='--', alpha=0.5)
plt.xlabel(test)
plt.ylabel("Max Error")
plt.title("Error in oscillating particle after 100 steps")
plt.minorticks_on()
plt.legend(loc='upper left')
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

# avgEn = np.average(tot)
# maxEn = np.max(tot)
# minEn = np.min(tot)
# absError = max(maxEn-avgEn,avgEn-minEn)
# relError = absError/avgEn;
# print "Relative error: %.2f%%\n"%(relError*100)

# plt.figure()
# plt.plot(pot,label='potential')
# plt.plot(kin,label='kinetic')
# plt.plot(tot,label='total')
# plt.legend(loc='lower left')
# plt.show()
