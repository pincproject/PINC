import h5py
import pylab as plt
import numpy as np
from itertools import izip as zip, count

import sys
sys.path.append('../framework')
from pincClass import *

ini = "../../local.ini"
pinc = PINC(iniPath = ini)
pinc.acc = "puAccND1KE"

trueSizes = np.array([2**x    for x in range(1,26)])
stepSizes = 2.0/trueSizes

omega = 1

error0 = np.zeros(trueSizes.shape)	# 0th order
error1 = np.zeros(trueSizes.shape)	# 1st order

for (i,trueSize,stepSize) in zip(count(),trueSizes,stepSizes):

	pinc.trueSize = np.array([trueSize])
	pinc.timeStep = 0.2
	pinc.nTimeSteps = 1

	exact = (0.112358)**2

	pinc.runCommand("rm *.h5")
	pinc.puErrorScaling()

	pop = h5py.File('test_pop.pop.h5','r')
	numerical = pop['/vel/specie 0/n=0.0']

	error0[i] = np.max(abs(numerical[0]-exact))
	error1[i] = np.max(abs(numerical[1]-exact))

	print("exact: %f, 0th: %f, 1st: %f"%(exact,numerical[0],numerical[1]))

lineSet = [-1,3]
order0 = np.log(error0[lineSet[1]]/error0[lineSet[0]])/np.log(stepSizes[lineSet[1]]/stepSizes[lineSet[0]])
order1 = np.log(error1[lineSet[1]]/error1[lineSet[0]])/np.log(stepSizes[lineSet[1]]/stepSizes[lineSet[0]])
print("0th order scheme error is of order %f"%order0)
print("1th order scheme error is of order %f"%order1)

plt.figure()
plt.loglog(stepSizes,error0,'o-b',label="puInterpND0")
plt.loglog(stepSizes,error1,'o-g',label="puInterpND1")
# plt.loglog(stepSizes[lineSet],error0[lineSet],'--b')
# plt.loglog(stepSizes[lineSet],error1[lineSet],'--g')
plt.loglog(stepSizes[lineSet],1e-1*stepSizes[lineSet]**1,'-b',label="O(dx^1)")
plt.loglog(stepSizes[lineSet],1e-1*stepSizes[lineSet]**2,'-g',label="O(dx^2)")
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='black', linestyle='--', alpha=0.5)
plt.xlabel('Step size')
plt.ylabel('Error')
plt.title('Interpolation of E(x)=x^2 to single point')
plt.minorticks_on()
plt.legend(loc='lower right')
plt.show()
