import h5py
import pylab as plt
import numpy as np
import sys
sys.path.append('../framework')
from pinc import *

assert sys.argv[1] in ['dx','dt']
mode = sys.argv[1]

if len(sys.argv) >= 3:
	ini = sys.argv[2]
else:
	ini = "langmuir.ini"

if len(sys.argv) >= 5:
	start = sys.argv[3]
	stop  = sys.argv[4]
else:
	start = 0
	stop = 0

pinc = Pinc(ini=ini)

if mode=='dt':
	if start==0: start = 3
	if stop ==0: stop  = 11

	nTimeStepss = np.array([2**x for x in range(start,stop)])
	trueSizes   = 128*np.ones(nTimeStepss.shape)

if mode=='dx':
	if start==0: start = 2
	if stop ==0: stop  = 20

	trueSizes   = np.array([2**x for x in range(start,stop)])
	nTimeStepss = 100*np.ones(trueSizes.shape)

stepSizes = 6.28/trueSizes
timeSteps = 15.0/nTimeStepss

pinc["grid:nDims"]=1
pinc["methods:acc"]='puAccND1KE'
pinc["methods:distr"]='puDistrND1'
pinc["methods:migrate"]='puExtractEmigrantsND'

error = np.zeros(timeSteps.shape)

for i in range(len(timeSteps)):

	timeStep = timeSteps[i]
	stepSize = stepSizes[i]
	trueSize = trueSizes[i]
	nTimeSteps = nTimeStepss[i]

	pinc["time:timeStep"] = timeStep
	pinc["time:nTimeSteps"] = nTimeSteps
	pinc["grid:trueSize"] = trueSize
	pinc["grid:stepSize"] = stepSize

	pinc.clean()
	pinc.run()

	hist = h5py.File('../../data/history.xy.h5','r')
	pot = hist['/energy/potential/total']
	kin = hist['/energy/kinetic/total']

	kin = kin[:,1];
	pot = pot[:,1];
	tot = pot+kin;

	error[i] = max(abs(tot-tot[0]))/tot[0]


xaxis = timeSteps if mode=='dt' else stepSizes

order = np.log(error[-1]/error[-2])/np.log(xaxis[-1]/xaxis[-2])
print order

plt.figure()
plt.loglog(xaxis,error,'o-b',label="Measured")
if mode=='dt':
	plt.loglog(xaxis,(error[-1]/xaxis[-1])*xaxis,'--b',label="O(dt)")
if mode=='dx':
	plt.loglog(xaxis,(error[-1]/xaxis[-1]**2)*xaxis**2,'--b',label="O(dx^2)")
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='black', linestyle='--', alpha=0.5)
plt.xlabel(mode)
plt.ylabel("Relative Error in Energy")
plt.title("Error in Parameter Sweep of %s"%ini)
plt.minorticks_on()
plt.legend(loc='upper left')
plt.show()
