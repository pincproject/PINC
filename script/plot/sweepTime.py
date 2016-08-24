import h5py
import pylab as plt
import numpy as np
import subprocess as sp

totalTime = 30.0
nTimeStepsList = [15,30,60]
timeStepList = []
errorList = []

prefix = "../../"
pincCmd = prefix + "mpinc.sh " + prefix + "input.ini "

for nTimeSteps in nTimeStepsList:

	timeStep = totalTime/nTimeSteps
	timeStepList.append(timeStep)
	pincArg = "time:nTimeSteps=%d time:timeStep=%f "%(nTimeSteps,timeStep)
	print pincCmd+pincArg
	sp.call("rm *.h5",shell=True)
	sp.call(pincCmd+pincArg,shell=True)


	hist = h5py.File('test_history.xy.h5','r')
	pot = hist['/energy/potential/total']
	kin = hist['/energy/kinetic/total']

	kin = kin[:,1];		# Extract y-axis
	pot = -pot[:,1];	# Extract y-axis and invert
	tot = pot+kin;		# Collect total energy

	avgEn = np.average(tot)
	maxEn = np.max(tot)
	minEn = np.min(tot)
	absError = max(maxEn-avgEn,avgEn-minEn)
	relError = absError/avgEn;
	percentError = relError*100;
	errorList.append(percentError)

for n in range(len(nTimeStepsList)):
	print "timeStep: %5f, error: %5.3f%%"%(timeStepList[n],errorList[n])
