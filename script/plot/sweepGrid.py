import h5py
import pylab as plt
import numpy as np
import subprocess as sp

length = 2.0*np.pi
trueSizeList = [16,32,64,128]
stepSizeList = []
errorList = []

prefix = "../../"
pincCmd = prefix + "mpinc.sh " + prefix + "input.ini "

for trueSize in trueSizeList:

	stepSize = length/trueSize
	stepSizeList.append(stepSize)
	pincArg = "grid:stepSize=%f,%f,%f grid:trueSize=%d,32,32 "%(stepSize,stepSize,stepSize,trueSize)
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

for n in range(len(trueSizeList)):
	print "trueSize: %5f, error: %5.3f%%"%(trueSizeList[n],errorList[n])
