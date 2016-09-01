#  @file		pincClass.c
#  @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
#  @copyright	University of Oslo, Norway
#  @brief		Framework PinC Class.
#  @date		26.10.15
#
#	This is a PINC class, that is used to control the PINC program. It has several
#	member functions that is handy to use when doing custom runs of the program.
#	It can clean up the folders as we as be used to run the program in succesion
#	with different settings.
#
#	This can be useful for dividing up the run into several parts, or to test the
#	performance of parts of the program with different settings.
#
#	(Sorry about the nonsophisticated script)
#

from pincClass import *
import subprocess
import h5py
import numpy as np

pinc = PINC(iniPath = "../../local.ini")

#Setting up wanted needed ini file
pinc.routine 		= "mgRun"
pinc.mgCycles 		= 1
pinc.startTime		= -1

class Settings:
	def __init__(self, nPre = 10, nPost = 10,
					nCoarse = 10, mgLevels = 3):
		self.nPre    = nPre
		self.nPost   = nPost
		self.nCoarse = nCoarse
		self.mgLevels	  = mgLevels
		#Store results
		self.mgCycles 	  = 0
		self.time 		  = float('Inf')

	def copy(self, copy):
		self.nPre 		= copy.nPre
		self.nPost		= copy.nPost
		self.nCoarse	= copy.nCoarse
		self.mgLevels 	= copy.mgLevels

	def setPinc(self, pinc):
		pinc.preCycles 			= self.nPre
		pinc.postCycles			= self.nPost
		pinc.coarseCycles		= self.nCoarse
		pinc.mgLevels 			= self.mgLevels
		pinc.startTime 			+= 1


def formatTimeCycles(fileName, nRuns):
	data = h5py.File(fileName,'r')
	time = np.array(data['time'][nRuns,1])
	mgCycles= np.array(data['cycles'][nRuns,1])
	data.close()

	return time, mgCycles

def nCoarse(nTries, nRun):
	nRun

bestRun = Settings()
currentRun = Settings(100,100,100,2)

pinc.clean()
nTries 	= 100
nRun	= 0
preInc 	= 1

for i in range(0, 3):	#mgLevels
	# prevPreTime = float('Inf')
	prevCoaTime = float('Inf')
	# prevPosTime	= float('Inf')
	for j in range(0,nTries): #nCoarse
		##Run, retrieve time and cycles used
		currentRun.setPinc(pinc)
		pinc.runMG()
		time, mgCycles = formatTimeCycles('test_timer.xy.h5',nRun)

		#Check if best run
		if(time < bestRun.time):
			bestRun.copy(currentRun)
			bestRun.time = time
			bestRun.mgCycles = mgCycles

		if(preInc == 1):
			if(time < prevCoaTime):
				prevCoaTime = time
				currentRun.nCoarse *= 2
			else:
				currentRun.nCoarse *= 0.25
				preInc = -1
		else:
			if(time < prevCoaTime):
				prevCoaTime = time
				currentRun.nCoarse *= 0.5
			else:
				prevCoaTime = float('Inf')
				preInce = 1
				break

		nRun += 1

	currentRun.mgLevels += 1


print "\nBest runtime \t= %e"	, bestRun.time, "ns"
print "\nProposed run:"
print "mgCycles \t= "		, bestRun.mgCycles
print "mgLevels \t= "		, bestRun.mgLevels
print "nPreSmooth \t= "		, bestRun.nPre
print "nPostSmooth \t= "	, bestRun.nPost
print "nCoarseSolve \t= "	, bestRun.nCoarse
