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

pinc = PINC()
# pinc.clean()
#Setting up wanted needed ini file
pinc.routine 		= "mgRun"
pinc.trueSize 		= [64,32,32]
pinc.mgCycles 		= 1
pinc.preCycles 		= 10
pinc.postCycles 	= 10
pinc.coarseCycles 	= 10
pinc.nSubdomains 	= [2,2,1]
pinc.preSmooth 		= 	"gaussSeidelRB"
pinc.postSmooth		= 	"gaussSeidelRB"
pinc.coarseSolver	=	"gaussSeidelRB"
pinc.startTime		= -1

class Settings:
	def __init__(self, preCycles = 10, postCycles = 10,
	coarseCycles = 10, mgLevels = 3):
		self.preCycles    = preCycles
		self.postCycles   = postCycles
		self.coarseCycles = coarseCycles
		self.mgLevels	  = mgLevels
		#Store results
		self.mgCycles 	  = 0
		self.time 		  = float('Inf')

	def copy(self, copy):
		self.preCycles 		= copy.preCycles
		self.postCycles		= copy.postCycles
		self.coarseCycles	= copy.coarseCycles
		self.mgLevels 		= copy.mgLevels

	def setPinc(self, pinc):
		pinc.preCycles 		= self.preCycles
		pinc.postCycles		= self.postCycles
		pinc.coarseCycles	= self.coarseCycles
		pinc.mgLevels 		= self.mgLevels
		pinc.startTime 		+= 1


def formatTimeCycles(fileName, nRuns):
	data = h5py.File(fileName,'r')
	time = np.array(data['time'][nRuns,1])
	mgCycles= np.array(data['cycles'][nRuns,1])
	data.close()

	return time, mgCycles

bestRun = Settings()
currentRun = Settings(2,2,2,3)


pinc.clean()
nTries 	= 100
nRun	= 0

for i in range(0, 1):	#mgLevels
	prevPreTime = float('Inf')
	prevCoaTime = float('Inf')
	prevPosTime	= float('Inf')
	currentRun.preCycles 	= 2
	currentRun.coarseCycles	= 2
	currentRun.postCycles	= 2
	for j in range(0,nTries): #preCycles
		##Run, retrieve time and cycles used
		for k in range(0,nTries):	#coarseCycles
			currentRun.setPinc(pinc)
			pinc.runMG()
			time, mgCycles = formatTimeCycles('test_timer.xy.h5',nRun)

			#Check if best run
			if(time < bestRun.time):
				bestRun.copy(currentRun)
				bestRun.time = time
				bestRun.mgCycles = mgCycles

			if(time < prevPreTime):
				prevPreTime = time
				currentRun.preCycles *= 2
			else:
				break

			nRun += 1

	currentRun.mgLevels += 1


print "\nBest runtime \t= "	, bestRun.time, "ns"
print "trueSize \t= "		, pinc.trueSize
print "nSubdomains \t= "	, pinc.nSubdomains
print "\nProposed run:"
print "mgCycles \t= "		, bestRun.mgCycles
print "mgLevels \t= "		, bestRun.mgLevels
print "nPreSmooth \t= "		, bestRun.preCycles
print "nPostSmooth \t= "	, bestRun.postCycles
print "nCoarseSolve \t= "	, bestRun.coarseCycles



#
# data = h5py.File('test_timer.xy.h5','r')
# time= data['time']
# cycles = data['cycles']
#
# print time
# print cycles


# for levels in mgLevels:
# 	for pre in preCycles:
# 		for post in postCycles:
# 			for coarse in coarseCycles:
# 				pinc.mgLevels = levels
# 				pinc.preCycles = pre
# 				pinc.postCycles = post
# 				pinc.coarseCycles = coarse
# 				pinc.runMG()
#
# 				pinc.startTime = pinc.startTime+1
# 				settings[run,:] = [levels, pre, post, coarse]
# 				run += 1
#
# data = h5py.File('test_timer.xy.h5','r')
# time= data['time']
#
# mgCycles= np.array([data['cycles'][:,1]])
#
# data = np.concatenate((time, mgCycles.T), axis = 1)
#
# print data.shape
# print settings.shape
#
# result = np.concatenate((data, settings), axis = 1)
#
# bestRun = np.argmin(result[:,1])
#
# print "\nBest runtime \t= ", result[bestRun,1], "ns"
# print "trueSize \t= ", pinc.trueSize
# print "nSubdomains \t= ", pinc.nSubdomains
# print "\nProposed run:"
# print "mgCycles \t= ", result[bestRun, 2]
# print "mgLevels \t= ", result[bestRun, 3]
# print "nPreSmooth \t= ", result[bestRun, 4]
# print "nPostSmooth \t= ", result[bestRun, 5]
# print "nCoarseSolve \t= ", result[bestRun, 6]
