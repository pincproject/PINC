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


preCycles 		= 2
postCycles 		= 2
coarseCycles 	= 2
mgLevels 		= np.arange(2,6)

size = mgLevels.shape[0]*preCycles.shape[0]*postCycles.shape[0]*coarseCycles.shape[0]

settings = np.empty([size, 4])

pinc.clean()
run = 0
for levels in mgLevels:
	for pre in preCycles:
		for post in postCycles:
			for coarse in coarseCycles:
				pinc.mgLevels = levels
				pinc.preCycles = pre
				pinc.postCycles = post
				pinc.coarseCycles = coarse
				pinc.runMG()

				pinc.startTime = pinc.startTime+1
				settings[run,:] = [levels, pre, post, coarse]
				run += 1

data = h5py.File('test_timer.xy.h5','r')
time= data['time']

mgCycles= np.array([data['cycles'][:,1]])

data = np.concatenate((time, mgCycles.T), axis = 1)

print data.shape
print settings.shape

result = np.concatenate((data, settings), axis = 1)

bestRun = np.argmin(result[:,1])

print "\nBest runtime \t= ", result[bestRun,1], "ns"
print "trueSize \t= ", pinc.trueSize
print "nSubdomains \t= ", pinc.nSubdomains
print "\nProposed run:"
print "mgCycles \t= ", result[bestRun, 2]
print "mgLevels \t= ", result[bestRun, 3]
print "nPreSmooth \t= ", result[bestRun, 4]
print "nPostSmooth \t= ", result[bestRun, 5]
print "nCoarseSolve \t= ", result[bestRun, 6]
