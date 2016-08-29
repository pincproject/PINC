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



from pincClass import *
import subprocess
import h5py
import numpy as np

pinc = PINC()
# pinc.clean()
#Setting up wanted needed ini file
pinc.routine = "mgRoutine"
pinc.trueSize = [16,4,4]
pinc.mgCycles = 1
pinc.preCycles = 10
pinc.postCycles = 10
pinc.coarseCycles = 100
pinc.nSubdomains = [2,1,1]
pinc.preSmooth = 	"jacobianND"
pinc.postSmooth= 	"jacobianND"
# pinc.coarseSolver=	"jacobianND"
# pinc.preSmooth = 	"gaussSeidelRB"
# pinc.postSmooth= 	"gaussSeidelRB"
pinc.coarseSolver=	"gaussSeidelRBND"





pinc.clean()

for levels in range(1,5):
	for cycles in range(1,100,10):
		pinc.mgLevels = levels
		pinc.preCycles = cycles
		pinc.postCycles = cycles
		pinc coarseCycles = cycles
		pinc.runMG()
		pinc.startTime = pinc.startTime+1



data = h5py.File('test_timer.xy.h5','r')
time= data['time']
cycles= np.array([data['cycles'][:,1]])

data = np.concatenate((time, cycles.T), axis = 1)

print data
