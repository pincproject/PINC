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
pinc.trueSize = [32,8,8]
pinc.preCycles = 100
pinc.postCycles = 100
pinc.coarseCycles = 100
pinc.nSubdomains = [4,1,1]
pinc.clean()

for i in range(2,6):
	pinc.mgLevels = i
	pinc.runMG()
	pinc.startTime = pinc.startTime+1


data = h5py.File('test_timer.xy.h5','r')
time= data['time']
cycles= data['cycles']

data = np.concatenate((time, cycles[:,0].T), axis = 1)

print data
