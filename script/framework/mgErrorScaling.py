#  @file		mgErrorScaling.py
#  @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
#  @copyright	University of Oslo, Norway
#  @brief		Investigate error scaling
#  @date		26.10.15
#
#
#

import sys, os, inspect

# Workaround to import from different folder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../plot")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import h5py
import numpy as np
import pylab as plt
from pinc import *
from utility import *

# if(len(sys.argv) > 1):
#     ini = "../../" + sys.argv[1]
# else:
#     ini = "../../local.ini"

pinc = Pinc(ini = "input/mgErrorScaling.ini")
path = '../../data/'

pinc['methods:mode'] 		= 'mgModeErrorScaling'
pinc['grid:nSubdomains']   	= np.array([2,1,1])
trueSize 					= np.array([8,8,8])
pinc['grid:trueSize']       = trueSize
pinc['grid:nDims']          = 3
pinc['time:startTime']      = 0


dim = 0

pinc.clean()

nTest = 4

meanE2Phi      = np.zeros(nTest)
meanE2E        = np.zeros(nTest)
stepSize       = np.zeros(nTest)

for n in range(nTest):
	pinc.run()
	dataPath            = path + 'error_' + str(n) + '.grid.h5'
	errorPhi            = transformData(dim, h5py.File(dataPath,'r'), 0., average=False)
	dataPath            = path + 'errorE_' + str(n) + '.grid.h5'
	errorE              = transformData(dim, h5py.File(dataPath,'r'), 0., average=False)
	#
	meanE2E[n]           = np.sqrt(np.sum(errorE*errorE)/trueSize[dim])
	meanE2Phi[n]         = np.sqrt(np.sum(errorPhi*errorPhi)/trueSize[dim])
	stepSize[n]          = 1./np.product(trueSize)
	trueSize[dim] 		*= 2
	pinc['grid:trueSize']  = trueSize
	pinc['time:startTime']      += 1
	del errorPhi
	del errorE

print np.log(meanE2E[nTest-1]/meanE2E[nTest-2])/np.log(stepSize[nTest-1]/stepSize[nTest-2])
print np.log(meanE2Phi[nTest-1]/meanE2Phi[nTest-2])/np.log(stepSize[nTest-1]/stepSize[nTest-2])
# plotScatterLogLog('Numerical Error in Electric Field', stepSize, meanE2E)
plotScatterLogLog('Numerical Error in Potential', stepSize, meanE2Phi)

plt.show()
