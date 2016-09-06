#  @file		mgErrorScaling.py
#  @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
#  @copyright	University of Oslo, Norway
#  @brief		Investigate error scaling
#  @date		26.10.15
#
#   Compares the scaling of the error of the solver, as a function of stepsize,
#   compared to an analytical solution. The test case is a sinusoidal function
#

import sys, os, inspect

# Workaround to import from different folder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../plot")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import h5py
import numpy as np
import pylab as plt
from pincClass import *
from utility import *

if(len(sys.argv) > 1):
    ini = "../../" + sys.argv[1]
else:
    ini = "../../local.ini"

pinc = PINC(iniPath = ini)

pinc.mode           = 'mgErrorScaling'
pinc.nSubdomains    = np.array([1])
pinc.trueSize       = np.array([8])
pinc.nDims          = 1
pinc.startTime      = 0

pinc.clean()

nTest = 3

maxError    = np.zeros(nTest)
meanE2      = np.zeros(nTest)
stepSize    = np.zeros(nTest)

for n in range(nTest):
    pinc.mgErrorScaling()
    dataPath = 'test_error_' + str(n) + '.grid.h5'
    error = transformData(0, h5py.File(dataPath,'r'), 0., average=False)
    maxError[n]         = np.max(np.abs(error))
    meanE2[n]           = np.sum(error*error)/pinc.trueSize[0]
    stepSize[n]         = 1./pinc.trueSize[0]
    pinc.trueSize[0]    *= 2
    pinc.startTime      += 1
    del error


plotScatterLogLog('$E_{max}$', stepSize, maxError)
plotScatterLogLog('$\\bar{E}^2$', stepSize, meanE2)

plt.show()
