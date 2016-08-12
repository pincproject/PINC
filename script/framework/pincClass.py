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

import subprocess


class PINC:
	def __init__(self, pincPath = "../../mpinc.sh", iniPath = "../../local.ini"):
		self.nTimeSteps = 100
		self.pincPath = pincPath
		self.iniPath = iniPath
		self.trueSize = [32,16,16]
		self.preCycles = 1
		self.postCycles = 1
		self.coarseCycles = 1
		self.mgLevels = 4
		self.mgCycles = 1
		self.routine = "regular"
		self.startTime = 0

	def runCommand(self, cmd):
		subprocess.call(cmd,shell=True)

	def clean(self):
		self.runCommand("rm *.h5")

	def runMG(self):
		cmd = self.pincPath + " " + self.iniPath
		cmd += " main:routine=" + self.routine
		cmd += " time:startTime=" + str(self.startTime)
		cmd += " grid:trueSize=" + str(self.trueSize[0]) + "," + str(self.trueSize[1]) + "," + str(self.trueSize[2])
		cmd += " multigrid:mgLevels=" + str(self.mgLevels)
		cmd += " multigrid:mgCycles=" + str(self.mgCycles)
		cmd += " multigrid:nPreSmooth=" + str(self.preCycles)
		cmd += " multigrid:nPostSmooth=" + str(self.postCycles)
		cmd += " multigrid:nCoarseSolve=" + str(self.coarseCycles)

		self.runCommand(cmd)
