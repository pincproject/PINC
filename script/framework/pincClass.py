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
		self.pincPath = pincPath
		self.iniPath = iniPath
		self.startTime = 0

	def runCommand(self, cmd):
		subprocess.call(cmd,shell=True)

	def clean(self):
		self.runCommand("rm *.h5")
		self.runCommand("rm *.txt")


	def arrToStr(self, array):
		string = str(array[0])
		if(len(array) > 1):
			for l in range(1,len(array)):
				string += "," + str(array[l])

		return string

	def runMG(self):
		cmd = self.pincPath + " " + self.iniPath
		cmd += " methods:mode=" + self.mode
		cmd += " time:startTime=" + str(self.startTime)
		cmd += " multigrid:mgLevels=" + str(self.mgLevels)
		cmd += " multigrid:mgCycles=" + str(self.mgCycles)
		cmd += " multigrid:nPreSmooth=" + str(self.preCycles)
		cmd += " multigrid:nPostSmooth=" + str(self.postCycles)
		cmd += " multigrid:nCoarseSolve=" + str(self.coarseCycles)

		self.runCommand(cmd)

	def mgErrorScaling(self):
		cmd = self.pincPath + " " + self.iniPath
		cmd += " methods:mode=" + self.mode
		cmd += " grid:nDims="	+ str(self.nDims)
		cmd += " grid:trueSize=" + self.arrToStr(self.trueSize)
		cmd += " grid:nSubdomains="	+self.arrToStr(self.nSubdomains)
		cmd += " time:startTime=" + str(self.startTime)

		self.runCommand(cmd)
