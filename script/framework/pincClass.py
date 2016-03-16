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
	def __init__(self, pincPath = "../../mpinc.sh", iniPath = "../../input.ini"):
		self.nTimeSteps = 100
		self.pincPath = pincPath
		self.iniPath = iniPath
		self.preCycles = 1
		self.postCycles = 1
		self.coarseCycles = 1
		self.mgLevels = 4
		self.mgCycles = 3

	def runCommand(self, cmd):
		subprocess.call(cmd,shell=True)

	def clean(self):
		self.runCommand("rm *.h5")

	def runMG(self):
		cmd = self.pincPath + " " + self.iniPath
		cmd += " multigrid:mgLevels=" + str(self.mgLevels)
		self.runCommand(cmd)
