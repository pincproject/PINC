# @file			pincClass.py
# @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
#  				Sigvald Marholm <sigvaldm@fys.uio.no>
# @copyright	University of Oslo, Norway
# @brief		Framework PinC Class.
# @date			26.10.15
#
# This is a PINC class, that is used to control the PINC program. It hasseveral
# member functions that is handy to use when doing custom runs of the program.
# It can clean up the folders as we as be used to run the program in succesion
# with different settings.
#
# This can be useful for dividing up the run into several parts, or to test the
# performance of parts of the program with different settings.
#

import subprocess
import numpy as np

class Pinc(dict):
	def __init__(self, pinc="./mpinc.sh", ini="langmuir.ini", path="../.."):
		# All commands will be executed from "path"

		self.pinc = pinc
		self.ini = ini
		self.path = path

	def run(self):
		cmd = self.pinc + " " + self.ini
		for key in self:
			cmd += " " + key + "=" + self.parse(key)
		self.runCommand(cmd)

	def runCommand(self, cmd):
		cmd = "cd " + self.path + "; " + cmd
		subprocess.call(cmd,shell=True)

	def clean(self):
		self.runCommand("rm -f data/*.h5")
		self.runCommand("rm -f data/*.txt")

	def parse(self, key):
		value = self[key]
		if isinstance(value,(list,np.ndarray)):
			string = str(value[0])
			if(len(value) > 1):
				for l in range(1,len(value)):
					string += "," + str(value[l])
			return string
		elif isinstance(value,(int,float)):
			return str(value)
		else:
			return value
