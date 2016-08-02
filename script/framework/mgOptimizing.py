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

pinc = PINC()

# pinc.clean()
#
# data = open("mgOptData.txt", "wr")
#
# data.write("Hello moron brothers")
#
# data.close()

for i in range(3):
	pinc.runMG()
	# pinc.clean()
