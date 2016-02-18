#!/bin/sh

#
# @file			mpinc.sh
# @author		Sigvald Marholm <sigvaldm@fys.uio.no>
# @copyright	University of Oslo, Norway
# @brief		PINC-wrapper for MPI
# @date			02.02.16
#
# Automatically determines number of MPI processes to use when
# calling PINC, based on feedback from PINC itself. Normally,
# The number of processes will be determined by the input file,
# but these settings can be overridden through arguments to PINC
# which is why it is better to fetch them from PINC itself.
#

NP=`./pinc "$@" getnp`
if [ "$NP" -eq 1 ]
then
	./pinc "$@"
else
	mpirun -np "$NP" ./pinc "$@"
fi
