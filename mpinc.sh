#!/bin/bash

#
# @file			mpinc.sh
# @brief		PINC-wrapper for MPI
# @author		Sigvald Marholm <sigvaldm@fys.uio.no>
#
# Automatically determines number of MPI processes to use when
# calling PINC, based on feedback from PINC itself. Normally,
# The number of processes will be determined by the input file,
# but these settings can be overridden through arguments to PINC
# which is why it is better to fetch them from PINC itself.
#

# Path of caller (this script's folder)
DIR="$( cd "$( dirname "$0" )" && pwd )"


# Get number of processes from PINC
NP=`$DIR/pinc "$@" getnp`

# Run
if [ "$NP" -eq 1 ]
then
	$DIR/pinc "$@"
else
	#mpirun --oversubscribe -np "$NP" $DIR/pinc "$@"
	mpirun -np "$NP" $DIR/pinc "$@"
fi
