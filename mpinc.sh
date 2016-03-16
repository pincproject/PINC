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

# Path of caller (this script's folder)
#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR="$( cd "$( dirname "$0" )" && pwd )"


# Get number of processes from PINC
NP=`$DIR/pinc "$@" getnp`

# Run
if [ "$NP" -eq 1 ]
then
	$DIR/pinc "$@"
else
	mpirun -np "$NP" $DIR/pinc "$@"
fi
