#!/bin/sh

#
# @file			check.sh
# @author		Sigvald Marholm <sigvaldm@fys.uio.no>
# @copyright	University of Oslo, Norway
# @brief		Custom sanity check on source files
# @date			12.10.15
#
# Custom sanity checks on source file.
# Currently only used to check line width.
#

NUM=1			# Number of current line
MAX=132			# Max line length
PREFMAX=80		# Max preferable line length
PREFFRACT=1		# Percentage of lines allowed to exceed max preferable
UNPREFLINES=0	# Number of lines exceeding preferrable

while read -r line; do
	if [ ${#line} -gt $MAX ] ; then
		echo $1:$NUM:${#line}: warning: line exceeds $MAX columns
	fi
	if [ ${#line} -gt $PREFMAX ] ; then
		UNPREFLINES=$((UNPREFLINES+1))
	fi
	NUM=$((NUM+1))
done < $1

if [ $((100*UNPREFLINES/NUM)) -gt $PREFFRACT ] ; then
	echo $1: warning: More than $PREFFRACT% of the lines exceeds $PREFMAX columns
fi
