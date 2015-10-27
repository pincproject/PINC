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

NUM=1
MAX=132
while read -r line; do
	if [ ${#line} -gt $MAX ] ; then
		echo $1:$NUM:${#line}: warning: line exceeds $MAX columns
	fi
	NUM=$((NUM+1))
done < $1
