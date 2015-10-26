#!/bin/bash
#
# Execution info script
# Sigvald Marholm, 06.10.15
#
# Description:
#   Run along with execution of simulations to store metadata
#   about the simulation necessary to claim reproducability.
#   Most importantly, this includes the exact version of the
#   software used to acheive the results, through Git tag.
#   Second comes versions of compiler, auxilliary software and
#   time of execution.

# Expected variables:
#   $EXECPATH - Path of somewhere in repository

# Name of info-files
INFOFILE="execinfo.txt"
OMPIFILE="ompiinfo.txt"

THISDIR=`pwd`

# Generate execinfo-file
cd $EXECPATH

TIME=`date +"%m.%d.%y %H:%M:%S"`
HASH=`git rev-parse HEAD`
TAG=`git describe --tags --abbrev=0 --always`

echo "Time: $TIME" > $INFOFILE
echo "Git:" >> $INFOFILE
echo "  Tag: $TAG" >> $INFOFILE
echo "  Hash: $HASH" >> $INFOFILE
echo "Additional flags: $ADDCFLAGS" >> $INFOFILE

mv $INFOFILE $THISDIR/

# Generate OMPI-info-file
ompi_info > $OMPIFILE
mv $OMPIFILE $THISDIR/

# Reset pwd
cd $THISDIR
