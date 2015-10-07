#!/bin/bash
#
# Plain (non-MPI) start script
# Sigvald Marholm, 06.10.15

# DiP3D folder relative to $SUBMITDIR
export EXECPATH="../../mn-fysrp-pic/DiP3D"
export SUBMITDIR=`pwd`
export SCRATCH=$SUBMITDIR/scratch	# Equivalent of $SCRATCH on Abel
mkdir $SCRATCH

# Store execution information
./execinfo.sh

# Copy source and input files to scratch area
cp -r $SUBMITDIR/$EXECPATH/src $SCRATCH/
cp $SUBMITDIR/input.txt $SCRATCH/
cp $SUBMITDIR/sphere.txt $SCRATCH/

# Create data folder
mkdir $SCRATCH/data

# Build DiP3D
cd $SCRATCH/src
make clean
make
cd ..

# Run simulation
./src/dust sphere.txt

# Copy data back and remove $SCRATCH
cp -r $SCRATCH/data/* $SUBMITDIR
rm -r $SCRATCH
