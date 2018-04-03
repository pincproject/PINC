#!/bin/bash

# Job name:
#SBATCH --job-name=PinC-FB1
#
# Project:
#SBATCH --account=nn9299k
#
# Wall clock limit:
#SBATCH --time=120:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=3936
#
#SBATCH --nodes=4 --ntasks-per-node=16
#
# Outfile
#SBATCH --output=pinc-%j.out

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

export LD_LIBRARY_PATH=/usit/abel/u1/sigvaldm/hdf5-1.8.17/lib:$LD_LIBRARY_PATH


module load openmpi.gnu
module load gsl/1.16
module load fftw/3.3.4

#Compile
#cd $SUBMITDIRECTORY/mn-fysrp-pic
make
#echo $PATH
#cd ..
## Copy input files to the work directory:
#cp $SUBMITDIRECTORY/mn-fysrp-pic/input/FBinstability_SI.ini $SCRATCH
#cp -R $SUBMITDIRECTORY/mn-fysrp-pic/data $SCRATCH/data
#cp $SUBMITDIRECTORY/mn-fysrp-pic/pinc $SCRATCH


## Make sure the results are copied back to the submit directory (see Work Directory below):
#chkfile /data/"*.h5" $SUBMITDIRECTORY/mn-fysrp-pic/data
#chkfile CollisionDump.txt $SUBMITDIRECTORY/mn-fysrp-pic

## Do some work:
#cd $SCRATCH
#mpirun -np 64 pinc FBinstability_SI.ini
