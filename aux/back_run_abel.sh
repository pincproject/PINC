#!/bin/bash

# Job name:
#SBATCH --job-name=PinC-FB1
#
# Project:
#SBATCH --account=nn9299k
#
# Wall clock limit:
#SBATCH --time=110:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=3936M
#
#SBATCH --ntasks=64
#
#SBATCH --ntasks-per-node=16
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


##echo $PATH
#Compile
echo "compiling"
make clean
make
## Copy input files to the work directory:
mkdir $SCRATCH/steffemb
cp $SUBMITDIR/pinc $SCRATCH/steffemb
cp $SUBMITDIR/FBinstability_SI.ini $SCRATCH/steffemb
mkdir $SCRATCH/steffemb/data

## Make sure the results are copied back to the submit directory (see Work Directory below):

chkfile $SCRATCH/steffemb/data
chkfile $SCRATCH/steffemb/CollisionDump.txt

## Do some work:
cd $SCRATCH/steffemb
##echo "compiling"
##make clean
##make
echo "Running pinc"
mpirun -np 64 pinc FBinstability_SI.ini
