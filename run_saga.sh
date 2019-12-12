#!/bin/bash

# Job name:
#SBATCH --job-name=PINC-test
#
# Project:
#SBATCH --account=nn9299k
#
# Wall clock limit:
#SBATCH --time=1:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=4G
#
#SBATCH --ntasks=16
#
# Outfile
#SBATCH --output=pinc-%j.out

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

###export LD_LIBRARY_PATH=/usit/abel/u1/sigvaldm/hdf5-1.8.17/lib:$LD_LIBRARY_PATH
module load FFTW/3.3.8-gompi-2019a
module load GSL/2.5-GCC-8.2.0-2.31.1
module load HDF5/1.10.5-gompi-2019a

module list 

#Compile
echo "compiling"
make clean
make
## Copy input files to the work directory:
mkdir $USERWORK/$SLURM_JOB_ID
cp $SUBMITDIR/pinc $USERWORK/$SLURM_JOB_ID
cp $SUBMITDIR/FBinstability_SI.ini $USERWORK/$SLURM_JOB_ID
mkdir $USERWORK/$SLURM_JOB_ID/data

## Make sure the results are copied back to the submit directory (see Work Directory below):

chkfile $USERWORK/$SLURM_JOB_ID/data
##chkfile $USERWORK/$SLURM_JOB_ID/CollisionDump.txt

## Do some work:
cd $USERWORK/$SLURM_JOB_ID
echo "Running pinc"
mpirun -np 64 pinc FBinstability_SI.ini
