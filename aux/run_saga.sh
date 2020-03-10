#!/bin/bash

# Job name:
#SBATCH --job-name=daedalus1-PINC
#
# Project:
#SBATCH --account=nn9299k
#
# Wall clock limit:
#SBATCH --time=7-00:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=4G
#
#SBATCH --ntasks=8
#
# Outfile
#SBATCH --output=pinc-%j.out

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

###export LD_LIBRARY_PATH=/usit/abel/u1/sigvaldm/hdf5-1.8.17/lib:$LD_LIBRARY_PATH
module load FFTW/3.3.8-gompi-2019a
##module load GSL/2.5-GCC-7.3.0-2.30
##module load FFTW/3.3.8-gompic-2018b
module load GSL/2.5-GCC-8.2.0-2.31.1
module load HDF5/1.10.5-gompi-2019a
##module swap OpenMPI/3.1.3-GCC-8.2.0-2.31.1 OpenMPI/3.1.4-GCC-8.3.0
module list

#Compile
echo "compiling"
make clean
make
## Copy input files to the work directory:
mkdir $SCRATCH//$SLURM_JOB_ID
cp $SUBMITDIR/pinc $SCRATCH/$SLURM_JOB_ID
cp $SUBMITDIR/daedalus1.ini $SCRATCH/$SLURM_JOB_ID
mkdir $SCRATCH/$SLURM_JOB_ID/data
cp $SUBMITDIR/object.grid.h5 $SCRATCH/$SLURM_JOB_ID/data

## Make sure the results are copied back to the submit directory (see Work Directory below):

savefile $SCRATCH/$SLURM_JOB_ID/data
##chkfile $USERWORK/$SLURM_JOB_ID/CollisionDump.txt

## Do some work:
cd $SCRATCH/$SLURM_JOB_ID
echo "local folder contents: "
ls
echo "data folder contents:"
ls data/
echo "Running pinc"
srun ./pinc daedalus1.ini


