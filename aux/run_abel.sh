#!/bin/bash

# NB:
# Using --mem-per-cpu higher than 3936 is penalized,
# i.e. you pay for two (or more) cores. It would be
# more efficient to instead increase the number of
# cores and make some use of what you pay for.
# C.f. Processor equivalents in Abel documentation.

#SBATCH --job-name=PINC
#SBATCH --account=nn9299k
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=3936
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4


## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

module load openmpi.gnu/1.10.2
module load gsl/1.16
module load fftw/3.3.4

export LD_LIBRARY_PATH=/usit/abel/u1/sigvaldm/hdf5-1.8.17/lib/:$LD_LIBRARY_PATH   

# Compile
# Strictly speaking it is not necessary to do this for every job
make veryclean
make

# Run
./mpinc.sh langmuirCold.ini &> pinc.log

