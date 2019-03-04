#!/bin/sh

#PBS -l select=1:ncpus=1:mem=8GB
#PBS -l walltime=20:00:00
#PBS -k oe
#PBS -J 1-7

# Load Modules
source /etc/profile.d/modules.sh
module load maple/2017

# Run the work
cd ${PBS_O_WORKDIR}

# Output the run time information for the grid.
./gridRuntimeInfo.sh

#time maple -q -c "LAMBDA:=${PBS_ARRAY_INDEX}.0;" plot-ti-star-2.mpl
time maple -q -c "N:=2;" -c "L:=${PBS_ARRAY_INDEX};" -c "NUMTHREADS:=1;" plot-ti-star.mpl
