#!/bin/sh

#PBS -l select=1:ncpus=36:mem=16GB
#PBS -l walltime=300:00:00           
#PBS -k oe
#PBS -J 1-6

# Load Modules
source /etc/profile.d/modules.sh
module load maple/2017

# Run the work
cd ${PBS_O_WORKDIR}

# Output the run time information for the grid.
./gridRuntimeInfo.sh

#time maple -q -c "NUMTHREADS:=36;" -c "LAMBDA:=${PBS_ARRAY_INDEX}.0;" plot-ti-star-3.mpl
time maple -q -c "N:=3;" -c "L:=${PBS_ARRAY_INDEX};" -c "NUMTHREADS:=36;" plot-ti-star.mpl
