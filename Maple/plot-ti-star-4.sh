#!/bin/sh

#PBS -l select=1:ncpus=36:mem=16GB
#PBS -l walltime=300:00:00           
#PBS -k oe
#PBS -J 1-3

# Load Modules
source /etc/profile.d/modules.sh
module load maple/2017

# Run the work
cd ${PBS_O_WORKDIR}

# Output the run time information for the grid.
./gridRuntimeInfo.sh

#time maple -q -c "NUMTHREADS:=32;" -c "LAMBDA:=0.5;" plot-ti-star-4.mpl
time maple -q -c "N:=4;" -c "L:=${PBS_ARRAY_INDEX};" -c "NUMTHREADS:=36;" plot-ti-star.mpl
