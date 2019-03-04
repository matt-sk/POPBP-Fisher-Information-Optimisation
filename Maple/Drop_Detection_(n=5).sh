#!/bin/sh

#PBS -l select=1:ncpus=36:mem=16GB
#PBS -l walltime=300:00:00
#PBS -k oe
#PBS -J 1-2

# Load Modules
source /etc/profile.d/modules.sh
module load maple/2017

# Run the work
cd ${PBS_O_WORKDIR}

# Output the run time information for the grid.
./gridRuntimeInfo.sh

# Run the drop detection.
time maple -q -c "NUMTHREADS:=36;" -c "THRESHOLD:=6;" -c "N:=5;" -c "L:=${PBS_ARRAY_INDEX};" Drop\ Detection.mpl
