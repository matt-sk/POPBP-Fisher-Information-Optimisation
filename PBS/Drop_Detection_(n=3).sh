#!/bin/sh

#PBS -l select=1:ncpus=36:mem=16GB
#PBS -l walltime=15:00:00
#PBS -k oe
#PBS -J 1-6

# Load Modules
source /etc/profile.d/modules.sh
module load maple/2017 gcc/7.2.0

# Output the run time information for the grid.
${PBS_O_WORKDIR}/gridRuntimeInfo.sh

# Run the drop detection.
cd ${PBS_O_WORKDIR}/../Maple
time maple -q -c "NUMTHREADS:=36;" -c "THRESHOLD:=6;" -c "N:=3;" -c "L:=${PBS_ARRAY_INDEX};" Drop\ Detection.mpl
