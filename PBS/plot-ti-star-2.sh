#!/bin/sh

#PBS -l select=1:ncpus=1:mem=8GB
#PBS -l walltime=20:00:00
#PBS -k oe
#PBS -J 1-7

# Load Modules
source /etc/profile.d/modules.sh
module load maple/2017 gcc/7.2.0

# Output the run time information for the grid.
${PBS_O_WORKDIR}/gridRuntimeInfo.sh

# Run the work
cd ${PBS_O_WORKDIR}../Maple
time maple -q -c "N:=2;" -c "L:=${PBS_ARRAY_INDEX};" -c "NUMTHREADS:=1;" plot-ti-star.mpl
