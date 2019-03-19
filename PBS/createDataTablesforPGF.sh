#!/bin/sh

#PBS -l select=1:ncpus=36:mem=16GB
#PBS -l walltime=5:00:00           
#PBS -k oe

# Load Modules
source /etc/profile.d/modules.sh
module load maple/2017 gcc/7.2.0

# Output the run time information for the grid.
cd ${PBS_O_WORKDIR}/../Maple
./gridRuntimeInfo.sh

# Run the work.
cd Optimisation
time maple -q -c "NUMTHREADS:=36;" createDataTablesforPGF.mpl
