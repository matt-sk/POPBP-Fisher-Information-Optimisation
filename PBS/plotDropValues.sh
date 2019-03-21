#!/bin/sh

#PBS -l select=1:ncpus=36:mem=16GB
#PBS -l walltime=300:00:00           
#PBS -k oe
#PBS -J 2-3

# Load Modules
source /etc/profile.d/modules.sh
module load maple/2017 gcc/7.2.0

# Array for lambda range. 
# We currently only have ranges for n=2 and n-3 cases.
LAMBDA_RANGE[2]=0..6
LAMBDA_RANGE[3]=0..4

# Output the run time information for the grid.
${PBS_O_WORKDIR}/gridRuntimeInfo.sh

# Run the work.
cd ${PBS_O_WORKDIR}/../Maple
N=${PBS_ARRAY_INDEX}
PLOT_RANGE=${LAMBDA_RANGE[${PBS_ARRAY_INDEX}]}
time maple -q -c "N:=${N};" -c "PLOT_RANGE:=${PLOT_RANGE};" -c "NUMTHREADS:=36;" plotDropValues.mpl
