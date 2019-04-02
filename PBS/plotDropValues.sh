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
unset LAMBDA_RANGE # Just in case we've somehow inherited a variable by that name externally.
LAMBDA_RANGE[2]=0..6
LAMBDA_RANGE[3]=0..4

# Set value of N from the PBS_ARRAY_INDEX
N=${PBS_ARRAY_INDEX}

# Make sure we have a valid case (that a LAMBDA_RANGE is assigned for our ${N}).
if [ -z ${LAMBDA_RANGE[${N}]+defined} ]; then 
	echo Invalid case: N=${PBS_ARRAY_INDEX} >&2
	exit 1
fi

# Output the run time information for the grid.
${PBS_O_WORKDIR}/gridRuntimeInfo.sh

# Run the work.
cd ${PBS_O_WORKDIR}/../Maple
PLOT_RANGE=${LAMBDA_RANGE[${N}]}
time maple -q -c "N:=${N};" -c "PLOT_RANGE:=${PLOT_RANGE};" -c "NUMTHREADS:=36;" plotDropValues.mpl
