#!/bin/bash

# Submit the PBS job script to the queue, catching the job ID as we do so.
if [ $# -eq 1 ] && [ $1 = 'express' ]
then
	JOB_ID=$(qsub -q express -P exp-XXXXX run_cilia_sim.pbs)
	echo "Using the express queue..."
else
	JOB_ID=$(qsub run_cilia_sim.pbs)
	echo "Using the normal queue..."
fi
JOB_ID=${JOB_ID:0:-4} # Throw away the .pbs extension

# Make the job-specific directory to store everything related to this job.
LOC=$JOB_ID
mkdir $LOC

# Copy the source files, along with any input files, to this directory.
# The executable will be compiled as part of the job, ensuring we have a preserved copy of the .par etc. used for a given simulation.
# Note: large input files should probably be copied as part of the job in case they take longer to copy than the job spends in the queue. 
cp -r src $LOC/
cp makefile $LOC/
cp config.hpp $LOC/
cp *.input $LOC/

# And we should be done!
echo "Job submitted!"
