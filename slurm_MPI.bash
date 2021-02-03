#!/bin/bash

#SBATCH -J MPI_TEST       # job name
#SBATCH -o log_slurm.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 16             # total number of tasks/cores requested
#SBATCH -N 1 		  # number of nodes you want to run on	
#SBATCH -p classroom      # queue (partition) -- defq, eduq, gpuq, shortq
#SBATCH -t 00:05:00       # run time (hh:mm:ss) - 12.0 hours in this example.

# Mail alert at start, end and abortion of execution
#SBATCH --mail-type=FAIL,END

# Send mail to this address
#SBATCH --mail-user=ravishankar@u.boisestate.edu

# Generally needed modules:
module load mpich/ge/gcc/64/3.2.1
module load slurm

# Execute the program
mpirun -np 16 ./bin/mm ./input/matrix4x4.m ./input/matrix4x4.m output4x4.m

# Exit if mpirun errored:
status=$?
if [ $status -ne 0 ]; then
    exit $status
fi

