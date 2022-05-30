#!/bin/bash
#BSUB -J petsc-problem1
#BSUB -q ser
#BSUB -n 1


module purge
module load intel/2018.4
module load mpi/intel/2018.4

/usr/bin/time -f "run time is: %e" mpirun -np 1 ./implicit.out > $LSB_JOBID.log 2>&1
