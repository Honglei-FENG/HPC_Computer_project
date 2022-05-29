#!/bin/bash
#BSUB -J valgrind-problem1
#BSUB -q ser
#BSUB -n 1


module purge
module load intel/2018.4
module load mpi/intel/2018.4
module load valgrind/3.14.0

valgrind mpirun ./implicit.out > $LSB_JOBID-valgrind.log 2>&1
