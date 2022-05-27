#!/bin/bash
#BSUB -J valgrind-problem1
#BSUB -q ser
#BSUB -n 1


module purge
module load intel/2018.4
module load mpi/intel/2018.4
module load valgrind/3.14.0

valgrind --tool=memcheck --leak-check=full \
  --show-leak-kinds=all --undef-value-errors=no \
  ./explicit.out > $LSB_JOBID.log 2>&1
