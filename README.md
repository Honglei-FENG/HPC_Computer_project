# HPC_Computer_Project

In the explicit equation folder, run make "explcit.out" to build the executable file.
In the implicit equation folder, run make "implcit.out" to build the executable file.

Use "-n " to modify the number of grids and "-dt " to modify the time step.

Use "-rread " to restart facility using HDF5.

If you want to run it directly, you can directly execute "bsub<script.sh" after the executable file is generated.

The output is in "emplicit.h5" or "implicit.h5", and the parameters of dx and dt are recorded in "emplicit.dat" or "implicit.dat".

The answer of each method is at the top of the performance analysis, and the implicit scheme will make the results at the bottom of the comparison.