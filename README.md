# HPC_Computer_Project

In the explicit equation folder, run make "explcit.out" to build the executable file.
In the implicit equation folder, run make "implcit.out" to build the executable file.

Use "-n " to modify the number of grids and "-dt " to modify the time step.

Use "-restart " to restart facility using HDF5.

If you want to run it directly, you can directly execute "bsub<script.sh" after the executable file is generated.

The output is in "emplicit.h5" or "implicit.h5".

The answer of each method is at the top of the performance analysis, and the implicit scheme will make the results at the bottom of the comparison.

There are some errors in the report due to the low accuracy of the results of the code submitted. Please refer to it as appropriate. The latest code has sufficient accuracy, but it does not match the report.
