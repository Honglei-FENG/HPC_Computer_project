# HPC_Computer_Project

中文版Readme

这个程序是用C语言写的，调用PETSc函数库，满足防御性编程的并行程序。

是南方科技大学的高性能计算的期末编程作业。

在explicit equation文件夹中，运行make“expl cit . out”来构建可执行文件。
在隐式方程文件夹中，运行make“implcit . out”构建可执行文件。

使用“-n”修改网格数，使用“-dt”修改时间步长。

使用“-restart”通过HDF5重启设备。

如果想直接运行，可以在生成可执行文件后直接执行“bsub<script.sh”。

输出在“emplicit.h5”或“implicit.h5”中，dx和dt的参数记录在“emplicit.dat”或“implicit.dat”中。

每种方法的答案都在性能分析的最上面，隐式方案会让结果在比较的最下面。

由于在写报告时的代码结果准确性较低，报告中会存在一些与代码不符合的错误。请酌情参考。最新代码具有较高的准确性（1e-6量级），但报告未随之更新。

另外，如果你对这个代码有什么问题，你可以发电子邮件给我，我的电子邮件地址是francis97@foxmail.com（或者南方科技大学的学校邮箱：12132391）。

English version of readme

This program is written in C language, calling PETSc function library, and is a parallel example of defensive programming.

It is the final programming project of SUSTech HPC.

In the explicit equation folder, run make "explcit.out" to build the executable file.
In the implicit equation folder, run make "implcit.out" to build the executable file.

Use "-n " to modify the number of grids and "-dt " to modify the time step.

Use "-restart " to restart facility using HDF5.

If you want to run it directly, you can directly execute "bsub<script.sh" after the executable file is generated.

The output is in "emplicit.h5" or "implicit.h5".

The answer of each method is at the top of the performance analysis, and the implicit scheme will make the results at the bottom of the comparison.

There are some errors in the report due to the low accuracy of the results of the code submitted. Please refer to it as appropriate. The latest code has sufficient accuracy, but it does not match the report.

Also if you have some problem about this code, you can send an e-mail to me, my e-mail address is francis97@foxmail.com.
