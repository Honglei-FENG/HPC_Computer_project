static char help[] = "Solves a tridiagonal linear system.\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <math.h>

#define pi acos(-1)    /*定义pi（3.1415926...）的值*/

int main(int argc,char **args)
{
  Vec            x, z, b;    /*设置所需向量*/
  Mat            A;    /*设置所需矩阵*/
  PetscErrorCode ierr;    /*检查错误信息*/
  PetscInt       i, ii, col[3], rstart, rend, nlocal, rank;
  /*其中i,ii是矩阵和向量的角标，col是三对角矩阵参数的位置，rstart和rend均为设置矩阵时需要的参数，nlocal和rank为程序并行化所需参数*/
  PetscInt       n = 128, start = 0, end;    /*这是将区域分成n块，start是起始边界，end是终止边界*/
  PetscReal      dx, dt = 0.00003, t = 0.0;    /*dx是空间步长，dt是时间步长，t是已经走过的时间*/
  PetscReal      p = 1.0, c = 1.0, k = 1.0;    /*设置初始的条件参数*/
  PetscReal      te = k/p/c, alpha;    /*通过dt和dx求解alpha，方便后续计算*/
  PetscScalar    zero = 0.0, value[3], u0 = 0.0;    /*value是设置三对角矩阵的参数，u0是初始条件*/
  PetscViewer    pv;


  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;    /*初始化Petsc*/
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);    /*从命令行读取n的值（若有）*/
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);    /*从命令行读取dt的值（若有）*/
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);    /*设置并行MPI参数*/
  ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);    /*将n的值打印出来，方便阅读输出文件时参考*/

  dx   = 1.0/n;    /*计算出每一小格的长度*/
  alpha= te*dt*n*n;    /*计算出CFL的值（显式格式中，CFL不能大于0.5）*/
  end  = n;    /*更新end的值*/
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dx = %f\n",dx);CHKERRQ(ierr);    /*将dx的值打印出来，方便阅读输出文件时参考*/
  ierr = PetscPrintf(PETSC_COMM_WORLD,"alpha = %f\n",alpha);CHKERRQ(ierr);    /*将alpha的值打印出来，方便阅读输出文件时参考*/

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);    /*创建一个并行空间*/
  ierr = VecSetSizes(x,PETSC_DECIDE,n+1);CHKERRQ(ierr);    /*创建一个长度n+1的矩阵*/
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);    /*从选项数据库中配置向量*/
  ierr = VecDuplicate(x,&z);CHKERRQ(ierr);    /*将x的格式赋给z*/
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);    /*将x的格式赋给b*/

  ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRQ(ierr);    /*设置并行x的起始终止点*/
  ierr = VecGetLocalSize(x,&nlocal);CHKERRQ(ierr);    /*设置并行x的位置点*/

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);    /*在并行空间创建一个矩阵*/
  ierr = MatSetSizes(A,nlocal,nlocal,n+1,n+1);CHKERRQ(ierr);    /*设置矩阵的行数和列数*/
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);    /*从选项数据库中配置矩阵*/
  ierr = MatSetUp(A);CHKERRQ(ierr);    /*开始建立矩阵*/

  if (!rstart)     /*若rstart为0时，即为首行*/
  {
    rstart = 1;    /*将rstart设为1*/
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1-2.0*alpha; value[1] = alpha;    /*设置要用到的参数*/
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);    /*设置三对角矩阵的第一行*/
  }
  
  if (rend == n+1)     /*最后一行*/
  {
    rend = n;    /*将rend设为n*/
    i    = n; col[0] = n-1; col[1] = n; value[0] = alpha; value[1] = 1-2.0*alpha;    /*设置要用到的参数*/
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);    /*设置三对角矩阵的最后一行*/
  }

  value[0] = alpha; value[1] = 1-2.0*alpha; value[2] = alpha;    /*设置三对角矩阵除首尾两行外的其余行的三个值*/
  for (i=rstart; i<rend; i++)    /*除首尾两行外的行*/
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;    /*设置要用到的参数*/
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);    /*设置三对角矩阵的三对角值*/
  }

  /* Assemble the matrix */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);    /*通知其余并行块将矩阵统一*/
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);    /*结束通知*/
  ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    /*打印矩阵，检查是否出错*/
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dx = %f\n",dx);CHKERRQ(ierr);    /*将dx的值打印出来，方便阅读输出文件时参考*/

  ierr = VecSet(z,zero);CHKERRQ(ierr);    /*设置初始向量z*/
  if(rank == 0)    /*开始设置初始条件*/
  {
    for(ii = 1; ii < n; ii++)    /*除首尾两个点外的其余点*/
    {
      u0   = exp(ii*dx);    /*根据当前位置来获取初始值*/
	    ierr = VecSetValues(z, 1, &ii, &u0, INSERT_VALUES);CHKERRQ(ierr);    /*将向量的对应位置的值进行修改*/
    }
  }
  
  ierr = VecAssemblyBegin(z);CHKERRQ(ierr);    /*通知其余并行块将向量统一*/
  ierr = VecAssemblyEnd(z);CHKERRQ(ierr);    /*结束通知*/

  
  ierr = VecSet(b,zero);CHKERRQ(ierr);    /*设置初始向量b*/
  if(rank == 0){    /*开始设置初始条件*/
    for(ii = 1; ii < n+1; ii++){    /*除首尾两个点外的其余点*/
    u0   = dt*sin(ii*dx*pi);    /*根据当前位置来获取传热值*/
	  ierr = VecSetValues(b, 1, &ii, &u0, INSERT_VALUES);CHKERRQ(ierr);    /*将向量的对应位置的值进行修改*/
    }
  }
  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);    /*通知其余并行块将向量统一*/
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);    /*结束通知*/
  
  while(PetscAbsReal(t)<=2.0){    /*计算0-2时间内的传播*/

     t += dt;    /*时间向前走*/
     
     ierr = MatMult(A,z,x);CHKERRQ(ierr);    /*求得下一个时刻的温度热量值*/
     ierr = VecAXPY(x,1.0,b);CHKERRQ(ierr);    /*将求得的值加上在时间步长内的*/

     ierr = VecSetValues(x, 1, &start, &zero, INSERT_VALUES);CHKERRQ(ierr);    /*设置边界条件*/
     ierr = VecSetValues(x, 1, &end, &zero, INSERT_VALUES);CHKERRQ(ierr);    /*设置边界条件*/

     ierr = VecAssemblyBegin(x);CHKERRQ(ierr);    /*统一向量更新*/
     ierr = VecAssemblyEnd(x);CHKERRQ(ierr);    /*结束更新*/
	   
     ierr = VecCopy(x,z);CHKERRQ(ierr);    /*将x的值赋给z*/
  }
  
  ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    /*打印向量，获得结束时显式方法的值*/

  /*Viewer to output in HDF5 format*/
  
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&pv);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"explicit.dat",&pv);CHKERRQ(ierr);
  ierr = VecView(z, pv);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&pv);CHKERRQ(ierr);
  
  ierr = VecDestroy(&x);CHKERRQ(ierr);    /*关闭向量x*/
  ierr = VecDestroy(&z);CHKERRQ(ierr);    /*关闭向量z*/
  ierr = VecDestroy(&b);CHKERRQ(ierr);    /*关闭向量b*/
  ierr = MatDestroy(&A);CHKERRQ(ierr);    /*关闭矩阵A*/

  ierr = PetscFinalize();    /*结束并行*/
  return ierr;
}

// EOF