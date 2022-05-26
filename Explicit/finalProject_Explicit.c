static char help[] = "Solves a tridiagonal linear system.\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <math.h>

#define pi acos(-1)

int main(int argc,char **args)
{
  Vec            x, z, b;          
  Mat            A;                /* linear system matrix */
  PetscReal      norm = 0.0, xi = 0.0;  /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       i, ii = 0, col[3], rstart, rend, nlocal, rank;
  PetscInt       n = 128, start = 0, end = n;    /*这是将区域分成n块*/
  PetscReal      dx = 0.0, dt = 0.00001, t = 0.0;    /*空间步长和时间步长*/
  PetscReal      p = 1.0, c = 1.0, k = 1.0;    /*设置初始的条件参数*/
  PetscReal      te = k/p/c, alpha = te*dt*n*n;      /*通过dt和dx求解alpha，方便后续计算*/
  PetscScalar    zero = 0.0, value[3], u0 = 0.0;


  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);

  dx = 1/(PetscReal)n;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dx = %f\n",dx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"alpha = %f\n",alpha);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n+1);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&z);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x,&nlocal);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n+1,n+1);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);


  if (!rstart) 
  {
    rstart = 1;
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1-2.0*alpha; value[1] = alpha;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  if (rend == n+1) 
  {
    rend = n;
    i    = n; col[0] = n-1; col[1] = n; value[0] = alpha; value[1] = 1-2.0*alpha;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  value[0] = alpha; value[1] = 1-2.0*alpha; value[2] = alpha;
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Assemble the matrix */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


  ierr = VecSet(z,zero);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dx = %f\n",dx);CHKERRQ(ierr);
  if(rank == 0)
  {
    for(ii = 1; ii < n; ii++)
    {
      xi = ii*dx;
      u0 = exp(xi);
	    ierr = VecSetValues(z, 1, &ii, &u0, INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  
  
  
  ierr = VecAssemblyBegin(z);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(z);CHKERRQ(ierr);
  ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  xi = 0.0;
  
  ierr = VecSet(b,zero);CHKERRQ(ierr);
  if(rank == 0){
    for(ii = 1; ii < n+1; ii++){
	  PetscReal inp;
    xi = ii*dx*pi;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"sin(x) = %f\n",sin(xi));CHKERRQ(ierr);
    inp = dt*sin(xi);
	  ierr = VecSetValues(b, 1, &ii, &inp, INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  
  while(PetscAbsReal(t)<3){

     t += dt;
     
     ierr = MatMult(A,z,x);CHKERRQ(ierr);
     ierr = VecAXPY(x,1.0,b);CHKERRQ(ierr);

     ierr = VecSetValues(x, 1, &start, &zero, INSERT_VALUES);CHKERRQ(ierr);
     ierr = VecSetValues(x, 1, &end, &zero, INSERT_VALUES);CHKERRQ(ierr);

     ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
     ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
	 
     ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
     ierr = VecScale(x,(PetscScalar)1.0/norm);CHKERRQ(ierr);
	 
     ierr = VecCopy(x,z);CHKERRQ(ierr);
  }
  
  ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
 
 
  ierr = VecDestroy(&x);CHKERRQ(ierr); 
  ierr = VecDestroy(&z);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}

// EOF