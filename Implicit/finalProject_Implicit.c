static char help[] = "Solves a tridiagonal linear system.\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <math.h>

#define pi acos(-1)

int main(int argc,char **args)
{
  Vec            x, b, u;          /* A*x=b+u */
  Mat            A;                /* linear system matrix */
  KSP            ksp;              /* linear solver context */
  PC             pc;               /* preconditioner context */
  PetscReal      norm = 0.0, normt = 1.0, tol = 1000.*PETSC_MACHINE_EPSILON;  /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       i, col[3], its, rstart, rend, nlocal, rank;
  PetscInt       n = 128;    /*这是将区域分成n块*/
  PetscReal      dx = 1/n, dt = 0.00001;    /*空间步长和时间步长*/
  PetscReal      p = 1.0, c = 1.0, k = 1.0;    /*设置初始的条件参数*/
  PetscReal      te = k/p/c, alpha = te*dt*n*n;      /*通过dt和dx求解alpha，方便后续计算*/
  PetscScalar    zero = 0.0, value[3];


  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);

  alpha = te*dt*n*n;

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x,&nlocal);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);


  if (!rstart) 
  {
    rstart = 1;
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1+2.0*alpha; value[1] = -alpha;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  if (rend == n) 
  {
    rend = n-1;
    i    = n-1; col[0] = n-2; col[1] = n-1; value[0] = -alpha; value[1] = 1+2.0*alpha;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  value[0] = -alpha; value[1] = 1+2.0*alpha; value[2] = -alpha;
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Assemble the matrix */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


  ierr = VecSet(b,zero);CHKERRQ(ierr);
  if(rank == 0){
    for(int i = 0; i < n; i++){
      PetscReal u0;
      u0 = exp((i+0.5)*dx);
      ierr = VecSetValues(b, 1, &i, &u0, INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  
  
  ierr = VecSet(u,zero);CHKERRQ(ierr);
  if(rank == 0){
    for(int i = 0; i < n; i++){
      PetscReal inp;
       inp = dt*sin(pi*(i+0.5)*dx);
      ierr = VecSetValues(u, 1, &i, &inp, INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  
  ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
  
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  while(PetscAbsReal(norm-normt)>tol){
    normt= norm;

    ierr = VecAXPY(b,1.0,u);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
      
    ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    ierr = VecScale(x,(PetscScalar)1.0/norm);CHKERRQ(ierr);
    ierr = VecCopy(x,b);CHKERRQ(ierr);
  }
  
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);CHKERRQ(ierr);

  ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"solution = %f\n",norm);CHKERRQ(ierr);
 
  ierr = VecDestroy(&x);CHKERRQ(ierr); 
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}

// EOF