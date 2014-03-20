#include <stdio.h>

#include "petscksp.h"

int main(int argc, char** argv)
{
  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */
  if( argc < 2 ) {
    printf("need a problem size\n");
    return 1;
  }

  int n  = atoi(argv[1]);
  int m  = n-1;

  // Initialize Petsc
  PetscInitialize(&argc,&argv,0,PETSC_NULL);

  // Create our vector
  Vec b;
  VecCreate(PETSC_COMM_WORLD,&b);
  VecSetSizes(b,PETSC_DECIDE,m*m);
  VecSetFromOptions(b);

  // Create our matrix
  Mat A;
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*m,m*m);
  MatSetUp(A);

  // Create linear solver object
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD,&ksp);

  // setup rhs
  double h    = 1./(double)n;
  for (int j=0; j < m*m; j++)
    VecSetValue(b,j,h*h,INSERT_VALUES);

  // setup matrix
  for (int i=0;i<m*m;++i) {
    MatSetValue(A,i,i,4.f,INSERT_VALUES);
    if (i%m != m-1)
      MatSetValue(A,i,i+1,-1.f,INSERT_VALUES);
    if (i%m)
      MatSetValue(A,i,i-1,-1.f,INSERT_VALUES);
    if (i > m)
      MatSetValue(A,i,i-m,-1.f,INSERT_VALUES);
    if (i < m*(m-1))
    MatSetValue(A,i,i+m,-1.f,INSERT_VALUES);
  }

  // sync processes
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  // solve
  KSPSetType(ksp,"cg");
  KSPSetTolerances(ksp,1.e-10,1.e-10,1.e6,10000);
  KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
  PC pc;
  KSPGetPC(ksp,&pc);
  PCSetFromOptions(pc);
  PCSetUp(pc);
  // setup solver
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
  Vec x;
  VecDuplicate(b,&x);
  KSPSolve(ksp,b,x);

  double val;
  VecNorm(x,NORM_INFINITY,&val);

  printf (" umax = %e \n",val);
  PetscFinalize();
  return 0;
}
