#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "blaslapack.h"
#include "common.h"
#include "solvers.h"
#include "poissoncommon.h"

double alpha=0.0;

double exact(double x)
{
  return x*(pow(x,5)-1.0);
}

double source(double x)
{
  return -30*pow(x,4);
}

int GaussJacobiPoisson1D(Vector u, double tol, int maxit)
{
  int it=0, i;
  Vector b = cloneVector(u);
  Vector e = cloneVector(u);
  copyVector(b, u);
  fillVector(u, 0.0);
  double max = tol+1;
  while (max > tol && ++it < maxit) {
    copyVector(e, u);
    collectVector(e);
    copyVector(u, b);
#pragma omp parallel for schedule(static)
    for (i=1;i<e->len-1;++i) {
      u->data[i] += e->data[i-1];
      u->data[i] += e->data[i+1];
      u->data[i] /= (2.0+alpha);
    }
    axpy(e, u, -1.0);
    e->data[0] = e->data[e->len-1] = 0.0;
    max = maxNorm(e);
  }
  freeVector(b);
  freeVector(e);

  return it;
}

int GaussSeidelPoisson1Drb(Vector u, double tol, int maxit)
{
  int it=0, i, j;
  double max = tol+1;
  Vector b = cloneVector(u);
  Vector r = cloneVector(u);
  Vector v = cloneVector(u);
  copyVector(b, u);
  fillVector(u, 0.0);
  while (max > tol && ++it < maxit) {
    copyVector(v, u);
    copyVector(u, b);
    collectVector(v);
    for (j=0;j<2;++j) {
#pragma omp parallel for schedule(static)
      for (i=1+j;i<r->len-1;i+=2) {
        u->data[i] += v->data[i-1];
        u->data[i] += v->data[i+1];
        r->data[i] = u->data[i]-(2.0+alpha)*v->data[i];
        u->data[i] /= (2.0+alpha);
        v->data[i] = u->data[i];
      }
      if (j == 0)
        collectVector(v);
    }
    max = maxNorm(r);
  }
  freeVector(b);
  freeVector(r);
  freeVector(v);

  return it;
}

void Poisson1D(Vector u, Vector v)
{
  int i;
  collectVector(v);
#pragma omp parallel for schedule(static)
  for (i=1;i<u->len-1;++i) {
    u->data[i] = (2.0+alpha)*v->data[i];
    u->data[i] -= v->data[i-1];
    u->data[i] -= v->data[i+1];
  }
  u->data[0] = u->data[u->len-1] = 0.0;
}

void Poisson1Dnoborder(Vector u, Vector v)
{
  int i;
#pragma omp parallel for schedule(static)
  for (i=0;i<u->len;++i) {
    u->data[i] = (2.0+alpha)*v->data[i];
    if (i > 0)
      u->data[i] -= v->data[i-1];
    if (i < u->len-1)
      u->data[i] -= v->data[i+1];
  }
}

//! \brief A struct describing a subdomain in additive schwarz
typedef struct {
  Matrix A;             //!< Poisson operator
  int Afactored;        //!< Whether or not the operator is factorized
  int size;             //!< Size of subdomain
  int from_disp;        //!< Left displacement when creating subdomain vector
  int to_source_disp;   //!< Left displacement from source recreating full vector
  int to_dest_disp;     //!< Left displacement in destination recreating full vector
  int to_dest_size;     //!< Datas to copy to destinatino recreating full vector
  Matrix Q;             //!< Eigenvectors of Poisson operator
  Vector lambda;        //!< Eigenvalues of Poisson operator
} SchwarzSubdomain;

SchwarzSubdomain subdomain;

Vector collectBeforePre(Vector u)
{
  collectVector(u);
  Vector result=createVector(subdomain.size);
  dcopy(&subdomain.size, u->data+subdomain.from_disp,
        &u->stride, result->data, &result->stride);

  return result;
}

void collectAfterPre(Vector u, const Vector v)
{
  int source, dest;

  dcopy(&subdomain.to_dest_size, v->data+subdomain.to_source_disp,
        &v->stride, u->data+subdomain.to_dest_disp, &u->stride);

  // west
  MPI_Cart_shift(*u->comm, 0,   -1, &source, &dest);
  MPI_Sendrecv(v->data,        1, MPI_DOUBLE, dest,   0,
               u->data, 1, MPI_DOUBLE, source, 0, *u->comm, MPI_STATUS_IGNORE);
  if (source > -1)
    u->data[u->len-2] += u->data[0];

  // east
  MPI_Cart_shift(*u->comm,  0,   1, &source, &dest);
  MPI_Sendrecv(v->data+v->len-1, 1, MPI_DOUBLE, dest,   1,
               u->data,          1, MPI_DOUBLE, source, 1, *u->comm, MPI_STATUS_IGNORE);
  if (source > -1)
    u->data[1] += u->data[0];

  u->data[0] = u->data[u->len-1] = 0.0;
}

void DiagonalizationPoisson1D(Vector u, const Vector lambda, const Matrix Q)
{
  Vector btilde = createVector(u->len);
  int i;
  MxV(btilde, Q, u, 1.0, 0.0, 'T');
  for (i=0;i<btilde->len;++i)
    btilde->data[i] /= (lambda->data[i]+alpha);
  MxV(u, Q, btilde, 1.0, 0.0, 'N');
  freeVector(btilde);
}

void Poisson1DPre(Vector u, Vector v)
{
  Vector tmp = collectBeforePre(v);
  if (subdomain.A) {
    llsolve(subdomain.A, tmp, subdomain.Afactored);
    subdomain.Afactored = 1;
  } else if (subdomain.Q)
    DiagonalizationPoisson1D(tmp, subdomain.lambda, subdomain.Q);
  else
    cgMatrixFree(Poisson1Dnoborder, tmp, 1e-10);
  collectAfterPre(u, tmp);
  freeVector(tmp);
}

int main(int argc, char** argv)
{
  int i, j, N, flag;
  Vector b, grid, e, lambda=NULL;
  double time, sum, h, tol=1e-6;
  int rank, size;
  int mpi_top_coords;
  int mpi_top_sizes;

  init_app(argc, argv, &rank, &size);

  if (argc < 3) {
    printf("need two parameters, N and flag [alpha] [tolerance]\n");
    printf(" - N is the problem size (in each direction\n");
    printf(" - flag = 1  -> Matrix-free Gauss-Jacobi iterations\n");
    printf(" - flag = 2  -> Matrix-free red-black Gauss-Seidel iterations\n");
    printf(" - flag = 3  -> Matrix-free CG iterations\n");
    printf(" - flag = 4  -> Matrix-free additive schwarz+Cholesky CG iterations\n");
    printf(" - flag = 5  -> Matrix-free additive schwarz+CG preconditioned CG iterations\n");
    printf(" - flag = 6  -> Matrix-free additive schwarz+diagonalization preconditioned CG iterations\n");
    printf(" - alpha is the Helmholtz scaling factor\n");
    printf(" - tolerance is the residual error tolerance in the iterative scheme\n");
    return 1;
  }
  N=atoi(argv[1]);
  flag=atoi(argv[2]);
  if (argc > 3)
    alpha = atof(argv[3]);
  if (argc > 4)
    tol=atof(argv[4]);

  if (N < 0) {
    if (rank == 0)
      printf("invalid problem size given\n");
    close_app();
    return 2;
  }

  if (flag < 0 || flag > 6) {
    if (rank == 0)
      printf("invalid flag given\n");
    close_app();
    return 3;
  }

  if (flag == 2 && (N-1)%2 != 0 && ((N-1)/size) % 2 != 0) {
    if (rank == 0)
      printf("need an even size (per process) for red-black iterations\n");
    close_app();
    return 4;
  }

  // setup topology
  mpi_top_coords = 0;
  mpi_top_sizes = 0;
  MPI_Dims_create(size, 1, &mpi_top_sizes);
  int periodic = 0;
  MPI_Comm comm;
  MPI_Cart_create(MPI_COMM_WORLD, 1, &mpi_top_sizes, &periodic, 0, &comm);
  MPI_Cart_coords(comm, rank, 1, &mpi_top_coords);

  b = createVectorMPI(N+1, &comm, 1, 1);
  e = createVectorMPI(N+1, &comm, 1, 1);

  // setup subdomain
  subdomain.A = subdomain.Q = NULL;
  subdomain.lambda = NULL;
  if (mpi_top_coords == 0) {
    subdomain.size = subdomain.to_dest_size = b->len-1;
    subdomain.from_disp = 1;
    subdomain.to_source_disp = 0;
    subdomain.to_dest_disp = 1;
  } else if (mpi_top_coords == mpi_top_sizes-1) {
    subdomain.size = b->len-1;
    subdomain.to_dest_size = b->len-2;
    subdomain.from_disp = 0;
    subdomain.to_source_disp = subdomain.to_dest_disp = 1;
  } else {
    subdomain.size = subdomain.to_dest_size = b->len;
    subdomain.from_disp = subdomain.to_source_disp = subdomain.to_dest_disp = 0;
  }

  grid = equidistantMesh(0.0, 1.0, N);
  h = 1.0/N;

  evalMeshDispl(b, grid, source);
  scaleVector(b, pow(h, 2));
  evalMeshDispl(e, grid, exact);
  axpy(b, e, alpha);
  b->data[0] = b->data[b->len-1] = 0.0;

  if (flag == 4) {
    subdomain.A = createMatrix(subdomain.size, subdomain.size);
    subdomain.Afactored = 0;
    diag(subdomain.A, -1, -1.0);
    diag(subdomain.A, 0, 2.0+alpha);
    diag(subdomain.A, 1, -1.0);
  }
  if (flag == 6) {
    subdomain.lambda = generateEigenValuesP1D(subdomain.size);
    subdomain.Q = generateEigenMatrixP1D(subdomain.size);
  }

  int its=-1;
  char method[128];
  time = WallTime();
  if (flag == 1) {
    its=GaussJacobiPoisson1D(b, tol, 1000000);
    sprintf(method,"Gauss-Jacobi");
  }
  if (flag == 2) {
    its=GaussSeidelPoisson1Drb(b, tol, 1000000);
    sprintf(method,"Gauss-Seidel");
  }
  if (flag == 3) {
    its=cgMatrixFree(Poisson1D, b, tol);
    sprintf(method,"CG");
  }
  if (flag == 4 || flag == 5 || flag == 6) {
    its=pcgMatrixFree(Poisson1D, Poisson1DPre, b, tol);
    sprintf(method,"PCG");
  }
  if (rank == 0) {
    printf("%s used %i iterations\n", method, its);
    printf("elapsed: %f\n", WallTime()-time);
  }

  evalMeshDispl(e, grid, exact);
  axpy(b,e,-1.0);
  b->data[0] = b->data[b->len-1] = 0.0;

  h = maxNorm(b);
  if (rank == 0)
    printf("max error: %e\n", h);
  
  freeVector(grid);
  freeVector(b);
  freeVector(e);
  if (lambda)
    freeVector(lambda);
  if (subdomain.A)
    freeMatrix(subdomain.A);
  if (subdomain.Q)
    freeMatrix(subdomain.Q);
  if (subdomain.lambda)
    freeVector(subdomain.lambda);

  MPI_Comm_free(&comm);

  close_app();
  return 0;
}
