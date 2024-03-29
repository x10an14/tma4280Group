#ifndef _EXERCISE_6_
#define _EXERCISE_6_

#include <mpi.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

//! \brief value of pi
#define M_PI 3.14159265358979323846

//The necessary prototyping of the fourier sine transform functions
void fst_(double *v, int *n, double *w, int *nn);
void fstinv_(double *v, int *n, double *w, int *nn);

typedef struct{
	double *data;
	int len;
} vector_t;
typedef vector_t *Vector;

typedef struct{
	double **data;
	Vector as_vec;
	Vector *col;
	MPI_Comm *comm;
	int comm_size;
	int comm_rank;
	int rows;
	int cols;
} matrix_t;
typedef matrix_t *Matrix;

// Function prototypes:
Vector createVector(int);
Matrix createMatrix(int, int);
void freeVector(Vector);
void freeMatrix(Matrix);
void splitVector(int , int , int , int **, int **, int **, int **);
int getMaxThreads();
double WallTime();
void sendArrange(double *, double *, int , int , int *, int , int *);
void recvArrange(double *recvbuf, double *outbuf, int collength, int );
void packTransp(Matrix, Matrix , int *, int *, int );
void fillWithNaturalNumbers(Matrix , int , int *, int );
void fillWithConst(Matrix, double );
void callFourier(Matrix, Matrix );
void callFourierInvrs(Matrix , Matrix);
void unpackTransp(Matrix, Matrix);
void printDoubleMatrix(double **, int, int);
void printDoubleVector(double *, int);
void printIntVector(int *, int);
double linearAverage();
double exactSolAppB(int, int, double);
void fillWithAppB(Matrix, int, int *, double);

#endif
