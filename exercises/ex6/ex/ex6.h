#ifndef _EXERCISE_6_
#define _EXERCISE_6_

#include <mpi.h>

//All the existing .h files I could think of with potential use.
//#include "../../../examples/common/common.h"
//#include "../../../examples/common/blaslapack.h"
//#include "../../../examples/poisson/solvers.h"
//#include "../../../examples/poisson/possioncommon.h"

//! \brief value of pi
#define M_PI 3.14159265358979323846

/*Whatever we want specific to our solution/implementation.*/

#ifdef HAVE_MPI
MPI_Comm WorldComm;
MPI_Comm SelfComm;
#endif

typedef struct{
	double *data;
	int len;
} vector_t;

typedef vector_t *Vector;

typedef struct{
	double **data;
	Vector as_vec;
	Vector *row;
	MPI_Comm *comm;
	int comm_size;
	int comm_rank;
	int rows;
	int cols;
} matrix_t;

typedef matrix_t *Matrix;

#endif
