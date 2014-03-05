#ifndef _EXERCISE_6_
#define _EXERCISE_6_

//All the existing .h files I could think of with potential use.
#include "../../examples/common/common.h"
//#include "../../examples/common/blaslapack.h"
#include "../../examples/poisson/solvers.h"
#include "../../examples/poisson/possioncommon.h"

/*Whatever we want specific to our solution/implementation.*/

MPI_Comm WorldComm;

typedef struct{
	double *data;
	int len;
} vector_t;

typedef vector_t *Vector;

typedef struct{
	double **data;
	Vector as_vec;
	MPI_Comm *comm;
	int comm_size;
	int comm_rank;
	int rows;
	int cols;
} matrix_t;

typedef matrix_t *Matrix;

#endif
