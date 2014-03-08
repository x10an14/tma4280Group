#include "ex6.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

Vector createVector(int len){
	Vector result = (Vector) calloc(1, sizeof(vector_t));
	result->len = len;
	result->data = NULL;

	return result;
}

Matrix createMatrix(int rows, int cols){
	Matrix result = (Matrix) calloc(1, sizeof(matrix_t));

	result->cols = cols;
	result->rows = rows;
	result->comm_size = 1;
	result->comm_rank = 0;

	result->as_vec = createVector(rows*cols);
	result->data = (double **) calloc(rows, sizeof(double*));
	result->data[0] = (double *) calloc(rows*cols, sizeof(double));
	result->as_vec->data = result->data[0];

	//Correctly set indices of a matrix with elements linearly stored in memory (c-standard, not fortran-standard regarding rows/cols)
	for (int i = 1; i < rows; ++i){
		result->data[i] = result->data[i-1] + cols;
	}

	//Make each row of data (which is linear in memory) accessible as a vector
	result->row = (Vector*) calloc(rows, sizeof(Vector));
	for (int i = 0; i < rows; ++i){
		result->row[i] = createVector(cols);
		result->row[i]->data = result->data[i];
	}

	return result;
}

//Copied from common.c, minor edits.
void splitVector(int globLen, int size, int** len, int** displ){
	*len = calloc(size,sizeof(int));
	*displ = calloc(size,sizeof(int));

	for (int i=0; i < size; ++i){
		(*len)[i] = globLen/size;

		if (globLen % size && i >= (size - (globLen % size))){
			(*len)[i]++;
		}

		if (i < size-1){
			(*displ)[i+1] = (*displ)[i]+(*len)[i];
		}
	}
}

/*DO NOT USE THE COMMONS LIBRARY!
*IF ANYTHING IN THE COMMONS LIBRARY IS OF USE, COPY IT OVER.
*THE COMMONS LIBRARY HAS TOO MUCH "UNNECESSARY DATA/JUNK"
*/

int main(int argc, char *argv[]){
	Matrix	/*The "work"-matrix*/matrix, \
			/*The transposed version of the matrix*/transpMat;
	Vector	/*The diagonal matrix*/diagMat, \
			/*The fourier-transform temp-store-matrix*/f_TempMat;

	//Lists keeping division of MPI labour-division of processes
	int *size, *displacement;

	//Just general variables used in each process.
	int matrixSize, matrixTempSize, n, rank, mpiSize, acquired;
	double h = 1.0;

	//General MPI startup/setup
	#ifdef HAVE_OPENMP
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &acquired);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#else
	MPI_init(&argc, &argv);
	#endif
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_dup(MPI_COMM_WORLD, &WorldComm);
	MPI_Comm_dup(MPI_COMM_SELF, &SelfComm);

	//Check for comandline argument
	if (argc < 2 && atoi(argv[1]) <= 4 && atoi(argv[1])%2 != 0){
		printf("Need a problem size! And it should be a power of 2 greater than 4!\n");
		//MPI_Comm_free(&WorldComm);
		//MPI_Comm_free(&SelfComm);
		MPI_Finalize();
		return -1;
	}

	/*			Initializing variables		*/

	n = atoi(argv[1]);
	matrixSize = n-1;
	matrixTempSize = matrixSize*4;
	h /= n;

	//We need to set up the rest of the variables and use the data given to us by splitVector in a smart manner.
	//Do something with splitVector and *size and *displacement
	splitVector(matrixSize, mpiSize, &size, &displacement);

	/*			Initializing structures		*/
	//Use a MPI communicator? We have two defined in the .h file. Maybe make createMatrix() set that for us.
	diagMat = createVector(matrixSize);
	f_TempMat = createVector(matrixTempSize);
	matrix = createMatrix(matrixSize, matrixSize);
	transpMat = createMatrix(matrixSize, matrixSize);
	diagMat->data = (double*) calloc(matrixSize, sizeof(double));
	f_TempMat->data = (double*) malloc(matrixTempSize*sizeof(double));

	/*	The rest of the initialization of structures, making sure that the matrix and diagMat variables only contain what each process needs	*/

	/*		Implementation of the fst_() call			*/

	/*		Implementation of the transpose				*/

	/*		Implementation of the fstinv_() call		*/

	/*		Implementation of the "tensor" operation		*/

	/*	//Namely this:
	*	bt[i][j] = bt[i][j] / (diag[i] + diag[j])
	*/

	/*		Implementation of the fst_() call			*/

	/*		Implementation of the transpose				*/

	/*		Implementation of the fstinv_() call		*/

	/*		Closing up and freeing variables			*/

	MPI_Comm_free(&WorldComm);
	MPI_Comm_free(&SelfComm);
	MPI_Finalize();
	return 0;
}
