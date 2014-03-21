#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "ex6.h"

/* This function does NOT allocate the vector memory data! */
Vector createVector(int len){
	Vector result = (Vector) calloc(1, sizeof(vector_t));
	result->len = len;
	result->data = NULL;

	return result;
}

/*While this function does allocate memory for the matrix data*/
Matrix createMatrix(int rows, int cols){
	Matrix result = (Matrix) calloc(1, sizeof(matrix_t));

	result->cols = cols;
	result->rows = rows;
	result->comm_size = 1;
	result->comm_rank = 0;

	result->as_vec = createVector(rows*cols);
	result->data = (double **) calloc(cols, sizeof(double*));
	result->data[0] = (double *) calloc(rows*cols, sizeof(double));
	result->as_vec->data = result->data[0];

	//Correctly set indices of a matrix with elements linearly stored in memory (c-standard, not fortran-standard regarding rows/cols)
	for (int i = 1; i < cols; ++i){
		result->data[i] = result->data[i-1] + rows;
	}

	//Make each row of data (which is linear in memory) accessible as a vector
	result->col = (Vector*) calloc(cols, sizeof(Vector));
	for (int i = 0; i < cols; ++i){
		result->col[i] = createVector(rows);
		result->col[i]->data = result->data[i];
	}

	return result;
}

//Copied from common.c, minor edits.
void splitVector(int globLen, int size, int** len, int** displ){
	*len = calloc(size,sizeof(int));
	*displ = calloc(size,sizeof(int));

	for (int i=0; i < size; ++i){
		(*len)[i] = globLen/size;

		if (globLen % size && i >= (size - globLen % size)){
			(*len)[i]++;
		}

		if (i < size-1){
			(*displ)[i+1] = (*displ)[i]+(*len)[i];
		}
	}
}

void printDoubleMatrix(double **ptr, int height, int width){
	for (int i = 0; i < height; ++i){
		printf("[%.2f", ptr[i][0]);
		for (int j = 1; j < width; ++j){
			printf(",\t%.2f", ptr[i][j]);
		}
		printf("\t]\n");
	}
}

void printDoubleVector(double *ptr, int length){
	printf("[%.2f", ptr[0]);
	for (int i = 1; i < length; ++i){
		printf(",\t%.2f", ptr[i]);
	}
	printf("]\n");
}

void printIntMatrix(int **ptr, int height, int width){
	for (int i = 0; i < height; ++i){
		printf("[%i", ptr[i][0]);
		for (int j = 1; j < width; ++j){
			printf(",\t%i", ptr[i][j]);
		}
		printf("\t]\n");
	}
}

void printIntVector(int *ptr, int length){
	printf("[%i", ptr[0]);
	for (int i = 1; i < length; ++i){
		printf(",\t%i", ptr[i]);
	}
	printf("]\n");
}

/* Arranges the sendbuffer properly before sending
 * Assumes the columns are arranged continually in the vector
 */
void sendArrange(double *sendbuf, double *vector, int collength, int colcnt, int *sizearr, int sizearrlength){
	int elements, counter, process, destoffset;

	for(int i = 0; i < colcnt; i++){	// One pass per col
		process = 0;
		elements = sizearr[process];	// Elements to go to the first process
		counter  = 0;					// Counter for how many elements have been moved
		destoffset = 0;					// Destination offset for elements (i.e. which process is this being sent to)

		for(int j = 0; j < collength; j++){	// Iterate through coloumn
			if((counter + 1 >= elements) && (process <= sizearrlength - 1)){
				counter = 0;
				process++;
				destoffset += elements*colcnt;
				elements = sizearr[process];
			}

			sendbuf[i*elements + counter + destoffset] = vector[i*collength + j];
			counter++;
		}
	}
}

int getMaxThreads(){
	#ifdef HAVE_OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
}

double WallTime(){
	/*#ifdef defined(HAVE_OPENMP)
		return omp_get_wtime();
	#else*/
		return MPI_Wtime();
	//#endif
}

void freeVector(Vector inpt){
	free(inpt->data);
	free(inpt);
}

void freeMatrix(Matrix inpt){
	for (int i = 0; i < inpt->cols; ++i){
		free(inpt->col[i]);
	}
	free(inpt->col);
	free(inpt->as_vec);
	free(inpt->data[0]);
	free(inpt->data);
	free(inpt);
}

void VariableInitz(int n, double *h, int *globColLen, int *tempMatSz){
	*h = 1.0/n;
	*h *= *h;
	*globColLen = n-1;
	*tempMatSz = n*4;
}

#define TEST 0
int print = 0;

int main(int argc, char *argv[]){
	Matrix	/*The "work"-matrix*/matrix, \
			/*The transposed version of the matrix*/transpMat;
	Vector	/*The diagonal matrix*/diagMat, \
			/*The fourier-transform temp-store-matrix*/tempMat;

	//Lists for holding how many rows per process (due to MPI division of labour), AKA MPI variables...
	int *size, *displacement, rank = 0, mpiSize = 1, acquired;

	//Global variables
	double h, time;
	int globColLen, n;

	//Process specific variables
	int tempMatSz, locMatSz, procColAmnt;

	/*		General MPI startup/setup		*/
	#ifdef HAVE_OPENMP
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &acquired);
	#else
	MPI_Init(&argc, &argv);
	#endif
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_dup(MPI_COMM_WORLD, &WorldComm);

	//Check for correct commandline argument
	//Remove the section comments in the if-test when we want to test "live"/production/release version.
	if (argc < 2 /*|| atoi(argv[1]) < 4 || atoi(argv[1])%2 != 0*/){
		printf("Need a problem size! And it should be a power of 2 greater than or equal to 4!\n");
		MPI_Finalize();
		return -1;
	}

	/*			Initializing variables		*/
	n = atoi(argv[1]);
	VariableInitz(n, &h, &globColLen, &tempMatSz);
	splitVector(globColLen, mpiSize, &size, &displacement);
	procColAmnt = size[rank];
	locMatSz = globColLen*procColAmnt;

	if(rank == TEST && print){
		printf("n: %i\nglobColLen: %i\nmpiSize: %i\n\n", n, globColLen, mpiSize);
		printf("Size vector:\n");
		printIntVector(size, mpiSize);
		printf("Displacement vector:\n");
		printIntVector(displacement, mpiSize);
		printf("\nSize of rank %i matrix: %ix%i\n", TEST, n-1, size[TEST]);
	}

	/*			Initializing structures		*/
	tempMat = createVector(tempMatSz);
	diagMat = createVector(globColLen);
	matrix = createMatrix(globColLen, procColAmnt);
	transpMat = createMatrix(globColLen, procColAmnt);
	diagMat->data = (double*) malloc(globColLen*sizeof(double));
	//Transpose temporary buffer
	double *sendbuf = calloc(globColLen*procColAmnt, sizeof(double));

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < diagMat->len; ++i){
		//Filling the diagonal matrix
		diagMat->data[i] = (double) 2.0*(1.0-cos(i + 1.0)*M_PI/(double)n);
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < locMatSz; ++i){
		//Filling up the work-matrix
		matrix->as_vec->data[i] = h;
	}

	time = WallTime();

	if(rank == TEST && print){
		printf("Diagonal matrix for rank == %i:\n", TEST);
		printDoubleVector(diagMat->data, globColLen);
	}

	/*#pragma omp parallel for schedule(static)
	for (int i = 0; i < procColAmnt; ++i){
		//Implementation of the first fst_() call
		tempMat->data = (double*) calloc(tempMatSz, sizeof(double));
		fst_(matrix->data[i], &globColLen, tempMat->data, &tempMatSz);
	}*/

	/*		Implementation of the first transpose			*/
	Matrix transpTest = createMatrix(3, 3);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			transpTest->data[j][i] = (double) (j+1) + 3*i;
		}
	}

	if (rank == 0){
		printf("Before transpose:\n");
		printDoubleMatrix(transpTest->data, 3, 3);
	}
	double *tst = calloc(9, sizeof(double));
	sendArrange(tst, transpTest->data[0], 3, 3, size, mpiSize);
	MPI_Alltoallv(&tst, size, displacement, MPI_DOUBLE, transpTest->data[0], size, displacement, MPI_DOUBLE, WorldComm);

	if(rank == 0){
		printf("\nAfter transpose:\n");
		printDoubleMatrix(transpTest->data, 3, 3);
		printf("\n");
	}

	/*sendArrange(sendbuf, matrix->data[0], globColLen,procColAmnt, size, mpiSize); // Arrange send buffer
	MPI_Alltoallv(&sendbuf, size, displacement, MPI_DOUBLE, matrix->data[0], size, displacement, MPI_DOUBLE, WorldComm);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procColAmnt; ++i){
		//Implementation of the first fstinv_() call
		tempMat->data = (double*) calloc(tempMatSz, sizeof(double));
		fstinv_(transpMat->data[i], &globColLen, tempMat->data, &tempMatSz);
	}*/

	/*		Implementation of the "tensor" operation		*/
	//Which for-loop level should get the open mp pragma?
	/*#pragma omp parallel for schedule(static)
	for (int i = 0; i < procColAmnt; ++i){
		for (int j = 0; j < globColLen; ++j){
			transpMat->data[i][j] /= diagMat->data[j] + diagMat->data[i];
		}
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procColAmnt; ++i){
		//Implementation of the second fst_() call
		tempMat->data = (double*) calloc(tempMatSz, sizeof(double));
	//fst_(transpMat->data[i], &globColLen, tempMat->data, &tempMatSz);
	}*/

	/*		Implementation of the second transpose			*/
	//Arrange send buffer(using same buffer as last time)
	/*sendArrange(sendbuf, matrix->data[0], globColLen,procColAmnt, size, mpiSize);
	MPI_Alltoallv(&sendbuf, size, displacement, MPI_DOUBLE, matrix->data[0], size, displacement, MPI_DOUBLE, WorldComm);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procColAmnt; ++i){
		//Implementation of the second fstinv_() call
		tempMat->data = (double*) calloc(tempMatSz, sizeof(double));
	//fstinv_(matrix->data[i], &globColLen, tempMat->data, &tempMatSz);
	}*/

	/*		Print time? (not yet implemented)				*/
	if(rank == 0){
		time = WallTime() - time;
		printf("t: %g\n", time);
	}

	/*		Closing up and freeing variables				*/
	freeMatrix(matrix);
	freeMatrix(transpMat);
	freeVector(tempMat);
	freeVector(diagMat);
	if(rank == 0){
		MPI_Comm_free(&WorldComm);
	}
	MPI_Finalize();

	return 0;
}
