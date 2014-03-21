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
 * Assumes the rows are arranged continually in the vector
 */
void sendArrange(double *sendbuf, double *vector, int rowlength, int rowcnt, int *sizearr, int sizearrlength)
{
	int elements, counter, process, destoffset;
	
	for(int i = 0; i < rowcnt; i++){	// One pass per row
		
		process = 0; 
		elements = sizearr[process];	// Elements to go to the first process
		counter  = 0;					// Counter for how many elements have been moved
		destoffset = 0;					// Destination offset for elements (i.e. which process is this being sent to)
		
		for(int j = 0; j < rowlength; j++){	// Iterate through row
		
			if ((counter + 1 >= elements) && !(process > sizearrlength -1)) 
			{ 
				counter = 0;
				process++;
				destoffset += elements*rowcnt;
				elements = sizearr[process];
			}
			
			sendbuf[i*elements + counter + destoffset] = vector[i*rowlength + j];
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
#ifdef HAVE_MPI
	return MPI_Wtime();
#elif defined(HAVE_OPENMP)
	return omp_get_wtime();
#else
	struct timeval tmpTime;
	gettimeofday(&tmpTime,NULL);
	return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
#endif
}

void freeVector(Vector inpt){
	free(inpt->data);
	// free(inpt);
}

void freeMatrix(Matrix inpt){
	for (int i = 0; i < inpt->rows; ++i){
		free(inpt->row[i]);
	}
	free(inpt->row);
	free(inpt->as_vec);
	free(inpt->data[0]);
	free(inpt->data);
	// free(inpt);
}

void VariableInitz(int n, double *h, int *globRowLen, int *tempMatSz){
	*h = 1.0/n;
	*h *= *h;
	*globRowLen = n-1;
	*tempMatSz = n*4;
}

#define TEST 0
int print = 0;

/*DO NOT USE THE COMMONS LIBRARY!
*IF ANYTHING IN THE COMMONS LIBRARY IS OF USE, COPY IT OVER.
*THE COMMONS LIBRARY HAS TOO MUCH "UNNECESSARY DATA/JUNK"
*/

int main(int argc, char *argv[]){
	Matrix	/*The "work"-matrix*/matrix, \
			/*The transposed version of the matrix*/transpMat;
	Vector	/*The diagonal matrix*/diagMat, \
			/*The fourier-transform temp-store-matrix*/tempMat;

	//Lists for holding how many rows per process (due to MPI division of labour), AKA MPI variables...
	int *size, *displacement, rank = 0, mpiSize = 1, acquired;

	//Global variables
	double h;
	int globRowLen, n;

	//Process specific variables
	int tempMatSz, locMatSz, procRowAmnt;

	/*		General MPI startup/setup		*/
	#ifdef HAVE_OPENMP
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &acquired);
	#else
	MPI_Init(&argc, &argv);
	#endif
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_dup(MPI_COMM_WORLD, &WorldComm);
	MPI_Comm_dup(MPI_COMM_SELF, &SelfComm);

	//Check for correct commandline argument
	//Remove the section comments in the if-test when we want to test "live"/production/release version.
	if (argc < 2 /*|| atoi(argv[1]) < 4 || atoi(argv[1])%2 != 0*/){
		printf("Need a problem size! And it should be a power of 2 greater than or equal to 4!\n");
		MPI_Finalize();
		return -1;
	}

	/*			Initializing variables		*/
	n = atoi(argv[1]);
	VariableInitz(n, &h, &globRowLen, &tempMatSz);
	splitVector(globRowLen, mpiSize, &size, &displacement);
	procRowAmnt = size[rank];
	locMatSz = globRowLen*procRowAmnt;

	if(rank == TEST && print){
		printf("n: %i\nglobRowLen: %i\nmpiSize: %i\n\n", n, globRowLen, mpiSize);
		printf("Size vector:\n");
		printIntVector(size, mpiSize);
		printf("Displacement vector:\n");
		printIntVector(displacement, mpiSize);
		printf("\nSize of rank %i matrix: %ix%i\n", TEST, n-1, size[TEST]);
	}

	/*			Initializing structures		*/
	tempMat = createVector(tempMatSz);
	diagMat = createVector(globRowLen);
	matrix = createMatrix(procRowAmnt, globRowLen);
	transpMat = createMatrix(procRowAmnt, globRowLen);
	diagMat->data = (double*) malloc(locMatSz*sizeof(double));

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < globRowLen; ++i){
		//Filling the diagonal matrix
		diagMat->data[i] = (double) 2.0*(1.0-cos(i + 1.0)*M_PI/(double)n);
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < locMatSz; ++i){
		//Filling up the work-matrix
		matrix->as_vec->data[i] = h;
	}

	if(rank == TEST && print){
		printf("Diagonal matrix for rank == %i:\n", TEST);
		printDoubleVector(diagMat->data, globRowLen);
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procRowAmnt; ++i){
		//Implementation of the first fst_() call
		tempMat->data = (double*) calloc(tempMatSz, sizeof(double));
		fst_(matrix->data[i], &globRowLen, tempMat->data, &tempMatSz);
	}

	/*		Implementation of the first transpose			*/
	// TODO: Call sendArrange to set up send buffer
	double * sendbuf;	// TODO: Allocate
	MPI_Alltoallv(sendbuf, /*splitvec result*/, /*send displacement array*/, MPI_DOUBLE, /*recvcount array. Splitvec result here too*/, /*recv displ array. same as send one*/, MPI_DOUBLE, /*communicator*/);
	
	// TODO: Invert matrix
	
	//ERLEND! =DDD

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procRowAmnt; ++i){
		//Implementation of the first fstinv_() call
		tempMat->data = (double*) calloc(tempMatSz, sizeof(double));
		fstinv_(transpMat->data[i], &globRowLen, tempMat->data, &tempMatSz);
	}

	/*		Implementation of the "tensor" operation		*/
	//Which for-loop level should get the open mp pragma?
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procRowAmnt; ++i){
		for (int j = 0; j < globRowLen; ++j){
			transpMat->data[i][j] /= diagMat->data[j] + diagMat->data[i];
		}
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procRowAmnt; ++i){
		//Implementation of the second fst_() call
		tempMat->data = (double*) calloc(tempMatSz, sizeof(double));
		fst_(transpMat->data[i], &globRowLen, tempMat->data, &tempMatSz);
	}

	/*		Implementation of the second transpose			*/
	
	// TODO: Same as prev transpose

	//ERLEND! =DDD

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procRowAmnt; ++i){
		//Implementation of the second fstinv_() call
		tempMat->data = (double*) calloc(tempMatSz, sizeof(double));
		fstinv_(matrix->data[i], &globRowLen, tempMat->data, &tempMatSz);
	}

	/*		Print time? (not yet implemented)				*/

	/*		Closing up and freeing variables				*/
	freeMatrix(matrix);
	freeMatrix(transpMat);
	freeVector(tempMat);
	freeVector(diagMat);

	if(rank == 0){
		MPI_Comm_free(&WorldComm);
		MPI_Comm_free(&SelfComm);
	}

	MPI_Finalize();
	return 0;
}
