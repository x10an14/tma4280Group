#include "ex6.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

/* THIS FUNCTION DOES NOT ALLOCATE ACTUAL VECTOR DATA MEMORY! */
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

		if (globLen % size && i >= (size - (globLen % size))){
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
void sendArrange(double *sendbuf, double *vector, double *rows, int rowlength, int rowcnt, int *sizearr, int sizearrlength)
{
	int elements, counter, process;
	
	for(int i = 0; i < rowcnt; i++){	// One pass per row
		
		process = 0; 
		elements = sizearr[process];	// Elements to go to the first process
		counter  = 0;					// Counter for how many elements have been moved
		
		for(int j = 0; j < rowlenght; j++){	// Iterate through row
			if ((counter + 1 > elements) && !(process > sizearrlength -1)) 
			{ 
				counter = 0;
				process++;
				elements = sizearr[process];
			}
			
			sendbuf[1*i + j] = vector[];
		}				
	}
}

/* Function to re-arrange the receive-buffer into a correctly continous piece of memory (vector) 
 * 
 */
void recvBufRearr(double *recvbuf, double *vector, double *rows, int rcvprrwprprc, int rowlength, int rowcnt, int processes)
{
	// TODO	
}

#define TEST 0

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
	int matrixSize, matrixTempSize, n = 0, rank = 0, mpiSize = 1, acquired;
	double h = 1.0;

	//General MPI startup/setup
	#ifdef HAVE_OPENMP
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &acquired);
	#else
	MPI_init(&argc, &argv);
	#endif
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_dup(MPI_COMM_WORLD, &WorldComm);
	MPI_Comm_dup(MPI_COMM_SELF, &SelfComm);

	//Check for comandline argument
	if (argc < 2 /*|| atoi(argv[1]) < 4 || atoi(argv[1])%2 != 0*/){
		printf("Need a problem size! And it should be a power of 2 greater than or equal to 4!\n");
		//MPI_Comm_free(&WorldComm);
		//MPI_Comm_free(&SelfComm);
		MPI_Finalize();
		return -1;
	} else{
		n = atoi(argv[1]);
	}

	/*			Initializing variables		*/
	h /= (double) n;
	matrixSize = n-1;
	splitVector(matrixSize, mpiSize, &size, &displacement);

	if(rank == TEST){
		printf("n: %i\nmatrixSize: %i\nmpiSize: %i\n\n", n, matrixSize, mpiSize);
		printf("Size vector:\n");
		printIntVector(size, mpiSize);
		printf("Displacement vector:\n");
		printIntVector(displacement, mpiSize);
	}

	matrixSize = size[rank];
	matrixTempSize = matrixSize*4;

	/*			Initializing structures		*/
	//Use a MPI communicator? We have two defined in the .h file. Maybe make createMatrix() set that for us.
	diagMat = createVector(matrixSize);
	f_TempMat = createVector(matrixTempSize);
	matrix = createMatrix(matrixSize, matrixSize);
	transpMat = createMatrix(matrixSize, matrixSize);
	diagMat->data = (double*) malloc(matrixSize*sizeof(double));
	f_TempMat->data = (double*) malloc(matrixTempSize*sizeof(double));


	if(rank == TEST){

	}


	/*	The rest of the initialization of structures, making sure that the matrix and diagMat variables only contain what each process needs	*/

	/*		Implementation of the fst_() call			*/

	/*		Implementation of the transpose				*/

	/*		Implementation of the fstinv_() call		*/

	/*		Implementation of the "tensor" operation		*/

	/*	//Namely this:
	*	transpMat->data[i][j] /= (diagMat->data[i] + diagMat->data[j]);
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
