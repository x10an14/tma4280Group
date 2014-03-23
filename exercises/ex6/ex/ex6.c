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

void printDoubleMatrix(double **ptr, int rows, int cols){
	for (int i = 0; i < rows; ++i){
		printf("[%.2f", ptr[i][0]);
		for (int j = 1; j < cols; ++j){
			printf(",\t%.2f", ptr[i][j]);
		}
		printf("]\n");
	}
}

void printDoubleVector(double *ptr, int length){
	printf("[%.2f", ptr[0]);
	for (int i = 1; i < length; ++i){
		printf(",\t%.2f", ptr[i]);
	}
	printf("]\n");
}

void printIntMatrix(int **ptr, int rows, int cols){
	for (int i = 0; i < rows; ++i){
		printf("[%i", ptr[i][0]);
		for (int j = 1; j < cols; ++j){
			printf(",\t%i", ptr[i][j]);
		}
		printf("]\n");
	}
}

void printIntVector(int *ptr, int length){
	printf("[%i", ptr[0]);
	for (int i = 1; i < length; ++i){
		printf(",\t%i", ptr[i]);
	}
	printf("]\n");
}

//Copied from common.c, minor edits.
void splitVector(int globLen, int size, int rank, int** len, int** displ, int **scount, int **sdisp){
	*len = (int*) malloc(size*sizeof(int));
	*displ = (int*) calloc(size, sizeof(int));
	*sdisp = (int*) malloc(size*sizeof(int));
	*scount = (int*) malloc(size*sizeof(int));

	//"Old" splitvector loop
	for (int i=0; i < size; ++i){
		(*len)[i] = globLen/size;

		if (globLen % size && i >= (size - globLen % size)){
			(*len)[i]++;
		}

		if (i < size-1){
			(*displ)[i+1] = (*displ)[i]+(*len)[i];
		}
	}

	//Added loop for transpose purposes
	for (int i = 0; i < size; ++i){
		(*scount)[i] = (*len)[rank]*(*len)[i];
		int tmp = 0;
		for (int j = 1; j <= i; ++j){
			tmp += (*len)[rank]*(*len)[j-1];
		}
		(*sdisp)[i] = tmp;
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

/* Arranges the sendbuffer properly before sending
 * Assumes the columns are arranged continually in the vector
 */
void sendArrange(double *sendbuf, double *vector, int collength, int colcnt, int *sizearr, int sizearrlength, int * displvec){
	int elements, counter, process, destoffset;

	for(int i = 0; i < colcnt; i++){	// One pass per col
		process = 0;
		elements = sizearr[process];	// Elements to go to the first process
		counter  = 0;					// Counter for how many elements have been moved
		destoffset = 0;					// Destination offset for elements (i.e. which process is this being sent to)

		for(int j = 0; j < collength; j++){	// Iterate through row
			if ((counter + 1 > elements) && !(process > sizearrlength -1))
			{
				counter = 0;
				process++;
				destoffset = displvec[process];
				elements = sizearr[process];
			}

			sendbuf[i*elements + counter + destoffset] = vector[i*collength + j];
			counter++;
		}
	}
}

void preTransp(Matrix inpt, Matrix outpt, int *size, int *scount, int *sdisp, int mpiSize, int rank){
	int cntr = 0, cols = size[rank];
	for (int p = 0; p < mpiSize; ++p){ //For each process
		/*if (rank == 0){
			printf("\nprocess:%i\n", p);
		}*/
		int j_start = sdisp[p]/cols;
		int j_end = (sdisp[p] + scount[p])/cols;
		/*if (rank == 0){
			printf("j_start:%i\nj_end:%i\n", j_start, j_end);
		}*/
		for (int i = 0; i < cols; ++i){ //For each coloumn
			/*if(rank == 0){
				printf("i:%i\n", i);
			}*/
			for (int j = j_start; j < j_end; ++j){ //Copy the section of the coloumn corresponding to each process into outpt
				/*if (rank == 0){
					printf("j:%i\n", j);
					printf("Copying value inpt->data[i][j]: %.2f into cntr:%i in outpt.\n", inpt->data[i][j], cntr);
				}*/
				outpt->as_vec->data[cntr] = inpt->data[i][j];
				cntr++;
			}
		}
	}
}

void fillWithNaturalNumbers(Matrix inpt, int rank, int *size, int totSize){
	double val = 1;

	if (rank > 0){
		int start = 0;
		for (int i = 0; i < rank; ++i){
			start += size[i];
		}
		val += start*inpt->rows;
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < totSize; ++i){
		//Filling up the work-matrix
		inpt->as_vec->data[i] = (double) val+i;
	}
}

void fillWithConstant(Matrix inpt, int totSize, double constant){
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < totSize; ++i){
		//Filling up the work-matrix
		inpt->as_vec->data[i] = constant;
	}
}

void callFourier(Matrix inpt, Matrix tmp){
	#ifdef HAVE_OPENMP
		int thrds = tmp->cols, runs = inpt->cols/thrds, \
					remains = inpt->cols%thrds, cntr = 0;
		for (int i = 0; i < runs; ++i){
			#pragma omp parallel for schedule(static) shared(cntr)
			for (int j = 0; j < thrds; ++j){
				fst_(inpt->data[cntr], &inpt->rows, tmp->data[j], &tmp->rows);
				++cntr;
			}
		}

		//If inpt->cols%threads != 0, then the following for-loop is necessary
		#pragma omp parallel for schedule(static) shared(cntr)
		for (int i = 0; i < remains; ++i){
			fst_(inpt->data[cntr], &inpt->rows, tmp->data[i], &tmp->rows);
			++cntr;
		}
	#else
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < inpt->cols; ++i){
			fst_(inpt->data[i], &inpt->rows, tmp->data[0], &tmp->rows);
		}
	#endif
}

void callFourierInvrs(Matrix inpt, Matrix tmp){
	#ifdef HAVE_OPENMP
		int thrds = tmp->cols, runs = inpt->cols/thrds, \
					remains = inpt->cols%thrds, cntr = 0;
		for (int i = 0; i < runs; ++i){
			#pragma omp parallel for schedule(static) shared(cntr)
			for (int j = 0; j < thrds; ++j){
				fstinv_(inpt->data[cntr], &inpt->rows, tmp->data[j], &tmp->rows);
				++cntr;
			}
		}

		//If inpt->cols%threads != 0, then the following for-loop is necessary
		#pragma omp parallel for schedule(static) shared(cntr)
		for (int i = 0; i < remains; ++i){
			fstinv_(inpt->data[cntr], &inpt->rows, tmp->data[i], &tmp->rows);
			++cntr;
		}
	#else
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < inpt->cols; ++i){
			fst_(inpt->data[i], &inpt->rows, tmp->data[0], &tmp->rows);
		}
	#endif
}

#define TEST 1
int print = 1;

int main(int argc, char *argv[]){
	Matrix	/*The "work"-matrix*/matrix, \
			/*The transposed version of the matrix*/transpMat, \
			/*The fourier-transform temp-store-matrix*/tempMat;
	Vector	/*The diagonal matrix*/diagMat;

	//Lists for holding how many cols per process (due to MPI division of labour), AKA MPI variables++...
	int *size, *displ, *scount, *sdisp, rank = 0, mpiSize = 1, acquired;

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
	h = (double) 1.0/(n*n); globColLen = n-1; tempMatSz = n*4;
	splitVector(globColLen, mpiSize, rank, &size, &displ, &scount, &sdisp);
	procColAmnt = size[rank]; locMatSz = globColLen*procColAmnt;

	if(rank == TEST && print){
		printf("\nn: %i\nglobColLen: %i\nmpiSize: %i\n\n", n, globColLen, mpiSize);
		printf("Size vector:\n");
		printIntVector(size, mpiSize);
		printf("Displacement vector:\n");
		printIntVector(displ, mpiSize);
		printf("\nSize of rank %i matrix: %ix%i\n", TEST, procColAmnt, globColLen);
		printf("Scount:\n");
		printIntVector(scount, mpiSize);
		printf("Sdisp:\n");
		printIntVector(sdisp, mpiSize);
	}

	/*			Initializing structures		*/
	diagMat = createVector(globColLen);
	matrix = createMatrix(globColLen, procColAmnt);
	transpMat = createMatrix(globColLen, procColAmnt);
	tempMat = createMatrix(getMaxThreads(), tempMatSz);
	diagMat->data = (double*) malloc(globColLen*sizeof(double));
	//Transpose temporary buffer
	double *sendbuf = calloc(locMatSz, sizeof(double));

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < diagMat->len; ++i){
		//Filling the diagonal matrix
		diagMat->data[i] = (double) 2.0*(1.0-cos(i + 1.0)*M_PI/(double)n);
	}

	fillWithNaturalNumbers(matrix, rank, size, locMatSz);
	//fillWithConstant(matrix, locMatSz, h);

	time = WallTime();

	if(rank == TEST && print){
		printf("\nDiagonal matrix for rank == %i:\n", TEST);
		printDoubleVector(diagMat->data, globColLen);
	}

	/*		Implementation of the first transpose			*/
	if (rank == TEST && print){
		printf("\nBefore transpose:\n");
		printf("matrix:\n");
		printDoubleMatrix(matrix->data, matrix->cols, matrix->rows);
		printf("transpMat:\n");
		printDoubleVector(transpMat->data[0], locMatSz);
	}

				/*Christians implementation*/
	preTransp(matrix, transpMat, size, scount, sdisp, mpiSize, rank);
	MPI_Alltoallv(transpMat->data[0], scount, sdisp, MPI_DOUBLE, matrix->data[0], scount, sdisp, MPI_DOUBLE, WorldComm);

				/*Erlends implementation*/
	//sendArrange(sendbuf, matrix->data[0], globColLen, procColAmnt, size, mpiSize, displ);
	//MPI_Alltoallv(sendbuf, scount, sdisp, MPI_DOUBLE, matrix->data[0], scount, sdisp, MPI_DOUBLE, WorldComm);

	if(rank == TEST && print){
		printf("\nAfter transpose:\n");
		printf("matrix:\n");
		printDoubleMatrix(matrix->data, matrix->cols, matrix->rows);
		printf("transpMat:\n");
		printDoubleVector(transpMat->data[0], locMatSz);
		printf("\n");
	}

	/*callFourierInvrs(matrix, tempMat);

	/*		Implementation of the "tensor" operation		*/
	//Which for-loop level should get the open mp pragma?
	/*#pragma omp parallel for schedule(static)
	for (int i = 0; i < procColAmnt; ++i){
		for (int j = 0; j < globColLen; ++j){
			matrix->data[i][j] /= diagMat->data[j] + diagMat->data[i];
		}
	}

	callFourier(matrix, tempMat);

	/*		Implementation of the second transpose			*/

	// Arrange send buffer(using same buffer as last time
	/*sendArrange(sendbuf, matrix->data[0], globRowLen,procRowAmnt, size, mpiSize, displacement);
	MPI_Alltoallv(&sendbuf, size, displacement, MPI_DOUBLE, matrix->data[0], size, displacement, MPI_DOUBLE, WorldComm);

	//Arrange send buffer(using same buffer as last time)
	/*preTransp(matrix, transpMat, size, scount, sdisp, mpiSize, rank);
	MPI_Alltoallv(transpMat->data[0], scount, sdisp, MPI_DOUBLE, matrix->data[0], scount, sdisp, MPI_DOUBLE, WorldComm);

	sendArrange(sendbuf, matrix->data[0], globColLen,procColAmnt, size, mpiSize);
	MPI_Alltoallv(&sendbuf, size, displ, MPI_DOUBLE, matrix->data[0], size, displ, MPI_DOUBLE, WorldComm);

	callFourierInvrs(matrix, tempMat);

	/*		Print time? (not yet implemented)				*/
	if(rank == TEST){
		time = WallTime() - time;
		printf("t: %g\n", time);
	}

	/*		Closing up and freeing variables				*/
	freeMatrix(matrix);
	freeMatrix(tempMat);
	freeMatrix(transpMat);
	freeVector(diagMat);
	if(rank == 0){
		MPI_Comm_free(&WorldComm);
	}
	MPI_Finalize();

	return 0;
}
