#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
		printf("[%f", ptr[i][0]);
		for (int j = 1; j < cols; ++j){
			printf(",\t%f", ptr[i][j]);
		}
		printf("]\n");
	}
}

void printDoubleVector(double *ptr, int length){
	printf("[%f", ptr[0]);
	for (int i = 1; i < length; ++i){
		printf(",\t%f", ptr[i]);
	}
	printf("]\n");
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
	//The process will always run on a separate processor than its threads, so it will always keep running while the threads do their work.
	return MPI_Wtime();
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
			if ((counter + 1 > elements) && !(process > sizearrlength -1)){
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

/* Arranger the receivebuffer back into column order because it will be in partial rows
 * afte receiving (see the pdf Parallelization of a fast Poisson solver)
 */
void recvArrange(double *recvbuf, double *outbuf, int collength, int colcnt){
	for(int i = 0; i < collength; i++){
		for(int j = 0; j < colcnt; j++){
			outbuf[j*collength + i] = recvbuf[i*colcnt+j];
		}
	}
}

void packTransp(Matrix inpt, Matrix outpt, int *scount, int *sdisp, int mpiSize){
	int cntr = 0, cols = inpt->cols;
	for (int p = 0; p < mpiSize; ++p){ //For each process
		int j_start = sdisp[p]/cols;
		int j_end = (sdisp[p] + scount[p])/cols;
		for (int i = 0; i < cols; ++i){ //For each coloumn
			for (int j = j_start; j < j_end; ++j){ //Copy the section of the coloumn corresponding to each process into outpt
				outpt->as_vec->data[cntr] = inpt->data[i][j];
				cntr++;
			}
		}
	}
}

void fillWithNaturalNumbers(Matrix inpt, int rank, int *displ, int totSize){
	double val = (displ[rank]*inpt->rows) + 1;
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < totSize; ++i){
		//Filling up the work-matrix
		inpt->as_vec->data[i] = (double) val+i;
	}
}

void fillWithConst(Matrix inpt, double constant){
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < inpt->cols; ++i){
		for (int j = 0; j < inpt->rows; ++j){
			//Filling up the work-matrix
			inpt->data[i][j] = constant;
		}

	}
}

void fillWithAppB(Matrix inpt, int rank, int *displ, double h){
	double offset, a = (5*M_PI*M_PI) * h;
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < inpt->cols; ++i){
		offset = (double) displ[rank];
		for (int j = 0; j < inpt->rows; ++j){
			inpt->data[i][j] = a * sin(M_PI*(i + offset + 1)) * sin(2*M_PI*(j + 1));
		}
	}
}

void callFourier(Matrix inpt, Matrix tmp){
	#ifdef HAVE_OPENMP
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < inpt->cols; ++i){
			fst_(inpt->data[i], &inpt->rows, tmp->data[omp_get_num_threads()], &tmp->rows);
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
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < inpt->cols; ++i){
			fstinv_(inpt->data[i], &inpt->rows, tmp->data[omp_get_num_threads()], &tmp->rows);
		}
	#else
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < inpt->cols; ++i){
			fstinv_(inpt->data[i], &inpt->rows, tmp->data[0], &tmp->rows);
		}
	#endif
}

void unpackTransp(Matrix outpt, Matrix inpt){
	int cntr = 0;
	for (int i = 0; i < outpt->rows; ++i){
		for (int j = 0; j < outpt->cols; ++j){
			outpt->data[j][i] = inpt->as_vec->data[cntr];
			++cntr;
		}
	}
}

double exactSolAppB(int col, int row, double h){
	//Return the exact solution to the x at coordinate col and row.
	return sin(M_PI*(col+1)*h)*sin(2*M_PI*(row+1)*h);
}

double linearAverage(Matrix inpt, double (*funcp)(int, int, double), double h){
	double avgErr = 0.0, col_err; int rows = inpt->rows, cols = inpt->cols;
	#pragma omp parallel for schedule(static) private(col_err) reduction(+:avgErr)
	for (int i = 0; i < cols; ++i){
		col_err = 0.0;
		for (int j = 0; j < rows; ++j){
			col_err += fabs(inpt->data[i][j] - (*funcp)(i, j, h));
		}
		avgErr += col_err/((double) rows);
	}
	avgErr /= ((double) cols);
	return avgErr;
}

#define TEST 0
#define PRINT 0
#define MPI_WTIME_IS_GLOBAL 1

int main(int argc, char *argv[]){
	Matrix	/*The "work"-matrix*/matrix, \
			/*The transposed version of the matrix*/transpMat, \
			/*The fourier-transform temp-store-matrix*/tempMat;
	Vector	/*The diagonal matrix*/diagMat;

	int /*Array holding amount of coloumns per process*/*size, \
		/*Array holding the offset so each process can know which coloumns reside in which process
		(when used with the above array)*/*displ, \
		/*Array telling each process how much data to send to all the others
		(potentially different values for each process in this array)*/ *scount, \
		/*Displacement array for the above array*/ *sdisp, \
		/*Typical MPI variables.*/ rank = 0, mpiSize = 1, acquired;

	double /*Global variables*/h, initTime = WallTime(), runTime, declareTime, singleFourier, \
		preTranspTime, singleTransp, halfWork, singleTensor, totRun, totTime;
	int /*Global variables*/globColLen, n, \
		/*Process specific variables*/tempMatSz, locMatSz, procColAmnt;

	/*		General MPI startup/setup		*/
	#ifdef HAVE_OPENMP
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &acquired);
	#else
	MPI_Init(&argc, &argv);
	#endif
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Check for correct commandline argument
	if (argc < 2 || atoi(argv[1]) < 4 || atoi(argv[1]) & (atoi(argv[1])-1) != 0){
		printf("Need a problem size! And it should be a power of 2 greater than or equal to 4!\n");
		MPI_Finalize();
		return -1;
	}

	declareTime = WallTime();

	/*			Initializing variables		*/
	n = atoi(argv[1]); h = (double) 1.0/(n*n); globColLen = n-1; tempMatSz = n*4;
	splitVector(globColLen, mpiSize, rank, &size, &displ, &scount, &sdisp);
	procColAmnt = size[rank]; locMatSz = globColLen*procColAmnt;

	if(rank == TEST && PRINT){
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
	/*				Erlends implementation					*/
	//double *sendbuf = calloc(locMatSz, sizeof(double));

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < diagMat->len; ++i){
		//Filling the diagonal matrix
		diagMat->data[i] = (double) 2.0*(1.0-cos(i + 1.0)*M_PI/(double)n);
	}

	//fillWithNaturalNumbers(matrix, rank, displ, locMatSz);
	//fillWithConst(matrix, locMatSz, h);
	//fillWithFunction() //For part e) of p6
	fillWithAppB(matrix, rank, displ, h);

	if(rank == TEST && PRINT){
		printf("\nDiagonal matrix for rank == %i:\n", TEST);
		printDoubleVector(diagMat->data, globColLen);
	}

	runTime = WallTime();
	callFourier(matrix, tempMat);
	singleFourier = WallTime() - runTime;
	/*		Implementation of the first transpose			*/
	if (rank == TEST && PRINT){
		printf("\nBefore transpose:\n");
		printf("matrix:\n");
		printDoubleMatrix(matrix->data, matrix->cols, matrix->rows);
		printf("transpMat:\n");
		printDoubleVector(transpMat->data[0], locMatSz);
	}

	preTranspTime = WallTime();
	/*				Christians implementation				*/
	packTransp(matrix, transpMat, scount, sdisp, mpiSize);
	MPI_Alltoallv(transpMat->data[0], scount, sdisp, MPI_DOUBLE, matrix->data[0], scount, sdisp, MPI_DOUBLE, MPI_COMM_WORLD);
	unpackTransp(transpMat, matrix);

	singleTransp = WallTime() - preTranspTime;

	/*				Erlends implementation					*/
	//sendArrange(sendbuf, matrix->data[0], globColLen, procColAmnt, size, mpiSize, displ);
	//double *recvbuf = malloc(sizeof(double)*globColLen*procColAmnt); //HAR FLYTTET DENNE LINJEN LENGERE OPP! NÅ BLIR DEN INITIALISERT DOBBELT!
	//MPI_Alltoallv(sendbuf, scount, sdisp, MPI_DOUBLE, recvbuf, scount, sdisp, MPI_DOUBLE, MPI_COMM_WORLD);
	//recvArrange(recvbuf, matrix->data[0], globColLen, procColAmnt);

	if(rank == TEST && PRINT){
		printf("\nAfter transpose:\n");
		printf("transpMat:\n");
		printDoubleMatrix(transpMat->data, transpMat->cols, transpMat->rows);
		printf("matrix:\n");
		printDoubleVector(matrix->data[0], locMatSz);
		printf("\n");
	}

	callFourierInvrs(transpMat, tempMat);
	halfWork = WallTime();
	/*		Implementation of the "tensor" operation		*/
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < procColAmnt; ++i){
		for (int j = 0; j < globColLen; ++j){
			transpMat->data[i][j] /= diagMat->data[j] + diagMat->data[i];
		}
	}
	singleTensor = WallTime() - halfWork;
	halfWork -= runTime;

	callFourier(transpMat, tempMat);
	/*		Implementation of the second transpose			*/
	/*				Christians implementation				*/
	packTransp(transpMat, matrix, scount, sdisp, mpiSize);
	MPI_Alltoallv(matrix->data[0], scount, sdisp, MPI_DOUBLE, transpMat->data[0], scount, sdisp, MPI_DOUBLE, MPI_COMM_WORLD);
	unpackTransp(matrix, transpMat);

	/*				Erlends implementation					*/
	//Arrange send buffer(using same buffer as last runTime)
	/*sendArrange(sendbuf, matrix->data[0], globColLen,procColAmnt, size, mpiSize);
	MPI_Alltoallv(&sendbuf, size, displ, MPI_DOUBLE, matrix->data[0], size, displ, MPI_DOUBLE, MPI_COMM_WORLD);*/

	callFourierInvrs(matrix, tempMat);
	/*				Print timings							*/
	totRun = WallTime() - runTime;
	totTime = WallTime() - initTime;
	if(rank == TEST){
		double decTime = (declareTime - initTime)*1000;
		/*printf("Initializatin + MPI_Init: %fms.\n", decTime);
		printf("Declaration and filling out variables: %fms.\n", ((runTime - initTime)*1000) - decTime);
		printf("Single first Fourier-call: %fms.\n", singleFourier*1000);
		printf("Single first transposition: %fms.\n", singleTransp*1000);
		printf("Tensor operation: %fms.\n", singleTensor*1000);
		printf("Half-way time: %fms.\n", halfWork*1000);
		printf("All tensor, transp, and fourier calls: %fms.\n", totRun*1000);*/

		//The time not including lines 278-312 in main() (from start of main, until declareTime = WallTime();)
		printf("time: %fms\n", (totTime - (declareTime - initTime))*1000);
	}

	/*					Error checking						*/
	double (*fp)(int, int, double), procAvgErr, globAvgErr;
	fp = exactSolAppB; procAvgErr = linearAverage(matrix, fp, h);
	MPI_Reduce(&procAvgErr, &globAvgErr, 1, MPI_DOUBLE, MPI_SUM, TEST, MPI_COMM_WORLD);
	if(rank == TEST){
		globAvgErr /= mpiSize;
		printf("error: %g\n\n", globAvgErr);
	}

	/*		Closing up and freeing variables				*/
	freeMatrix(matrix); freeMatrix(tempMat); freeMatrix(transpMat);
	freeVector(diagMat); free(size); free(displ); free(scount); free(sdisp);
	MPI_Finalize();

	return 0;
}
