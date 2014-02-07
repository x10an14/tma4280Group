#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "../../examples/common/common.h"
#include "ex4.h"

int rank,					// Rank of this process
	size,					// Total number of processes
	vecLength;				// Vector length

//MPI communicator
MPI_comm comm;

//MPI datatype
MPI_datatype vector;

//Function for creating and committing MPI datatypes
void create_types(){
	//Creating vector type
	MPI_Type_vector(vecLength, 1, vecLength, MPI_FLOAT, &vector);
	MPI_Type_commit(&vector);
}

void fillVectorNumerically(Vector inpt){
	for (int i = 0; i < inpt->glob_len; ++i){
		inpt->data[i] = 1.0/pow((i+1), 2.0);
	}
}

double getVectorSum(Vector inpt){
	double sum = 0.0;

	#pragma omp parallel for schedule(dynamic, 5) reduction(+:sum) //Add if-defs for later
	for (int i = (inpt->glob_len - 1); i >= 0; --i){
		sum += inpt->data[i];
	}

	return sum;
}

int main(int argc, char const *argv[]){
	vecLength = 0;

	if(argc <= 1){
		printf("Too few arguments given!\n\tProgram aborted.\n");
		return -1;
	} else{
		vecLength = atoi(argv[1]);
	}

	//Initialize MPI, get rank and size
	MPI_Init(&argc, &argv); //argc: number of args, argv: arg-vector //TODO: Change?
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//TODO: We don't need a communicator for 1D topology?

	// Create data-type(s)
	create_types();

	//TODO: Only do if rank 0
	//Generate vector "v"
	Vector numericV = createVector(vecLength);
	fillVectorNumerically(numericV);

	//TODO: Scatter data to MPI ranks

	//Compute sum of "v" on processor(s).
	//double vSum = getVectorSum(numericV);

	// TODO: Gather sums to rank 0 (preferably by binary tree for efficiency)

	//Set up vectors and "help-vectors" for computing the difference with different k-values
	Vector kValues = createVector(12);
	Vector difference = createVector(12);
	Vector *vectorList = (Vector*) malloc(12*sizeof(Vector));
	for (int i = 0; i < kValues->glob_len; ++i){
		kValues->data[i] = 3+i;
		vectorList[i] = createVector(3+i);
		fillVectorNumerically(vectorList[i]);
	}

	//Compute the differences
	for (int i = 0; i < kValues->glob_len; ++i){
		difference->data[i] = (pow(PI, 2.0)/6.0) - getVectorSum(vectorList[i]);
	}

	//Print the differences
	printf("\nBelow are the differences of the sums of vectors with the values: v[i] = 1/i^2,\nwhere n = 2^k, and k has been given different values:\n");
	for (int i = 0; i < kValues->glob_len; ++i){
		printf("With k = %d, the difference is: %.2f\n",
			(int) kValues->data[i], difference->data[i]);
	}

	//MPI cleanup
	MPI_Type_free(&vector);
	MPI_Finalize();

	return 0;
}
