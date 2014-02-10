#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "../../examples/common/common.h"
#include "ex4.h"

int rank,					// Rank of this process
	size,					// Total number of processes
	vecLength;				// Vector vecLength

//MPI datatype
MPI_Datatype vector;

//Function for creating and committing MPI datatypes
void create_types(){
	//Creating vector type
	MPI_Type_vector(vecLength, 1, vecLength, MPI_FLOAT, &vector);
	MPI_Type_commit(&vector);
}

void fillVectorNumerically(Vector inpt){
	for (int i = 0; i < inpt->len; ++i){
		inpt->data[i] = 1.0/pow((i+1), 2.0);
	}
}

double getVectorSum(Vector inpt){
	double sum = 0.0;

	#pragma omp parallel for schedule(dynamic, 5) reduction(+:sum) //Add if-defs for later
	for (int i = (inpt->len - 1); i >= 0; --i){
		sum += inpt->data[i];
	}

	return sum;
}

// TODO: TOP PRIORITY! Use common.h convenience functions for MPI
// TODO: Add at argv[2] for number of MPI ranks
int main(int argc, char *argv[]){
	vecLength = 0;

	if(argc <= 1){
		printf("Too few arguments given!\n\tProgram aborted.\n");
		return -1;
	} else{
		vecLength = atoi(argv[1]);
	}

	// TODO: IMPORTANT! Make all MPI code dependent on presence of MPI. E.g. make IFDEFs
	// Initialize MPI, get rank and size
	MPI_Init(&argc, &argv); // argc: number of args, argv: arg-vector // TODO: Change?

	MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create data-type(s)
    create_types();

	// Let rank 0 generate vector "v"
	if (rank == 0)
	{
		Vector numericV = createVector(vecLength);
		fillVectorNumerically(numericV);
	}

	// TODO: Complete function call, Scatter data to MPI ranks
	/* void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, 
               MPI_Comm com */
	int scatter_res = MPI_Scatter(void *sendbuf, int sendcnt, vector,
               void *recvbuf, int recvcnt, vector, 0, 
               MPI_COMM_WORLD);
	
	// TODO: Convert to summing on local vector-piece if MPI is in use.
	//Compute sum of "v" on processor(s).
	//double vSum = getVectorSum(numericV);
	
	// TODO: Complete function call, Gather sums to rank 0 (preferably by binary tree for efficiency)
	/* void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm */
	int gather_res = MPI_Gather(void *sendbuf, int sendcnt, vector, 
               void *recvbuf, int recvcnt, vector, 
               0, MPI_COMM_WORLD);

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