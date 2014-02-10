#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "../../examples/common/common.h"
#include "ex4.h"

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
	int vecLength = 0, rank = 0, size = 0;
	double wTime = 0.0;

#ifdef HAVE_MPI
	//We can safely call it at this point in time because only rank 0 should be initialized afaik?
	init_app(argc, argv, &rank, &size);
	wTime = WallTime();
#endif

	if(rank == 0){
		if(argc <= 1){
			printf("Too few arguments given!\n\tProgram aborted.\n");
			return -1;
		} else{
			vecLength = (int) pow(2.0, (double) atoi(argv[1]));
			// Let rank 0 generate vector "v"
#ifdef HAVE_MPI
			Vector numericV = createVectorMPI(vecLength, &WorldComm/*Is this correct? Should it be SelfComm?*/, 0 /*No idea what this variable is for...*/);
			(numericV);
#else
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
#endif
		}
	}


	//MPI cleanup
	close_app();

	return 0;
}
