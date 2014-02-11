#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h> //Already included in common.h (below)

#include "../../examples/common/common.h"
#include "ex4.h"

void fillVectorNumerically(Vector inpt){
	for (int i = 0; i < inpt->len; ++i){
		inpt->data[i] = 1.0/pow((i+1), 2.0);
	}
}

double getVectorSum(Vector inpt, int start, int end){
	double sum = 0.0;

	//Let the compiler decide whether OpenMP is run or not.
	#pragma omp parallel for schedule(dynamic, 5) reduction(+:sum)
	for (int i = end; i >= start; --i){
		sum += inpt->data[i];
	}

	return sum;
}

// TODO: TOP PRIORITY! Use common.h convenience functions for MPI
//First attempt/iteration of above TODO has been done.
// TODO: Add at argv[2] for number of MPI ranks
int main(int argc, char *argv[]){
	//Initialization of variables
	int vecLength = 0, rank = 0, size = 0, useOpenMP = 0;
	Vector kValues, difference, *vectorList;

	//Initialization of time
	double wTime = WallTime();

	//Initialization of MPI
	init_app(&argc, argv, &rank, &size); //We can safely call it at this point in time because only rank 0 should be initialized afaik?

	//Pre-work to be done by the "master node", AKA node with rank==0:
	if(rank == 0){
		if(argc <= 2){ //Check that we got vecLength and useOpenMP params
			//The above check _should_(?) work if I have understood how MPI_Init() does its job. My impression is that it should receive argc and argv as call-by-reference, and "clean them up", such that they correspond with the expected arguments given at runtime.
			//(AFAIK they do not at normally correspond as such due to MPI complexities).
			printf("Too few arguments given!\n\tProgram aborted.\n");
			return -1;
		} else{
			vecLength = (int) pow(2.0, (double) atoi(argv[1]));
			useOpenMP = atoi(argv[2]);

			// Let rank 0 generate vector "v"
			#ifdef HAVE_MPI //Create and fill data structures
				Vector numericV = createVectorMPI(vecLength, &WorldComm/*Is this correct? Should it be SelfComm?*/, 0 /*No idea what this variable is for...*/);
				fillVectorNumerically(numericV);
			#else
				//Set up vectors and "help-vectors" for computing the difference with different k-values
				kValues = createVector(12);
				difference = createVector(12);
				vectorList = (Vector*) malloc(12*sizeof(Vector));
				for (int i = 0; i < kValues->glob_len; ++i){
					kValues->data[i] = 3+i;
					vectorList[i] = createVector(3+i);
					fillVectorNumerically(vectorList[i]);
				}
			#endif
			//TODO: Ensure that the above code works as intended/it should, and then split (scatter?) the work to each node
		}
	}

	//TODO: Make each node execute its job. No clue how... as of yet...
	#ifdef HAVE_MPI
		// TODO: Complete function call, Scatter data to MPI ranks
		/* void *sendbuf, int sendcnt, MPI_Datatype sendtype,
				void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
				MPI_Comm com */
		/*int scatter_res = MPI_Scatter(void *sendbuf, int sendcnt, vector, void *recvbuf, int recvcnt, vector, 0, MPI_COMM_WORLD);*/

	// TODO: Convert to summing on local vector-piece if MPI is in use.
	//Compute sum of "v" on processor(s).
	//double vSum = getVectorSum(numericV);

	// TODO: Complete function call, Gather sums to rank 0 (preferably by binary tree for efficiency)
	/* void *sendbuf, int sendcnt, MPI_Datatype sendtype,
			void *recvbuf, int recvcnt, MPI_Datatype recvtype,
			int root, MPI_Comm comm */
	/*int gather_res = MPI_Gather(void *sendbuf, int sendcnt, vector,
			void *recvbuf, int recvcnt, vector,
			0, MPI_COMM_WORLD);*/
	#else

	#endif

	if(rank == 0){ //If no MPI, just execute the code as usual.
		//Compute the differences
		for (int i = 0; i < kValues->glob_len; ++i){
			difference->data[i] = (pow(PI, 2.0)/6.0) - getVectorSum(vectorList[i]);
		}
	}

	//TODO: Collect and add up all sums from all nodes into the node with rank == 0 again, so that it can report the result. All MPI implementation of this is lacking.
	#ifdef HAVE_MPI

	#else
	//Perhaps unnecessary clause? Maybe the "else" can be dropped. I suspect so...
	#endif

	if(rank == 0){ //Enclosure of all "non-MPI" calls and variable-use with if(rank == 0), so as to make sure to avoid wrong use of variables during runtime.
		//Get the elapsed time
		//TODO: Confirm this is the correct way
		wTime = WallTime() - wTime;
		//TODO: Find out how to report elapsed time in human-readable format.

		//Print the differences
		printf("\nBelow are the differences of the sums of vectors with the values: v[i] = 1/i^2,\nwhere n = 2^k, and k has been given different values:\n");
		for (int i = 0; i < kValues->glob_len; ++i){
			printf("With k = %d, the difference is: %.2f\n",
				(int) kValues->data[i], difference->data[i]);
		}
	}

	//MPI cleanup
	close_app();

	return 0;
}
