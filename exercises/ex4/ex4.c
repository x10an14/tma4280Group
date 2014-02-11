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
	int vecLength = 0, rank = 0, size = 1, k = 1;
	double glob_sum  = 0.0, loc_sum = 0.0, actual_sum = (pow(PI, 2.0)/6.0);
	Vector kVector;

	//Initialization of time
	double wTime = WallTime();

	//Initialization of MPI
	init_app(&argc, argv, &rank, &size);

	if(argc <= 1){ //Check that we got vecLength and useOpenMP params
		//The above check _should_(?) work if I have understood how MPI_Init() does its job. My impression is that it should receive argc and argv as call-by-reference, and "clean them up", such that they correspond with the expected arguments given at runtime.
		//(AFAIK they do not at normally correspond as such due to MPI complexities).
		printf("Too few arguments given!\n\tProgram aborted.\n");
		return -1;
	}
	k = atoi(argv[1]);
	vecLength = (int) pow(2.0, (double) k);

	//Pre-work to be done by the "master node", AKA node with rank==0:
	if(rank == 0){
		//Let rank 0 generate vector
		kVector = createVector(vecLength);
		fillVectorNumerically(kVector);

		#ifdef HAVE_MPI
			//Send vector pieces to all nodes (ranks)
			/* void *sendbuf, int sendcnt, MPI_Datatype sendtype,
					void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
					MPI_Comm com */
			/*int scatter_res = MPI_Scatter(void *sendbuf, int sendcnt, vector, void *recvbuf, int recvcnt, vector, 0, MPI_COMM_WORLD);*/
		#else
			//Progress as if there's no extra nodes.

		#endif
	}

	//Figure out the smart way to have a buffer locally and universally, so that each process (rank) can compute it's local sum
	if(size == 0){
		glob_sum = getVectorSum(kVector, 0, kVector->glob_len);
	} else{
		//For making the code more readable, I make each rank calculate its steps.
		int step = pow(2.0, 1.0*rank*k/size), next_step = (pow(2.0, 1.0*(rank+1)*k/size) - 1);

		//Probably won't work, we need to malloc each local sum?
		if(rank == 0){
			loc_sum = getVectorSum(kVector, 0, (int) pow(2.0, 1.0*k/size)-1.0);
		} else{
			loc_sum = getVectorSum(kVector, step, next_step);
		}
	}

	//TODO: Collect and add up all sums from all nodes into the node with rank == 0 again, so that it can report the result. All MPI implementation of this is lacking.
	#ifdef HAVE_MPI //Receive the data from all MPI processes
		//TODO: Figure out how MPI_Receive works.
		//TODO: This section of code should collect all the loc_sums into glob_sum
		/*int gather_res = MPI_Gather(void *sendbuf, int sendcnt, vector,
			void *recvbuf, int recvcnt, vector,
			0, MPI_COMM_WORLD);*/
	#endif

	if(rank == 0){ //Enclosure of all "non-MPI" calls and variable-use with if(rank == 0), so as to make sure to avoid wrong use of variables during runtime.
		//Get the elapsed time
		//TODO: Confirm this is the correct way
		wTime = WallTime() - wTime;
		//TODO: Find out how to report elapsed time in human-readable format.

		//If we don't have MPI, the sum should be in the above buffer.
		glob_sum = actual_sum - glob_sum;

		//Print the differences
		printf("\nBelow is the difference of the sum between the vector of length n, with n of vectors with the values: v[i] = 1/i^2,\nwhere n = 2^k, and k has been given different values:\n");
		printf("With k = %d, the difference is: %.2f\n", k, glob_sum);
		printf("Time it took to complete: %f seconds.\n", wTime);
	}

	//MPI cleanup
	close_app();

	return 0;
}
