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

int main(int argc, char *argv[]){
	//Initialization of variables
	Vector kVector, fullVector;
	double wTime = WallTime(); //Initialization of time
	int vecLength = 0, rank = 0, size = 1, k = 3, scatter_res = 0, receive_res = 0;
	double glob_sum = 0.0, loc_sum = 0, actual_sum = (pow(PI, 2.0)/6.0);

	//Initialization of MPI
	init_app(&argc, argv, &rank, &size);

	if(argc <= 1){ //Check that we got our k parameter
		printf("Too few arguments given!\n\tProgram aborted.\n");
		return -1;
	}

	//Utilize main()-parameter(s)
	k = atoi(argv[1]);
	vecLength = (int) pow(2.0, (double) k);

	//Pre-work to be done by the "master node", AKA node with rank==0:
	if(rank == 0){
		//Let rank 0 generate vector
		#ifdef HAVE_MPI
			fullVector = createVectorMPI(vecLength, &WorldComm, 1);
		#else
			fullVector = createVector(vecLength);
		#endif
		fillVectorNumerically(fullVector);
	}

	//How we enabled the program to both run with and without MPI
	#ifdef HAVE_MPI
		kVector = createVectorMPI(vecLength/size, &WorldComm, 1); //If we have MPI, make each process makes its own kVector (in addition to fullVector in rank == 0)
		kVector->len = vecLength/size;
		
		//ERROR TESTING:
		int err_rnk1 = 0;
		int err_rnk2 = 0;
		int err_rnk3 = 0;
		
		MPI_Request req[3];
		
		if (rank == 0)
		{
			MPI_Irecv(&err_rnk1, 1, MPI_INT, 1, 0, WorldComm, req);
			MPI_Irecv(&err_rnk2, 1, MPI_INT, 1, 0, WorldComm, req + 1);
			MPI_Irecv(&err_rnk3, 1, MPI_INT, 1, 0, WorldComm,req + 2);
		}
		else
		{
			int err_status = 2;
			// Check pointers
			if (&(kVector->data) == NULL || kVector == NULL || WorldComm == NULL || &(kVector->len) == NULL || &(fullVector->data) == NULL)
			{
				err_status = 1;
			}
			// Send status to rank 0
			MPI_Send(&err_status, 1, MPI_INT, 0, 0, WorldComm);
			
		}
		
		MPI_Status stat[3];
		
		MPI_Wait(req, stat);
		MPI_Wait(req + 1, stat + 1);
		MPI_Wait(req + 2, stat + 2);
		
		if (rank == 0)
		{
			printf("err_rnk1: %d \r\n",err_rnk1);
			printf("err_rnk2: %d \r\n",err_rnk2);
			printf("err_rnk3: %d \r\n",err_rnk3);
			
		}
		// END ERROR testing
		
		//Send vector pieces to all nodes (ranks)
		scatter_res = MPI_Scatter(fullVector->data, kVector->len, MPI_DOUBLE, kVector->data, kVector->len, MPI_DOUBLE, 0, WorldComm);
		if (rank == 0)
			printf("Scatter result code: %d \r\n", scatter_res);
		//Calculate sum in each node/rank
		loc_sum = getVectorSum(kVector, 0, vecLength/size);
		//Receive all sums into the glob_sum variable on rank == 0 node
		receive_res = MPI_Reduce(&loc_sum, &glob_sum, 1, MPI_DOUBLE, MPI_SUM, 0, WorldComm);
		//kVector is of no more use to this program
		freeVector(kVector);
		if (rank == 0)
			freeVector(fullVector);
	#else
		glob_sum = getVectorSum(fullVector, 0, vecLength);
		freeVector(fullVector);
	#endif

	//Finish program
	if(rank == 0){
		//Get the elapsed time
		wTime = WallTime() - wTime;

		//Error reporting
		if(scatter_res != MPI_SUCCESS || receive_res != MPI_SUCCESS){
			printf("scatter_res: %d\nreceive_res: %d\n", scatter_res, receive_res);
		}

		//Print the time and difference
		printf("diff:\t%f\n", (actual_sum - glob_sum)*1000.0);	//The difference with the k given as parameter, multiplied by 1000.
		printf("ms:\t%f\n\n", wTime*1000.0);					//Milliseconds.
	}

	//MPI cleanup
	close_app();

	return 0;
}
