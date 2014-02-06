#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../examples/common/common.h"
#include "ex4.h"

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
	int vecLength = 0;

	if(argc <= 1){
		printf("Too few arguments given!\n\tProgram aborted.\n");
		return -1;
	} else{
		vecLength = atoi(argv[1]);
	}

	//Generate vector "v"
	Vector numericV = createVector(vecLength);
	fillVectorNumerically(numericV);

	//Compute sum of "v" on one processor.
	//double vSum = getVectorSum(numericV);

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

	return 0;
}