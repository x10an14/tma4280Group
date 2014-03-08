#include "ex6.h"

Vector createVector(int len){
	Vector result = (Vector) calloc(1, sizof(vector_t));
	result->len = len;
	result->data = NULL;
}

Matrix createMatrix(int rows, int cols){
	Matrix result = (Matrix) calloc(1, sizof(matrix_t));

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

	return result:
}

//Copied from common.c, minor edits.
void splitVector(int globLen, int size, int** len, int** displ){
	*len = calloc(size,sizeof(int));
	*displ = calloc(size,sizeof(int));

	for (i=0; i < size; ++i){
		(*len)[i] = globLen/size;

		if (globLen % size && i >= (size - (globLen % size))){
			(*len)[i]++;
		}

		if (i < size-1){
			(*displ)[i+1] = (*displ)[i]+(*len)[i];
		}
	}
}

/*DO NOT USE THE COMMONS LIBRARY!
*IF ANYTHING IN THE COMMONS LIBRARY IS OF USE, COPY IT OVER.
*THE COMMONS LIBRARY HAS TOO MUCH "UNNECESSARY DATA/JUNK"
*/

int int main(int argc, char *argv[]){
	Matrix matrix, tempMat, transpMat, f_Mat, f_TempMat, diagMat;
	int *size, *displacement;
	int matrixSize, matrixTempSize, n, rank, mpiSize;
	double h = 1.0;
	init_app(&argc, argv, &rank, &mpiSize);

	if (argc < 2 && argv[1] > 4 && argv[1]%2 == 0){
		printf("Need a problem size! And it should be a power of 2 greater than 4!\n");
		return;
	}

	/*			Initializing variables		*/

	//We need to set up the rest of the variables and use the data given to us by splitVector in a smart manner.
	//Do something with splitVector and *size and *displacement
	splitVector(matrixSize, mpiSize, &size, &displacement);

	n = atoi(argv[1]);
	matrixSize = n-1;
	matrixTempSize = matrixSize*4;
	h /= n;

	/*			Initializing structures		*/
	//Use a MPI communicator? We have two defined in the .h file. Maybe make createMatrix() set that for us.
	diagMat = createVector(m);
	matrix = createMatrix(m, m);
	transpMat = createMatrix(m, m);
	f_TempMat = createVector(matrixTempSize);
	diagMat->data = (double*) calloc(m, sizeof(double));
	f_TempMat->data = (double*) malloc(matrixTempSize*sizeof(double));


	return 0;
}
