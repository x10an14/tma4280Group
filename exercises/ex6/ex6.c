#include "ex6.h"

Vector createVector(int len){
	Vector result = (Vector) calloc(1, sizof(vector_t));
	result->len = len;
	result->data = (double*) calloc(len, sizeof(double));
}

Matrix createMatrix(int rows, int cols){
	Matrix result = (Matrix) calloc(1, sizof(matrix_t));
	result->as_vec = createVector(rows*cols);

}

//Copied from common.c, minor edits.
void splitVector(int globLen, int size, int** len, int** displ){
	*len = calloc(size,sizeof(int));
	*displ = calloc(size,sizeof(int));

	for (i=0;i<size;++i){
		(*len)[i] = globLen/size;

		if (globLen % size && i >= (size - globLen % size)){
			(*len)[i]++;
		}

		if (i < size-1){
			(*displ)[i+1] = (*displ)[i]+(*len)[i];
		}
	}
}

/* DO NOT USE THE COMMONS LIBRARY!
*IF ANYTHING IN THE COMMONS LIBRARY IS OF USE, COPY IT OVER.
*THE COMMONS LIBRARY HAS TOO MUCH "UNNECESSARY DATA/JUNK"
*/

int int main(int argc, char *argv[]){
	Matrix matrix, tempMatrix, transposedMatrix;
	int *size, *displacement;
	int matrixSize, matrixTempSize, n, rank = 0, mpiSize = 1;
	init_app(&argc, argv);

	if (argc < 2 && argv[1] > 4 && argv[1]%2 == 0){
		printf("Need a problem size! And it has to be a power of 2 greater than 4!\n");
		return;
	}

	//Initializing variables
	n = argv[1]
	matrixSize = n-1;
	matrixTempSize = matrixSize*4;

	//Do something with splitVector and *size and *displacement
	splitVector(matrixSize, mpiSize, &size, &displacement);

	//We need to set up the rest of the variables and use the data given to us by splitVector in a smart manner.

	//Use the MPI Communicatior "WorldComm"(maybe?), allocate data for the matrix, and whether we pad the data
	//^-- Explanation of the last three inputs of the below calls.
	matrix = createMatrixMPI(matrixSize, &WorldComm, 1, 0);
	tempMatrix = createMatrixMPI(matrixTempSize, &WorldComm, 1, 0);
	transposedMatrix = createMatrixMPI(matrixSize, &WorldComm, 1, 0);

	return 0;
}
