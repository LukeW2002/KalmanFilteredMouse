#ifndef LINEAR_ALG
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <errno.h>


		/* MATHS STRUCTURES DEFINITIONS */


typedef struct {
 double *data;
 int rows;
 int cols; 
}Matrix;

typedef struct {
 double *data;
 int dimension;
}Vector;

void *allocateMemory( size_t size);
Vector *allocateVector(int dimension);
Matrix *allocateMatrix(int rows, int cols);
void freeMatrix(Matrix *matrix);
void freeVector(Vector *vector);
void printMatrix(Matrix *matrix);
void printVector(Vector *vector);

		/* MATRIX OPERATIONS */
	
void vectorSubtract(Vector *vector1, Vector *vector2, Vector *result); // Define with vector2 - vector1
void vectorAdd(Vector *vector1, Vector *vector2, Vector *result);
void setZeroVector(Vector *vector);

void vectorMatrixMultiply(Matrix *matrix, Vector *vector, Vector *result);
void matrixMultiply( Matrix *matrix1, Matrix *matrix2, Matrix *result); // Defined with right hand multiplication of matrix1 to matrix2
void matrixAdd( Matrix *matrix, Matrix *matrix2, Matrix *result);
void matrixSubtract( Matrix *matrix, Matrix *matrix2, Matrix *result);
void matrixTranspose(Matrix *matrix, Matrix *result);
void matrixCopy(Matrix *matrix, Matrix *result);
void setIdentity(Matrix *matrix);
void setZero(Matrix *matrix);

			/* MATRIX INVERSION */
void choleskyDecomp(Matrix *matrix, Matrix *result);
void forwardSubsitution(Matrix *matrix, Vector *e, Vector *X);
void backSubsitution(Matrix *matrix, Vector *e, Vector *X);
void matrixInverse(Matrix *A, Matrix *result);


#endif

