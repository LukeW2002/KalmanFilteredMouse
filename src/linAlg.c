#include "header/linAlg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void printMatrix(Matrix *matrix) 
{
	for (int i = 0; i < matrix->rows; i++)
	{
		for (int j = 0; j < matrix->cols; j++)
		{
			printf("%lf\t",matrix->data[i*(matrix->cols) + j] );
		}
		printf("\n");
	}
	printf("\n");
}

void printVector(Vector *vector)
{
	for (int j = 0; j < vector->dimension; j++)
	{
		printf("%lf\n", vector->data[j]);
	}
	printf("\n");
}



void *allocateMemory( size_t size)
{
	void *ptr = (double *)malloc(size);
	if (ptr == NULL)
	{
		perror("Failed to allocate memory");
		exit(-1);
	}
	return ptr;
}


Vector *allocateVector(int dimension)
{
	Vector *vector = (Vector *)malloc(sizeof(Vector));
	vector->dimension = dimension;	
	vector->data = allocateMemory(dimension * sizeof(double)); //Allocates an array to memory
	//for (int i = 0; i< dimension; i++)
	//{
	//	vector->data[i] = 0;
	//}
	return vector;
}


Matrix *allocateMatrix(int rows, int cols)
{	
	Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
	matrix->rows = rows;
	matrix->cols = cols;
	matrix->data = allocateMemory(rows*cols*sizeof(double));
	//for (int i = 0; i< rows*cols; i++)
	//{
	//	matrix->data[i] = 0;
	//}
	return matrix;
}

void freeMatrix(Matrix *matrix)
{
	free(matrix->data);
	free(matrix);
}
void freeVector(Vector *vector)
{
	free(vector->data);
	free(vector);
}

		/* MATRIX OPERATIONS */
void setIdentity(Matrix *matrix)
{
	for (int i = 0; i< matrix->rows; i++)
	{
		for (int j = 0; j < matrix->cols; j++ )
		{
			matrix->data[i * (matrix->cols) + j] = 0;
			matrix->data[i * (matrix->cols) + i] = 1;
		}
	}
}

void setZero(Matrix *matrix)
{
	for (int i = 0; i< matrix->rows; i++)
	{
		for (int j = 0; j < matrix->cols; j++ )
		{
			matrix->data[i * (matrix->cols) + j] = 0;
		}
	}
}

void matrixMultiply(Matrix *matrix1, Matrix *matrix2, Matrix *result )
{
	if (matrix1->cols != matrix2->rows)
	{
		printf("(cols, rows)\n");
		printf("Matrix1:\t (%i, %i)\n Matrix2:\t (%i, %i)", matrix1->cols, matrix1->rows, matrix2->cols, matrix2->rows);
		printf("\n ERROR matrixMultiply \n INCORRECT DIMENSION \n");
		exit(1);
	}
	for (int i = 0; i < matrix1->rows; ++i)
	{
		for (int j = 0; j < matrix2->cols; ++j) 
		{
			result->data[i*(matrix2->cols) +j] = 0;
			for (int l = 0; l < matrix1->cols; ++l)
			{
				result->data[i*(matrix2->cols) + j] += matrix1->data[i*(matrix1->cols) + l] * matrix2->data[l*(matrix2->cols)+j];

			}
		}
	}
}

void setZeroVector(Vector *vector)
{
	for (int i = 0; i< vector->dimension; i++)
	{
		vector->data[i] = 0;
	}
}

void vectorMatrixMultiply( Matrix *matrix, Vector *vector, Vector *result )
{
	if (matrix->cols != vector->dimension)
	{
		printf("\n ERROR vectorMultiply \n INCORRECT DIMENSION \n");
		exit(1);
	}
	for (int i = 0; i < matrix->rows; i++)
	{
		result->data[i] = 0;
		for(int j = 0; j < matrix->cols; j++)
		{
			result->data[i] += matrix->data[i*(matrix->cols) + j] * vector->data[j];
		}
	}
}

void vectorAdd(Vector *vector1, Vector *vector2, Vector *result)
{
	if (vector1->dimension != vector2->dimension)
	{
		printf("\n ERROR vectorAdd \n INCORRECT DIMENSION  \n");
		printf("vector1 dim %i\n vector2 dim %i\n", vector1->dimension, vector2->dimension);
		exit(1);
	}
	for (int i = 0; i < vector1->dimension; i++)
	{
		result->data[i] = vector2->data[i] + vector1->data[i];
	}
}

void vectorSubtract(Vector *vector1, Vector *vector2, Vector *result)
{
	if (vector1->dimension != vector2->dimension)
	{
		printf("\n ERROR vectorSubtract \n INCORRECT DIMENSION  \n");
		printf("vector1 dim %i\n vector2 dim %i\n", vector1->dimension, vector2->dimension);
		exit(1);
	}
	for (int i = 0; i < vector1->dimension; i++)
	{
		result->data[i] = vector2->data[i] - vector1->data[i];
	}
}

void matrixAdd( Matrix *matrix1, Matrix *matrix2, Matrix *result)
{
	if (matrix1->rows != matrix2->rows && matrix1->cols != matrix2->cols)
	{
		printf("\n ERROR matrixAdd \n INCORRECT DIMENSION \n");
		exit(1);
	}

	for (int i = 0; i< matrix1->rows; ++i)
	{
		for (int j = 0; j< matrix2->cols; ++j)
		{
			result->data[i*(matrix1->rows) +j] = matrix1->data[i* (matrix1->rows) +j ] + matrix2->data[i*(matrix2->rows) + j];
		}
	}
}


void matrixSubtract( Matrix *matrix1, Matrix *matrix2, Matrix *result)
{
	if (matrix1->rows != matrix2->rows && matrix1->cols != matrix2->cols)
	{
		printf("\n ERROR matrixAdd \n INCORRECT DIMENSION \n");
		exit(1);
	}

	for (int i = 0; i< matrix1->rows; ++i)
	{
		for (int j = 0; j< matrix2->cols; ++j)
		{
			result->data[i*(matrix1->rows) +j] = matrix1->data[i* (matrix1->rows) +j ] - matrix2->data[i*(matrix2->rows) + j];
		}
	}
}

void matrixCopy(Matrix *matrix, Matrix *result)
{
	for (int i = 0; i < matrix->rows; i++)
	{
		for (int j = 0; j < matrix->cols; j++)
		{
			result->data[i * (matrix->cols) + j] = matrix->data[i * (matrix->cols) + j];
		}
	}
}

void matrixTranspose(Matrix *matrix, Matrix *result)
{
	if (matrix->cols != result->rows || matrix->rows != result->cols)
	{
		printf("error matrix transpose\n");
	}
	for ( int i = 0; i < result->rows; i++)
	{
		for (int  j = 0; j < result->cols; j++)
		{
			result->data[i * result->cols + j] = matrix->data[j * matrix->cols + i];
		}
	}
}

void matrixTransposeSqaure(Matrix *matrix, Matrix *result)
{
	double *transpose;
	double *ptrMatrix = matrix->data;

	for ( int i =0; i < matrix->rows; i++)
	{
		transpose = &result->data[i];
		for (int  j = 0; j < matrix->cols; j++)
		{
			*transpose = *ptrMatrix;
			ptrMatrix++;
			transpose += matrix->rows;
		}
	}
}

			/* MATRIX INVERSION */
void choleskyDecomp(Matrix *matrix, Matrix *result)
{
	if (matrix->rows != matrix-> cols)
	{
		printf("\n ERROR choleskyDecomp \n INCORRECT DIMENSION \n");
		exit(1);
	}
	double sum;
	int n = matrix->rows;
	setZero(result);
	for(int i = 0; i < n; i++)
	{
		for (int j = 0; j < (i+1) ; j++)
		{
			sum = 0.0;
			for (int k = 0; k < j; k++)
			{
				sum += result->data[ i * n +k] * result->data[j * n + k];
			}
			if ( i == j)
			{
				result->data[i * n + j] = sqrt(matrix->data[i * n + i] - sum);
			}
			else
			{
				result->data[i * n + j] = (1.0 / result->data[j * n + j]) * (matrix->data[ i * n + j] - sum);
			}
		}
	}
}


void forwardSubsitution(Matrix *matrix, Vector *e, Vector *X)
{
	if (matrix->rows != matrix-> cols)
	{
		printf("\n ERROR forwardSubsitution \n INCORRECT DIMENSION \n");
		exit(1);
	}
	for (int i = 0; i < matrix->rows; i++)
	{
		double sum = e->data[i];
		for ( int j = 0; j <= (i-1); j++)
		{
			sum -= matrix->data[i*(matrix->cols) +j ] * X->data[j];
		}
		X->data[i] = sum / matrix->data[i*(matrix->cols)+i];
	}	
}

void backSubsitution(Matrix *matrix, Vector *e, Vector *X)
{
	if (matrix->rows != matrix-> cols)
	{
		printf("\n ERROR backSubsitution \n INCORRECT DIMENSION \n");
		exit(1);
	}
	for (int i = (matrix->rows) - 1; i >= 0; i--)
	{
		double sum = e->data[i];
		for (int j = i+1; j < (matrix->cols); j++)
		{
			sum -= matrix->data[i*(matrix->cols) + j] * X->data[j];
		}
		X->data[i] = sum / matrix->data[i*(matrix->cols) + i];
	}	
}

void matrixInverse(Matrix *matrix, Matrix *result)
{
	if (matrix->rows != matrix-> cols)
	{
		printf("\n ERROR matrixInverse \n INCORRECT DIMENSION \n");
		printf("Matrix:\t (%i, %i)\n", matrix->cols, matrix->rows);
		exit(1);
	}
	Matrix  *L = allocateMatrix(matrix->rows, matrix->cols);
	Matrix *transposeL = allocateMatrix(matrix->rows, matrix->cols);

	choleskyDecomp(matrix, L);
	matrixTranspose(L, transposeL);
	//printf("L\n");
	//printMatrix(L);
	//printf("Lt \n");
	//printMatrix(transposeL);

	for (int j = 0; j < matrix->rows; j++)
	{
		Vector *e = allocateVector(matrix->rows); 
		Vector *x = allocateVector(matrix->rows); 
		Vector *y = allocateVector(matrix->rows); 
		setZeroVector(e);
		setZeroVector(x);
		setZeroVector(y);

		e->data[j] = 1.0;
		forwardSubsitution(L, e, x);
		backSubsitution(transposeL, x, y);
		for (int k = (matrix->rows)-1; k >=0 ; k--)
		{
			result->data[k*(matrix->rows)+j] = y->data[k];
			//printf("result: %fl\n y: %fl\n", result->data[k*(matrix->rows)+j], y->data[k]);
		}
		freeVector(e);
		freeVector(x);
		freeVector(y);
	}
	freeMatrix(L);
	freeMatrix(transposeL);
}

