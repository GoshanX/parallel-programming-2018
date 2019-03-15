#include <stdio.h>
#include <stdlib.h>

#define MTRX_SIZE 8
#define BLCK_SIZE 2

typedef double matrix[MTRX_SIZE][MTRX_SIZE];
typedef double block[BLCK_SIZE][BLCK_SIZE];



void initialization_matrix(matrix m1, matrix m2, matrix m3, matrix m4)
{
	for (int i = 0; i < MTRX_SIZE; i++)
		for (int j = 0; j < MTRX_SIZE; j++)
		{
			m1[i][j] = rand() % 9;
			m2[i][j] = rand() % 9;
			m3[i][j] = 0;
			m4[i][j] = 0;
		}
}
void basic_multiplication(matrix m1, matrix m2, matrix m4)
{
	int i, j, k;
	for (i = 0; i < MTRX_SIZE; i++) {
		for (j = 0; j < MTRX_SIZE; j++) {
			double sum = 0;
			for (k = 0; k < MTRX_SIZE; k++) {
				sum += m1[i][k] * m2[k][j];
			}
			m4[i][j] += sum;
		}
	}
}

int matrix_equal(matrix m1, matrix m2)
{
	int equal = 1;
	for (int i = 0; i < MTRX_SIZE; i++)
		for (int j = 0; j < MTRX_SIZE; j++)
			if (abs(m1[i][j] - m2[i][j]) > 0.00001)
				equal = 0;
	return equal;
}

void matrix_print(matrix m) {
	int i, j;
	for (i = 0; i < MTRX_SIZE; i++) {
		for (j = 0; j < MTRX_SIZE; j++) {
			printf("%f ", m[i][j]);
			if (j == MTRX_SIZE - 1) {
				printf("\n");
			}
		}
	}
	printf("\n");
}

void block_add(block b1, block b2) {
	int i, j;
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			b1[i][j] += b2[i][j];
		}
	}
}

void block_mul(block b1, block b2, block b3) {
	int i, j, k;
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			b3[i][j] = 0;
			for (k = 0; k < BLCK_SIZE; k++) {
				b3[i][j] += b1[i][k] * b2[k][j];
			}
		}
	}
}

void block_get(block b, matrix m, int k, int l) {
	int i, j;
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			b[i][j] = m[BLCK_SIZE * k + i][BLCK_SIZE * l + j];
		}
	}
	/*
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			printf("%f", b[i][j]);
			if (j == BLCK_SIZE - 1) {
				printf("\n");
			}
		}
	}*/
}

void block_set(block b, matrix m, int k, int l) {
	int i, j;
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			m[BLCK_SIZE * k + i][BLCK_SIZE * l + j] = b[i][j];
		}
	}
	
}



void block_reset(block b) {
	int i, j;
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			b[i][j] = 0;
		}
	}
}

void matrix_mul(matrix m1, matrix m2, matrix m3) {
	int i, j, k;
	block b1, b2, b_mul, b_sum;
	for (i = 0; i < MTRX_SIZE / BLCK_SIZE; i++) {
		for (j = 0; j < MTRX_SIZE / BLCK_SIZE; j++) {
			block_reset(b_sum);
			for (k = 0; k < MTRX_SIZE / BLCK_SIZE; k++) {
				block_get(b1, m1, i, k);
				block_get(b2, m2, k, j);
				block_mul(b1, b2, b_mul);
				block_add(b_sum, b_mul);
			}
			block_set(b_sum, m3, i, j);
		}
	}
}

matrix *matrix_new(void) {
	return malloc(sizeof(matrix));
}



int main() {

	if (MTRX_SIZE%BLCK_SIZE != 0) {
		printf("Invalid values\n");
		system("pause");
		return 1;
	}
	
	matrix *m1, *m2, *m3, *m4;

	m1 = matrix_new();
	m2 = matrix_new();
	m3 = matrix_new();
	m4 = matrix_new();
	initialization_matrix(m1, m2, m3, m4);

	
	matrix_mul(*m1, *m2, *m3);
	basic_multiplication(m1, m2, m4);
	if (MTRX_SIZE < 10) {
		printf(" Matrix 1\n");
		matrix_print(*m1);
		printf(" Matrix 2\n");
		matrix_print(*m2);
		printf("Fox\n");
		matrix_print(*m3);
		printf("Basic\n");
		matrix_print(*m4);
	}
	if (matrix_equal(m3, m4)) printf("Answer is correct\n");
	else printf("Answer is incorrect\n");

	free(m1);
	free(m2);
	free(m3);
	system("pause");
}