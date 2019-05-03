//подключение библиотеки TBB
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include <tbb/tick_count.h>
using namespace tbb;

#include "time.h"
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
using namespace std;

#define MTRX_SIZE 512
#define BLCK_SIZE 256





void initialization_matrix(double **m1, double **m2, double **m3, double **m4, double **m5, double **b1, double **b2, double **b_mul, double **b_sum)
{
	for (int i = 0; i < MTRX_SIZE; i++)
		for (int j = 0; j < MTRX_SIZE; j++)
		{
			m1[i][j] = double(rand()) / RAND_MAX;
			m2[i][j] = double(rand()) / RAND_MAX;
			m3[i][j] = 0;
			m4[i][j] = 0;
			m5[i][j] = 0;
			
		}
	for (int i = 0; i < BLCK_SIZE; i++)
		for (int j = 0; j < BLCK_SIZE; j++)
		{
			b1[i][j] = 0;
			b2[i][j] = 0;
			b_mul[i][j] = 0;
			b_sum[i][j] = 0;

		}

}
void basic_multiplication(double **m1, double **m2, double **m4)
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

int matrix_equal(double **m1, double **m2)
{
	int equal = 1;
	for (int i = 0; i < MTRX_SIZE; i++)
		for (int j = 0; j < MTRX_SIZE; j++)
			if (abs((double)(m1[i][j] - m2[i][j])) > 100)
				equal = 0;
	return equal;
}

void matrix_print(double **m) {
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

void block_add(double **b1, double **b2) {
	int i, j;
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			b1[i][j] += b2[i][j];
		}
	}
}

void block_mul(double **b1, double **b2, double **b3) {
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

void block_get(double **b, double **m, int k, int l) {
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

void block_set(double **b, double **m, int k, int l) {
	int i, j;
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			m[BLCK_SIZE * k + i][BLCK_SIZE * l + j] = b[i][j];
		}
	}

}
void block_reset(double **b) {
	int i, j;
	for (i = 0; i < BLCK_SIZE; i++) {
		for (j = 0; j < BLCK_SIZE; j++) {
			b[i][j] = 0;
		}
	}
}




void matrix_mul(double **m1, double **m2, double **m3,double **b1, double **b2, double **b_mul, double **b_sum) {
	int i, j, k;
	
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
class TBB_Fox {
	double ** const m1;
	double ** const m2;
	double ** const m3;
	double ** const b1;
	double ** const b2;
	double ** const b_mul;
	double ** const b_sum;
public:
	TBB_Fox(double **_m1, double **_m2, double **_m3, double **_b1, double **_b2, double **_b_mul, double **_b_sum) :m1(_m1), m2(_m2), m3(_m3), b1(_b1), b2(_b2), b_mul(_b_mul), b_sum(_b_sum)
	{}


	void operator()(const blocked_range<int>& r)  const {
		int i, j, k;


		int begin = r.begin();
		int end = r.end();


		for (i = begin; i < end; i++) {
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

};




int main() {

	if (MTRX_SIZE%BLCK_SIZE != 0) {
		printf("Invalid values\n");
		system("pause");
		return 1;
	}



	double** m1 = new double*[MTRX_SIZE];
	for (size_t index = 0; index < MTRX_SIZE; ++index) {
		m1[index] = new double[MTRX_SIZE];
	}
	double** m2 = new double*[MTRX_SIZE];
	for (size_t index = 0; index < MTRX_SIZE; ++index) {
		m2[index] = new double[MTRX_SIZE];
	}
	double** m3 = new double*[MTRX_SIZE];
	for (size_t index = 0; index < MTRX_SIZE; ++index) {
		m3[index] = new double[MTRX_SIZE];
	}
	double** m4 = new double*[MTRX_SIZE];
	for (size_t index = 0; index < MTRX_SIZE; ++index) {
		m4[index] = new double[MTRX_SIZE];
	}
	double** m5 = new double*[MTRX_SIZE];
	for (size_t index = 0; index < MTRX_SIZE; ++index) {
		m5[index] = new double[MTRX_SIZE];
	}
	double** b1 = new double*[BLCK_SIZE];
	for (size_t index = 0; index < BLCK_SIZE; ++index) {
		b1[index] = new double[BLCK_SIZE];
	}
	double** b2 = new double*[BLCK_SIZE];
	for (size_t index = 0; index < BLCK_SIZE; ++index) {
		b2[index] = new double[BLCK_SIZE];
	}
	double** b_mul = new double*[BLCK_SIZE];
	for (size_t index = 0; index < BLCK_SIZE; ++index) {
		b_mul[index] = new double[BLCK_SIZE];
	}
	double** b_sum = new double*[BLCK_SIZE];
	for (size_t index = 0; index < BLCK_SIZE; ++index) {
		b_sum[index] = new double[BLCK_SIZE];
	}
	initialization_matrix(m1, m2, m3, m4,m5,b1,b2,b_mul,b_sum);
	
	 
	TBB_Fox s(m1, m2, m3,b1,b2,b_mul,b_sum);

	task_scheduler_init init(2);

	tick_count start = tick_count::now();
	parallel_for(blocked_range<int>(0, (MTRX_SIZE / BLCK_SIZE), (MTRX_SIZE / BLCK_SIZE)/2), s);
	tick_count finish = tick_count::now();

	
	printf("TBB Time = %f seconds\n", (finish - start).seconds());


	tick_count start1 = tick_count::now();
	matrix_mul(m1, m2, m4, b1, b2, b_mul, b_sum);
	tick_count finish1 = tick_count::now();
	basic_multiplication(m1, m2, m5);
	init.terminate();
	printf("Seq Time = %f seconds\n", (finish1 - start1).seconds());
	printf("Acceleration = %f seconds\n", (finish1 - start1).seconds()/(finish - start).seconds());

	if (MTRX_SIZE < 10) {
		printf(" Matrix 1\n");
		matrix_print(m1);
		printf(" Matrix 2\n");
		matrix_print(m2);
		printf("FoxTBB\n");
		matrix_print(m3);
		printf("Fox\n");
		matrix_print(m4);
		printf("Basic\n");
		matrix_print(m5);
	}
	if (matrix_equal(m3, m5)) printf("Answer is correct\n");
	else printf("Answer is incorrect\n");

	free(m1);
	free(m2);
	free(m3);
	free(m4);
	free(b1);
	free(b2);
	free(b_mul);
	free(b_sum);
	system("pause");
}
