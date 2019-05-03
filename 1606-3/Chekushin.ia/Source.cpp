////подключение библиотеки TBB
//#include "tbb/task_scheduler_init.h"
//#include "tbb/parallel_for.h"
//#include "tbb/parallel_reduce.h"
//#include "tbb/blocked_range.h"
//#include <tbb/tick_count.h>
//using namespace tbb;
//
//#include "time.h"
//#include <vector>
//#include <random>
//#include <ctime>
//#include <iostream>
//using namespace std;
//
//#define MTRX_SIZE 512
//#define BLCK_SIZE 256
//
//
//
//
//
//void initialization_matrix(double **m1, double **m2, double **m3, double **m4, double **m5, double **b1, double **b2, double **b_mul, double **b_sum)
//{
//	for (int i = 0; i < MTRX_SIZE; i++)
//		for (int j = 0; j < MTRX_SIZE; j++)
//		{
//			m1[i][j] = double(rand()) / RAND_MAX;
//			m2[i][j] = double(rand()) / RAND_MAX;
//			m3[i][j] = 0;
//			m4[i][j] = 0;
//			m5[i][j] = 0;
//			
//		}
//	for (int i = 0; i < BLCK_SIZE; i++)
//		for (int j = 0; j < BLCK_SIZE; j++)
//		{
//			b1[i][j] = 0;
//			b2[i][j] = 0;
//			b_mul[i][j] = 0;
//			b_sum[i][j] = 0;
//
//		}
//
//}
//void basic_multiplication(double **m1, double **m2, double **m4)
//{
//	int i, j, k;
//	for (i = 0; i < MTRX_SIZE; i++) {
//		for (j = 0; j < MTRX_SIZE; j++) {
//			double sum = 0;
//			for (k = 0; k < MTRX_SIZE; k++) {
//				sum += m1[i][k] * m2[k][j];
//			}
//			m4[i][j] += sum;
//		}
//	}
//}
//
//int matrix_equal(double **m1, double **m2)
//{
//	int equal = 1;
//	for (int i = 0; i < MTRX_SIZE; i++)
//		for (int j = 0; j < MTRX_SIZE; j++)
//			if (abs((double)(m1[i][j] - m2[i][j])) > 100)
//				equal = 0;
//	return equal;
//}
//
//void matrix_print(double **m) {
//	int i, j;
//	for (i = 0; i < MTRX_SIZE; i++) {
//		for (j = 0; j < MTRX_SIZE; j++) {
//			printf("%f ", m[i][j]);
//			if (j == MTRX_SIZE - 1) {
//				printf("\n");
//			}
//		}
//	}
//	printf("\n");
//}
//
//void block_add(double **b1, double **b2) {
//	int i, j;
//	for (i = 0; i < BLCK_SIZE; i++) {
//		for (j = 0; j < BLCK_SIZE; j++) {
//			b1[i][j] += b2[i][j];
//		}
//	}
//}
//
//void block_mul(double **b1, double **b2, double **b3) {
//	int i, j, k;
//	for (i = 0; i < BLCK_SIZE; i++) {
//		for (j = 0; j < BLCK_SIZE; j++) {
//			b3[i][j] = 0;
//			for (k = 0; k < BLCK_SIZE; k++) {
//				b3[i][j] += b1[i][k] * b2[k][j];
//			}
//		}
//	}
//}
//
//void block_get(double **b, double **m, int k, int l) {
//	int i, j;
//	for (i = 0; i < BLCK_SIZE; i++) {
//		for (j = 0; j < BLCK_SIZE; j++) {
//			b[i][j] = m[BLCK_SIZE * k + i][BLCK_SIZE * l + j];
//		}
//	}
//	/*
//	for (i = 0; i < BLCK_SIZE; i++) {
//		for (j = 0; j < BLCK_SIZE; j++) {
//			printf("%f", b[i][j]);
//			if (j == BLCK_SIZE - 1) {
//				printf("\n");
//			}
//		}
//	}*/
//}
//
//void block_set(double **b, double **m, int k, int l) {
//	int i, j;
//	for (i = 0; i < BLCK_SIZE; i++) {
//		for (j = 0; j < BLCK_SIZE; j++) {
//			m[BLCK_SIZE * k + i][BLCK_SIZE * l + j] = b[i][j];
//		}
//	}
//
//}
//void block_reset(double **b) {
//	int i, j;
//	for (i = 0; i < BLCK_SIZE; i++) {
//		for (j = 0; j < BLCK_SIZE; j++) {
//			b[i][j] = 0;
//		}
//	}
//}
//
//
//
//
//void matrix_mul(double **m1, double **m2, double **m3,double **b1, double **b2, double **b_mul, double **b_sum) {
//	int i, j, k;
//	
//	for (i = 0; i < MTRX_SIZE / BLCK_SIZE; i++) {
//		for (j = 0; j < MTRX_SIZE / BLCK_SIZE; j++) {
//			block_reset(b_sum);
//			for (k = 0; k < MTRX_SIZE / BLCK_SIZE; k++) {
//				block_get(b1, m1, i, k);
//				block_get(b2, m2, k, j);
//				block_mul(b1, b2, b_mul);
//				block_add(b_sum, b_mul);
//			}
//			block_set(b_sum, m3, i, j);
//		}
//	}
//}
//class TBB_Fox {
//	double ** const m1;
//	double ** const m2;
//	double ** const m3;
//	double ** const b1;
//	double ** const b2;
//	double ** const b_mul;
//	double ** const b_sum;
//public:
//	TBB_Fox(double **_m1, double **_m2, double **_m3, double **_b1, double **_b2, double **_b_mul, double **_b_sum) :m1(_m1), m2(_m2), m3(_m3), b1(_b1), b2(_b2), b_mul(_b_mul), b_sum(_b_sum)
//	{}
//
//
//	void operator()(const blocked_range<int>& r)  const {
//		int i, j, k;
//
//
//		int begin = r.begin();
//		int end = r.end();
//
//
//		for (i = begin; i < end; i++) {
//			for (j = 0; j < MTRX_SIZE / BLCK_SIZE; j++) {
//				block_reset(b_sum);
//				for (k = 0; k < MTRX_SIZE / BLCK_SIZE; k++) {
//					block_get(b1, m1, i, k);
//					block_get(b2, m2, k, j);
//					block_mul(b1, b2, b_mul);
//					block_add(b_sum, b_mul);
//				}
//				block_set(b_sum, m3, i, j);
//			}
//		}
//
//	}
//
//};
//
//
//
//
//int main() {
//
//	if (MTRX_SIZE%BLCK_SIZE != 0) {
//		printf("Invalid values\n");
//		system("pause");
//		return 1;
//	}
//
//
//
//	double** m1 = new double*[MTRX_SIZE];
//	for (size_t index = 0; index < MTRX_SIZE; ++index) {
//		m1[index] = new double[MTRX_SIZE];
//	}
//	double** m2 = new double*[MTRX_SIZE];
//	for (size_t index = 0; index < MTRX_SIZE; ++index) {
//		m2[index] = new double[MTRX_SIZE];
//	}
//	double** m3 = new double*[MTRX_SIZE];
//	for (size_t index = 0; index < MTRX_SIZE; ++index) {
//		m3[index] = new double[MTRX_SIZE];
//	}
//	double** m4 = new double*[MTRX_SIZE];
//	for (size_t index = 0; index < MTRX_SIZE; ++index) {
//		m4[index] = new double[MTRX_SIZE];
//	}
//	double** m5 = new double*[MTRX_SIZE];
//	for (size_t index = 0; index < MTRX_SIZE; ++index) {
//		m5[index] = new double[MTRX_SIZE];
//	}
//	double** b1 = new double*[BLCK_SIZE];
//	for (size_t index = 0; index < BLCK_SIZE; ++index) {
//		b1[index] = new double[BLCK_SIZE];
//	}
//	double** b2 = new double*[BLCK_SIZE];
//	for (size_t index = 0; index < BLCK_SIZE; ++index) {
//		b2[index] = new double[BLCK_SIZE];
//	}
//	double** b_mul = new double*[BLCK_SIZE];
//	for (size_t index = 0; index < BLCK_SIZE; ++index) {
//		b_mul[index] = new double[BLCK_SIZE];
//	}
//	double** b_sum = new double*[BLCK_SIZE];
//	for (size_t index = 0; index < BLCK_SIZE; ++index) {
//		b_sum[index] = new double[BLCK_SIZE];
//	}
//	initialization_matrix(m1, m2, m3, m4,m5,b1,b2,b_mul,b_sum);
//	
//	 
//	TBB_Fox s(m1, m2, m3,b1,b2,b_mul,b_sum);
//
//	task_scheduler_init init(2);
//
//	tick_count start = tick_count::now();
//	parallel_for(blocked_range<int>(0, (MTRX_SIZE / BLCK_SIZE), (MTRX_SIZE / BLCK_SIZE)/2), s);
//	tick_count finish = tick_count::now();
//
//	
//	printf("TBB Time = %f seconds\n", (finish - start).seconds());
//
//
//	tick_count start1 = tick_count::now();
//	matrix_mul(m1, m2, m4, b1, b2, b_mul, b_sum);
//	tick_count finish1 = tick_count::now();
//	basic_multiplication(m1, m2, m5);
//	init.terminate();
//	printf("Seq Time = %f seconds\n", (finish1 - start1).seconds());
//	printf("Acceleration = %f seconds\n", (finish1 - start1).seconds()/(finish - start).seconds());
//
//	if (MTRX_SIZE < 10) {
//		printf(" Matrix 1\n");
//		matrix_print(m1);
//		printf(" Matrix 2\n");
//		matrix_print(m2);
//		printf("FoxTBB\n");
//		matrix_print(m3);
//		printf("Fox\n");
//		matrix_print(m4);
//		printf("Basic\n");
//		matrix_print(m5);
//	}
//	if (matrix_equal(m3, m5)) printf("Answer is correct\n");
//	else printf("Answer is incorrect\n");
//
//	free(m1);
//	free(m2);
//	free(m3);
//	free(m4);
//	free(b1);
//	free(b2);
//	free(b_mul);
//	free(b_sum);
//	system("pause");
//}
//
////подключение библиотеки TBB
////#include "tbb/task_scheduler_init.h"
////#include "tbb/parallel_for.h"
////#include "tbb/parallel_reduce.h"
////#include "tbb/blocked_range.h"
////#include <tbb/tick_count.h>
////
////#include "time.h"
////#include <vector>
////#include <random>
////#include <ctime>
////#include <iostream>
////using namespace std;
////using namespace tbb;
////
////#define MTRX_SIZE 512//1024//2048
////#define BLCK_SIZE 128//128//256
////
////void initialization_matrix(int* matrA, int* matrB, int* matrC, int* matrCchecked)
////{
////	for (int i = 0; i < MTRX_SIZE; i++)
////		for (int j = 0; j < MTRX_SIZE; j++)
////		{
////			matrA[i*MTRX_SIZE + j] = rand() % 9;
////			matrB[i*MTRX_SIZE + j] = rand() % 9;
////			matrC[i*MTRX_SIZE + j] = 0;
////			matrCchecked[i*MTRX_SIZE + j] = 0;
////		}
////}
////
////int matrix_equal(int *matrA, int *matrB)
////{
////	int equal = 1;
////	for (int i = 0; i < MTRX_SIZE; i++)
////		for (int j = 0; j < MTRX_SIZE; j++)
////			if (abs((int)(matrA[i*MTRX_SIZE + j] - matrB[i*MTRX_SIZE + j])) > 0.01)
////				equal = 0;
////	return equal;
////}
////
////void print_matrix(int *matr, int rows, int columns) {
////	for (int i = 0; i < rows; i++) {
////		for (int j = 0; j < columns; j++) {
////			cout << matr[i * columns + j] << " ";
////		}
////		cout << endl;
////	}
////}
////
////class TBB_Fox {
////	int *matrA;
////	int *matrB;
////	int *matrC;
////public:
////	TBB_Fox(int *matrA, int *matrB, int *matrC) :matrA(matrA), matrB(matrB), matrC(matrC)
////	{}
////
////	void operator()(const blocked_range<int>& r)  const {
////		int begin = r.begin();
////		int end = r.end();
////		cout << begin << " " << end << endl;
////
////		for (int jj = 0; jj < MTRX_SIZE; jj += BLCK_SIZE)
////			for (int kk = 0; kk < MTRX_SIZE; kk += BLCK_SIZE)
////				for (int ii = 0; ii < MTRX_SIZE; ii += BLCK_SIZE)
////					for (int j = 0; j < BLCK_SIZE; j++)
////						for (int k = 0; k < BLCK_SIZE; k++)
////							for (int i = 0; i < BLCK_SIZE; i++)
////								matrC[(j + jj) * MTRX_SIZE + (ii + i)] += matrA[(jj + j) * MTRX_SIZE + (kk + k)] * matrB[(kk + k)*MTRX_SIZE + (ii + i)];
////	}
////};
////
////void matrix_mul(int* matrA, int* matrB, int* matrC)
////{
////	for (int jj = 0; jj < MTRX_SIZE; jj += BLCK_SIZE)
////		for (int kk = 0; kk < MTRX_SIZE; kk += BLCK_SIZE)
////			for (int ii = 0; ii < MTRX_SIZE; ii += BLCK_SIZE)
////				for (int j = 0; j < BLCK_SIZE; j++)
////					for (int k = 0; k < BLCK_SIZE; k++)
////						for (int i = 0; i < BLCK_SIZE; i++)
////							matrC[(j + jj) * MTRX_SIZE + (ii + i)] += matrA[(jj + j) * MTRX_SIZE + (kk + k)] * matrB[(kk + k)*MTRX_SIZE + (ii + i)];
////}
////
////int main() {
////
////	if (MTRX_SIZE % BLCK_SIZE != 0) {
////		cout << "Error: wrong Size and Block Size!" << endl;
////		system("pause");
////		return 1;
////	}
////
////	int* matrA = new int[MTRX_SIZE*MTRX_SIZE];
////	int* matrB = new int[MTRX_SIZE*MTRX_SIZE];
////	int* matrC = new int[MTRX_SIZE*MTRX_SIZE];
////	int* matrCchecked = new int[MTRX_SIZE*MTRX_SIZE];
////
////	initialization_matrix(matrA, matrB, matrC, matrCchecked);
////
////	TBB_Fox s(matrA, matrB, matrC);
////
////	task_scheduler_init init;
////
////	tick_count start = tick_count::now();
////	parallel_for(blocked_range<int>(0, (MTRX_SIZE / BLCK_SIZE), (MTRX_SIZE / BLCK_SIZE)/2), s);// /4  (MTRX_SIZE / BLCK_SIZE)
////	tick_count finish = tick_count::now();
////	cout << "TBB Time = " << (finish - start).seconds() << " seconds\n";
////	tick_count start1 = tick_count::now();
////	matrix_mul(matrA, matrB, matrCchecked);
////	tick_count finish1 = tick_count::now();
////	init.terminate();
////	cout << "Seq Time = " << (finish1 - start1).seconds() << " seconds\n";
////	cout << "Acceleration = " << ((finish1 - start1).seconds()) / ((finish - start).seconds()) << endl;
////
////	if (MTRX_SIZE < 11) {
////		cout << "Matrix 1" << endl;
////		print_matrix(matrA, MTRX_SIZE, MTRX_SIZE);
////		cout << endl << "Matrix 2" << endl;
////		print_matrix(matrB, MTRX_SIZE, MTRX_SIZE);
////		cout << endl;
////		cout << endl << "Fox" << endl;
////		print_matrix(matrC, MTRX_SIZE, MTRX_SIZE);
////		cout << endl;
////		cout << endl << "Basic" << endl;
////		print_matrix(matrCchecked, MTRX_SIZE, MTRX_SIZE);
////		cout << endl;
////	}
////
////	if (matrix_equal(matrC, matrCchecked)) cout << "Answer is correct\n";
////	else cout << "Answer is incorrect\n";
////
////	free(matrA);
////	free(matrB);
////	free(matrC);
////	free(matrCchecked);
////	system("pause");
////}

#include <tbb/tbb.h>

#include <iostream>
#include <ctime>
#include <random>

double* CreateRandMatrix(int size) {
	double* tmp = new double[size*size];
	for (int i = 0; i < size*size; i++)
		tmp[i] = (std::rand() % 10000) / 1000.0f;
	return tmp;
}
void CheckEqual(double* A, double*B, int n) {
	for (int i = 0; i < n*n; i++) {
		if (std::fabs(A[i] - B[i]) > 0.000001) {
			std::cout << "NOT equal" << std::endl;
			return;
		}
	}
	std::cout << "equal" << std::endl;
}
int PrintMatrix(double* A, int N) {
	for (int i = 0; i < N*N; i += N) {
		for (int j = 0; j < N; j++)
			std::cout << A[i + j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	return 1;
}
void MultiMatrix(double* A, double* B, double* C, int n, int bSize) {
	for (int i = 0; i < bSize; ++i)
		for (int j = 0; j < bSize; ++j)
			for (int k = 0; k < bSize; ++k) {
				C[i * n + j] += A[i * n + k] * B[k * n + j];
			}
}
class Multiplicator {
	double *A, *B, *C;
	int n, q;
public:
	Multiplicator(double *a, double *b, double *c,
		int N, int Q) : A(a), B(b), C(c), n(N), q(Q) {}
	void operator()(const tbb::blocked_range2d<int>& r) const {
		int bSize = n / q;
		for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
			for (int j = r.cols().begin(); j < r.cols().end(); ++j) {
				for (int k = 0; k < q; ++k) {
					MultiMatrix(&A[(i*n + k)*bSize], &B[(k*n + j)*bSize], &C[(i*n + j)*bSize], n, bSize);
				}
			}
		}
	}
};

void Foxs(double* A, double* B, double* C, int n, int q) {
	int bSize = n / q;
	for (int i = 0; i < q; i++)
		for (int j = 0; j < q; j++) {
			for (int k = 0; k < q; k++) {
				MultiMatrix(&A[(i*n + j)*bSize], &B[(j*n + k)*bSize], &C[(i*n + k)*bSize], n, bSize);
			}
		}
}
void FoxsTBB(double* A, double* B, double* C, int n, int q, int grainSize) {
	tbb::task_scheduler_init init(4);
	tbb::parallel_for(tbb::blocked_range2d<int>(0, q, grainSize, 0, q, grainSize),
		Multiplicator(A, B, C, n, q));
}
int main(int argc, char** argv) {
	double *A, *B, *C, *C1;
	int N = 1024, q = 256;
	if (argc == 3) {
		N = atoi(argv[1]);
		q = atoi(argv[2]);
	}
	srand((unsigned int)time(0));
	A = CreateRandMatrix(N);
	B = CreateRandMatrix(N);
	//  PrintMatrix(A, N);
	C = new double[N*N];
	C1 = new double[N*N];
	for (int i = 0; i < N*N; i++) {
		C[i] = 0.0;
		C1[i] = 0.0;
	}
	int grainSize = 1;
	tbb::tick_count t1 = tbb::tick_count::now();
	Foxs(A, B, C, N, q);
	tbb::tick_count t2 = tbb::tick_count::now();
	std::cout << "Time consistent with Fox algorithm is " << (t2 - t1).seconds() << std::endl;
	//  PrintMatrix(C, N);
	tbb::tick_count t1tbb = tbb::tick_count::now();
	FoxsTBB(A, B, C1, N, q, grainSize);
	tbb::tick_count t2tbb = tbb::tick_count::now();
	std::cout << "Time parallel is " << (t2tbb - t1tbb).seconds() << std::endl;
	//  PrintMatrix(C1, N);
	CheckEqual(C, C1, N);

	std::cout << "Acceleration - " << ((t2 - t1).seconds()) / ((t2tbb - t1tbb).seconds()) << std::endl;

	delete[] A;
	delete[] B;
	delete[] C;
	delete[] C1;
	system("pause");
	return 0;
}