// mpi3.cpp : main project file.


#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <thread>
#include <windows.h>
using namespace std;

int ProcNum, Rank, Grid, GridComm, ColComm, RowComm;
int GridC[2];


// инициализация матриц
void InitializationMatrix(double* MatrixA, double* MatrixB, int Size)
{
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size; j++)
		{
			MatrixA[i*Size + j] = rand() % 9;
			MatrixB[i*Size + j] = rand() % 9;
		}
}


void Scatter(double* Matrix, double* MatrixBlock, int Size, int BSize)
{
	double* MatrixRow = new double[BSize*Size];

	if (GridC[1] == 0)  MPI_Scatter(Matrix, BSize*Size, MPI_DOUBLE, MatrixRow, BSize*Size, MPI_DOUBLE, 0, ColComm);

	for (int i = 0; i < BSize; i++)
		MPI_Scatter(&MatrixRow[i*Size], BSize, MPI_DOUBLE, &(MatrixBlock[i*BSize]), BSize, MPI_DOUBLE, 0, RowComm);

	delete[] MatrixRow;
}

// разделение данных на блоки  для умножения
void Distribution(double* MatrixA, double* MatrixB, double* Block, double* BlockB, int Size, int BSize)
{
	Scatter(MatrixA, Block, Size, BSize);
	Scatter(MatrixB, BlockB, Size, BSize);
}

// создание решетки
void CreateGird()
{
	int DimSize[2];
	int Periodic[2];
	int Subdims[2];


	DimSize[0] = Grid;
	DimSize[1] = Grid;
	Periodic[0] = 0;
	Periodic[1] = 0;

	MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);

	MPI_Cart_coords(GridComm, Rank, 2, GridC);

	Subdims[0] = 0;
	Subdims[1] = 1;

	MPI_Cart_sub(GridComm, Subdims, &RowComm);

	Subdims[0] = 1;
	Subdims[1] = 0;

	MPI_Cart_sub(GridComm, Subdims, &ColComm);
}

void BlockAComm(int i, double *BlockA, double* Block, int BSize)
{
	int tmp = (GridC[0] + i) % Grid;

	if (GridC[1] == tmp)
	{
		for (int j = 0; j < BSize*BSize; j++) BlockA[j] = Block[j];
	}

	MPI_Bcast(BlockA, BSize*BSize, MPI_DOUBLE, tmp, RowComm);
}

void BlockBComm(double *BlockB, int BSize)
{
	MPI_Status Status;
	int Next = GridC[0] + 1;
	if (GridC[0] == Grid - 1) Next = 0;
	int Prev = GridC[0] - 1;
	if (GridC[0] == 0) Prev = Grid - 1;

	MPI_Sendrecv_replace(BlockB, BSize*BSize, MPI_DOUBLE, Next, 0, Prev, 0, ColComm, &Status);
}

// последовательное перемножение матриц

void Mult(double* MatrixA, double* MatrixB, double* MatrixR, int Size)
{
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size; j++)
			for (int k = 0; k < Size; k++)
				MatrixR[i*Size + j] += MatrixA[i*Size + k] * MatrixB[k*Size + j];
}

// метод Фокса

void Method(double* BlockA, double* Block, double* BlockB, double* BlockR, int BSize)
{
	for (int i = 0; i < Grid; i++)
	{
		BlockAComm(i, BlockA, Block, BSize);
		Mult(BlockA, BlockB, BlockR, BSize);
		BlockBComm(BlockB, BSize);
	}
}

// соединение данных

void CollectionResult(double *Result, double* BlockR, int Size, int BSize)
{
	double *Row = new double[Size*BSize];

	for (int i = 0; i < BSize; i++) MPI_Gather(&BlockR[i*BSize], BSize, MPI_DOUBLE, &Row[i*Size], BSize, MPI_DOUBLE, 0, RowComm);

	if (GridC[1] == 0) MPI_Gather(Row, BSize*Size, MPI_DOUBLE, Result, BSize*Size, MPI_DOUBLE, 0, ColComm);
	delete[]Row;
}

void PrintMatrix(double *matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << matrix[i * cols + j] << " ";
		}
		cout << endl;
	}
}

int main(int argc, char *argv[])
{
	double *MatrixA = NULL, *BlockA = NULL; // матрица №1
	double *MatrixB = NULL, *BlockB = NULL; // матрица №2 
	double *Result = NULL, *BlockR = NULL; // результат
	double *Block = NULL;
	int Size, BSize; // размер пространства
	double t1, t2;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

	Grid = sqrt((double)ProcNum);

	if (ProcNum != Grid*Grid) { if (Rank == 0) cout << "Wrong count of Processes! Think about it..." << endl; }
	else
	{
		if (Rank == 0) printf("Matrix multiplication: \n");

		CreateGird();

		Size = std::stoi(std::string(argv[1]));

		//if (Rank == 0)
		//{
		//	do
		//	{
		//		printf("Enter size: ");
		//		scanf_s("%d", &Size);

		//		if (Size%Grid != 0) printf("...\n");
		//	} while (Size%Grid != 0);
		//}

		MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

		BSize = Size / Grid;

		BlockA = new double[BSize*BSize];
		BlockB = new double[BSize*BSize];
		BlockR = new double[BSize*BSize];
		Block = new double[BSize*BSize];

		for (int i = 0; i < BSize*BSize; i++) BlockR[i] = 0;

		if (Rank == 0)
		{
			MatrixA = new double[Size*Size];
			MatrixB = new double[Size*Size];
			Result = new double[Size*Size];
			InitializationMatrix(MatrixA, MatrixB, Size);
			if (Size < 11)
			{
				cout << "Matrix A" << endl;
				PrintMatrix(MatrixA, Size, Size);
				cout << endl << "Matrix B" << endl;
				PrintMatrix(MatrixB, Size, Size);
				cout << endl;
			}
		}
		
		t1 = MPI_Wtime();

		Distribution(MatrixA, MatrixB, Block, BlockB, Size, BSize);

		Method(BlockA, Block, BlockB, BlockR, BSize);

		CollectionResult(Result, BlockR, Size, BSize);

		t2 = MPI_Wtime();
		
		if (Rank == 0)
		{
			if (Size < 11)
			{
				cout << endl <<  "Matrix C" << endl;
				PrintMatrix(Result, Size, Size);
				cout << endl;
			}
			cout << "Parallel method: " << t2 - t1 << endl;
		}
		
		/*
		t1 = MPI_Wtime();
		Mult(MatrixA, MatrixB, Result, Size);
		t2 = MPI_Wtime();
		if (Rank == 0)
		{
			cout << "Sequential method: " << t2 - t1 << endl;
			if (Size < 11)
			{
				cout << "Matrix C" << endl;
				PrintMatrix(Result, Size, Size);
				cout << endl;
			}
		}
		*/
		delete[] MatrixA;
		delete[] MatrixB;
		delete[] Result;
		delete[] BlockA;
		delete[] BlockB;
		delete[] BlockR;
	}

	MPI_Finalize();

	//_getch();//Будет ждать пока пользователь не нажмёт любую клавишу

	return 0;
}