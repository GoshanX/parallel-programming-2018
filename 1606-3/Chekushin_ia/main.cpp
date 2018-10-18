#include "mpi.h" 
#include <time.h>
#include <math.h>
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;
int main(int argc, char *argv[])

{
	int ProcRank, ProcNum;
	int n,k;
	double avg;
	int proc_sum = 0;
	int sumpar = 0;
	double tstart, tfinish, t1, t2, t3;
	double sumlin = 0, avglin;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status Status;
	int *a,*tmp;
	n = stoi(argv[1]);
	a = new int[n];
	for (int i = 0; i < n; i++)
	{
		a[i] = rand() % 10;
	}
	k = n / ProcNum;
	tmp = new int[k];

	if (ProcRank == 0)
	{
		
		tstart = MPI_Wtime();
		for (int i = 1; i < ProcNum; i++)
		{
			tmp = a + k*(i - 1);
			MPI_Send(tmp, k, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		for (int i = (ProcNum - 1)*k; i < n; i++)
			sumpar = sumpar + a[i];

		for (int i = 1; i<ProcNum; i++)
		{
			MPI_Recv(&proc_sum, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Status);
			sumpar += proc_sum;
		}
		avg = (double)sumpar / n;
		tfinish = MPI_Wtime();

		t1 = MPI_Wtime();
		for (int i = 0; i < n; i++)
			sumlin += a[i];
		avglin = (double)sumlin / n;
		t2 = MPI_Wtime();
		printf("\n Parallel prog");
		printf("\n\n AVG SumPar = ");
		printf("%f", avg);
		printf("\n\n Time = %10.10f", tfinish - tstart);
		printf("\n\n Sequential prog");
		printf("\n\n AVG SumLin = ");
		printf("%f", avglin);
	    printf("\n\n Time = %10.10f", t2 - t1);
		printf("\n");
		printf("Acceleration = %f", (t2 - t1) / (tfinish - tstart));

	}
	if (ProcRank != 0)
	{
		MPI_Recv(tmp, k, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
		for (int i = 0; i<k; i++)
			proc_sum += tmp[i];
		MPI_Send(&proc_sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}