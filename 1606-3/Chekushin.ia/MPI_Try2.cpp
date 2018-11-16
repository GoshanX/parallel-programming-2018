

#//#include "stdafx.h"
#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include<iostream>
using namespace std;
int NEW_log2(const int x)
{
	int y=x;
	int i=0;
	while(y%2==0)
	{
		i++;
		y=y/2;
	}
	if((int)(pow(2, (float)(i)))==x)
	{
		return i;
	}
	else
	{
		return 0;
	}
}
int _to(int r, int procNum)
{
	return (r  + procNum) % procNum + 1;
}

int _from(int r, int procNum)
{
	return (r - 1) % procNum;
}

void NEW_SUM( void *a, void *a1, int n, MPI_Datatype t)
{
	int i;
	if(t==MPI_INT)
	{
		for(i=0; i<n; i++)
		{
			((int *)a1)[i]+=((int *)a)[i];
		}
	}
	if(t==MPI_FLOAT)
	{
		for(i=0; i<n; i++)
		{
			((float *)a1)[i]=((float *)a)[i] + ((float *)a1)[i];
		}
	}
	if(t==MPI_DOUBLE)
	{
		for(i=0; i<n; i++)
		{
			((double *)a1)[i]=((double *)a)[i] + ((double *)a1)[i];
		}
	}
}


void NEW_MIN(const void *a, void *a1, int n, MPI_Datatype t)
{
	int i;
	if(t==MPI_INT)
	{
		for(i=0; i<n; i++)
		{
			if(((int *)a1)[i]>((int *)a)[i])
			{
				((int *)a1)[i]=((int *)a)[i];
			}
		}
	}
	if(t==MPI_FLOAT)
	{
		for(i=0; i<n; i++)
		{
			if(((float *)a1)[i]>((float *)a)[i])
			{
				((float *)a1)[i]=((float *)a)[i];
			}
		}
	}
	if(t==MPI_DOUBLE)
	{
		for(i=0; i<n; i++)
		{
			if(((double *)a1)[i]>((double *)a)[i])
			{
				((double *)a1)[i]=((double *)a)[i];
			}
		}
	}
}


int Tree_Allreduce(void* a, void* a1, int n, MPI_Datatype t, MPI_Op o, MPI_Comm c)
{
	int ProcRank, ProcNum, h;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    
	h=(int)(log((float)ProcNum)/log(2.0));
	if(NEW_log2(ProcNum)==0) h++;
	//============================Сбор на 0 процесс============================
	for(int i = 0; i<h;i++)
	{	
			int j=0;
			while(j<ProcNum)
			{
				int k =j+(int)pow(2, (float)(i));
		
				if(k<ProcNum)
				{
					if(ProcRank==j)
					{
						MPI_Status status;
						MPI_Recv(a1,n,t,k,0,MPI_COMM_WORLD,&status);
						if(o==MPI_MIN) NEW_MIN(a,a1,n, t);
						if(o==MPI_SUM) NEW_SUM(a1,a,n,t );
					}
					if(ProcRank==k)
					{
						MPI_Send(a,n,t,j,0,MPI_COMM_WORLD);
					}
				}
				j=j+(int)pow(2, (float)(i+1));
				
			}
	}



	

	//============================Рассылка обратно============================
	int procNum, procRank;
	MPI_Comm_size(c, &procNum);
	MPI_Comm_rank(c, &procRank);

	int _r = _to(procRank, procNum);
	
	if (procRank != 0)
		MPI_Recv(a, n, t, _from(_r / 2, procNum), 0, c, MPI_STATUS_IGNORE);

	if (2 * _r <= procNum)
	{
		MPI_Send(a, n, t, _from(2 * _r, procNum), 0, c);
	}
	
	if (2 * _r + 1 <= procNum)
	{
		MPI_Send(a, n, t, _from(2 * _r + 1, procNum), 0, c);
	}
	for(int i=0; i < n;i++)
			((int*)a1)[i]=((int*)a)[i];
	
	return 0;
}


int main(int argc, char *argv[])
{
	int ProcRank, ProcNum;
	double TStart, TFinish;
	
	MPI_Init(&argc, &argv);	
	
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int m = atoi(argv[0]);
	
	int *a = new int[10];
	int *b = new int[10];
	
	//============================Пример суммирования============================

	for (int i = 0; i < 10; i++)
	{
		a[i] = i;
	}
	for (int i = 0; i < 10; i++)
	{
		b[i] = 0;
	}
	TStart = MPI_Wtime();
	Tree_Allreduce(a,b,10,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	//MPI_Allreduce(a,b,10,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	TFinish = MPI_Wtime();
	for (int k = 0; k < 10; k++)	
		if(ProcRank==m)
		{
			cout<<b[k]<<" ";
			
		}

	//============================Пример минимизирования============================
/*
	int coef=ProcNum-ProcRank-1;
	for (int i = 0; i < 100; i++)
	{
		a[i] = i*coef;
	}
	for (int i = 0; i < 100; i++)
	{
		b[i] = 0;
	}
	
	Tree_Allreduce(a,b,100,ProcRank,ProcNum,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
	
	for (int k = 0; k < 100; k++)	
		if(ProcRank==1)
		{
			cout<<b[k]<<" ";
		}
	
	*/	

	
	if(ProcRank==m) cout<<endl<<"Time of work = "<<TFinish-TStart;
	MPI_Finalize();
	
	
}
