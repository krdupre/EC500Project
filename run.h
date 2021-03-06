#include <mpi.h>
#include <chrono>

#include "saint_venant.h"
#include "array_copy.h"
#include "write.h"

void runSerial(double time, double dt, int N_write, double** init)
{
    int i,j;

    double** eta_o = new double* [N];
    double** u_o = new double* [N];
    double** eta_n = new double* [N];
    double** u_n = new double* [N];
    double** v_o = new double* [N];
    double** v_n = new double* [N];

    for (i = 0; i < N; i++)
    {
	eta_o[i] = new double [N];
	u_o[i] = new double [N];
	eta_n[i] = new double [N];
	v_o[i] = new double [N];
	u_n[i] = new double [N];
	v_n[i] = new double [N];
	for (j = 0; j < N; j++)
	{
	    eta_o[i][j] = init[2][0];
	}
    }
    for(i = 1; i < (int) init[0][0]+1; i++)    
    {
	eta_o[(int)init[0][i]][(int)init[1][i]] += init[2][i];
    }
    for (i = 0; i < N; i++)
    {
	for (j = 0; j < N; j++)
	{
	    u_n[i][j] = 0.0;
	    v_n[i][j] = 0.0;
	}
    }

    int Nt = (int) time/dt;

    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now(); 

    for (i = 0; i < Nt+1; i++)
    {
	momentum(eta_o,eta_n,u_o,v_o,dt,0,N);
	saint_venant(u_o,u_n,eta_o,eta_n,v_o,v_n,dt,0,N);

	a_copy(eta_n,eta_o);
	a_copy(u_n,u_o);
	a_copy(v_n,v_o);

	if (i%N_write == 0)
	{
	    writet(eta_o, dt*i*10);
	}
    }

    chrono::time_point<chrono::steady_clock>end_time = chrono::steady_clock::now();
    chrono::duration<double>difference_in_time = end_time - begin_time;
    double runTime = difference_in_time.count();

    printf("Run took %f seconds.\n", runTime);
}

void runParallel(double time, double dt, int N_write, int argc, char** argv, double** init)
{
    int i,j;
    int N_proc,PID;

    double max_time;

    double** eta_o = new double* [N];
    double** u_o = new double* [N];
    double** eta_l = new double* [N];
    double** eta_g = new double* [N];
    double** u_n = new double* [N];
    double** v_o = new double* [N];
    double** v_n = new double* [N];

    for (i = 0; i < N; i++)
    {
	eta_o[i] = new double [N];
	u_o[i] = new double [N];
	eta_l[i] = new double [N];
	eta_g[i] = new double [N];
	v_o[i] = new double [N];
	u_n[i] = new double [N];
	v_n[i] = new double [N];
	for (j = 0; j < N; j++)
	{
	    eta_o[i][j] = init[2][0];
	}
    }
    for(i = 1; i < (int) init[0][0]+1; i++)    
    {
	eta_o[(int)init[0][i]][(int)init[1][i]] += init[2][i];
    }
    for (i = 0; i < N; i++)
    {
	for (j = 0; j < N; j++)
	{
	    u_o[i][j] = 0.0;
	    v_o[i][j] = 0.0;
	}
    }

    int Nt = (int) time/dt;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &N_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &PID);

    int M = N/N_proc;
    int n_ex  = N - M*N_proc;
    int m1,m2;
    if (PID < n_ex) 
    { 
        m1 = PID*(M+1);
	m2 = (PID+1)*(M+1);
    }
    else
    {
	m1 = n_ex + PID*M;
	m2 = n_ex + (PID+1)*M;
    }

    MPI_l0(eta_l,m1,m2);
    MPI_l0(u_n,m1,m2);
    MPI_l0(v_n,m1,m2);
    MPI_Barrier(MPI_COMM_WORLD);

    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now(); 

    for (i = 0; i < Nt+1; i++)
    {
	momentum(eta_o,eta_l,u_o,v_o,dt,m1,m2);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_copy(eta_l,eta_g);

	saint_venant(u_o,u_n,eta_o,eta_g,v_o,v_n,dt,m1,m2);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_copy(u_n, u_o);
	MPI_copy(v_n, v_o);

	a_copy(eta_g,eta_o);
    }

    chrono::time_point<chrono::steady_clock>end_time = chrono::steady_clock::now();
    chrono::duration<double>difference_in_time = end_time - begin_time;
    double runTime = difference_in_time.count();

    MPI_Allreduce(&runTime, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (PID == 0)
    {
	printf("Number of processes: %d, Time: %f seconds\n", N_proc, max_time);
    }

    MPI_Finalize();
}
