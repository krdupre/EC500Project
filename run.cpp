#define g 9.81
#define N 256

#include "saint_venant.h"
#include "array_copy.h"
#include "write.h"

int main(int argv, char** argc)
{
    int i,j;

    double time = 1.0;
    double dt = 0.001;
    int N_write = 100;

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
	    eta_o[i][j] = 1.0;
	}
    }
    eta_o[N/2][N/2] += 0.5;
    for (i = 0; i < N; i++)
    {
	for (j = 0; j < N; j++)
	{
	    u_n[i][j] = 0.0;
	    v_n[i][j] = 0.0;
	}
    }

    int Nt = (int) time/dt;

    //writet(eta_o,0.0);

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

    return 0;
}
