#include <cmath>

#ifndef g
#define g 9.81
#endif

#ifndef h
#define h 1.0
#endif

void saint_venant(double** u_old, double** u_new, double** eta1, double** eta2, double** v_old, double** v_new,double dt, int start, int end)
{
    int i,j;

    for (i = start; i < end; i++)
    {
	for (j = 0; j < N; j++)
	{
 	    u_new[i][j] = (1/eta2[i][j])*(eta1[i][j]*u_old[i][j] 
		- dt/(2*h)*(eta1[(i+1)%N][j]*u_old[(i+1)%N][j]*u_old[(i+1)%N][j] 
		+ g/2*eta1[(i+1)%N][j]*eta1[(i+1)%N][j]
	        - eta1[(N+i-1)%N][j]*u_old[(N+i-1)%N][j]*u_old[(N+i-1)%N][j] 
		- g/2*eta1[(N+i-1)%N][j]*eta1[(N+i-1)%N][j])
		- dt/(2*h)*(eta1[i][(j+1)%N]*u_old[i][(j+1)%N]*v_old[i][(j+1)%N]
		- eta1[i][(N+j-1)%N]*u_old[i][(N+j-1)%N]*v_old[i][(N+j-1)%N]));

	    v_new[i][j] = (1/eta2[i][j])*(eta1[i][j]*v_old[i][j] 
		- dt/(2*h)*(eta1[i][(j+1)%N]*v_old[i][(j+1)%N]*v_old[i][(j+1)%N] 
		+ g/2*eta1[i][(j+1)%N]*eta1[i][(j+1)%N]
	        - eta1[i][(N+j-1)%N]*v_old[i][(N+j-1)%N]*v_old[i][(N+j-1)%N] 
		- g/2*eta1[i][(N+j-1)%N]*eta1[i][(N+j-1)%N])
		- dt/(2*h)*(eta1[(i+1)%N][j]*v_old[(i+1)%N][j]*u_old[(i+1)%N][j]
		- eta1[(N+i-1)%N][j]*v_old[(N+i-1)%N][j]*u_old[(N+i-1)%N][j]));
	}
    }
}

void momentum(double** eta_old, double** eta_new, double** u, double** v,  double dt, int start, int end)
{
    int i,j;

    for (i = start; i < end; i++)
    {
	for (j = 0; j < N; j++)
	{
	    eta_new[i][j] = eta_old[i][j] - dt/2*(eta_old[(i+1)%N][j]*u[(i+1)%N][j]
					  -eta_old[(N+i-1)%N][j]*u[(N+i-1)%N][j])
		                          - dt/2*(eta_old[i][(j+1)%N]*v[i][(j+1)%N]
					  -eta_old[i][(N+j-1)%N]*v[i][(N+j-1)%N]);
	}
    }
}
