#include <cmath>

#ifndef g
#define g 9.81
#endif

void saint_venant(double** u_old, double** u_new, double** eta1, double** eta2, double** v_old, double** v_new,double dt)
{
    int i,j;

    for (i = 1; i < N-1; i++)
    {
	for (j = 1; j < N-1; j++)
	{
 	    u_new[i][j] = (1/eta2[i][j])*(eta1[i][j]*u_old[i][j] 
		        - dt/2*(eta1[i+1][j]*u_old[i+1][j]*u_old[i+1][j] + g/2*eta1[i+1][j]*eta1[i+1][j]
	                      - eta1[i-1][j]*u_old[i-1][j]*u_old[i-1][j] - g/2*eta1[i-1][j]*eta1[i-1][j])
			- dt/2*(eta1[i][j+1]*u_old[i][j+1]*v_old[i][j+1]
			      - eta1[i][j-1]*u_old[i][j-1]*v_old[i][j-1]));

	    v_new[i][j] = (1/eta2[i][j])*(eta1[i][j]*v_old[i][j] 
		        - dt/2*(eta1[i][j+1]*v_old[i][j+1]*v_old[i][j+1] + g/2*eta1[i][j+1]*eta1[i][j+1]
	                      - eta1[i][j-1]*v_old[i][j-1]*v_old[i][j-1] - g/2*eta1[i][j-1]*eta1[i][j-1])
			- dt/2*(eta1[i+1][j]*v_old[i+1][j]*u_old[i+1][j]
			      - eta1[i-1][j]*v_old[i-1][j]*u_old[i-1][j]));
	}
    }

    u_new[0][0] = (1/eta2[0][0])*(eta1[0][0]*u_old[0][0] 
	     - dt/2*(eta1[1][0]*u_old[1][0]*u_old[1][0] + g/2*eta1[1][0]*eta1[1][0]
		   - eta1[N-1][0]*u_old[N-1][0]*u_old[N-1][0] - g/2*eta1[N-1][0]*eta1[N-1][0])
	     - dt/2*(eta1[0][1]*u_old[0][1]*v_old[0][1]
		   - eta1[0][N-1]*u_old[0][N-1]*v_old[0][N-1]));

    u_new[0][N-1] = (1/eta2[0][N-1])*(eta1[0][N-1]*u_old[0][N-1] 
		- dt/2*(eta1[1][N-1]*u_old[1][N-1]*u_old[1][N-1] + g/2*eta1[1][N-1]*eta1[1][N-1]
	              - eta1[N-1][N-1]*u_old[N-1][N-1]*u_old[N-1][N-1] - g/2*eta1[N-1][N-1]*eta1[N-1][N-1])
	 	- dt/2*(eta1[0][0]*u_old[0][0]*v_old[0][0]
		      - eta1[0][N-2]*u_old[0][N-2]*v_old[0][N-2]));

    u_new[N-1][0] = (1/eta2[N-1][0])*(eta1[N-1][0]*u_old[N-1][0] 
		  - dt/2*(eta1[0][0]*u_old[0][0]*u_old[0][0] + g/2*eta1[0][0]*eta1[0][0]
	                - eta1[N-2][0]*u_old[N-2][0]*u_old[N-2][0] - g/2*eta1[N-2][0]*eta1[N-2][0])
		  - dt/2*(eta1[N-1][1]*u_old[N-1][1]*v_old[N-1][1]
			- eta1[N-1][N-1]*u_old[N-1][N-1]*v_old[N-1][N-1]));

    u_new[N-1][N-1] = (1/eta2[N-1][N-1])*(eta1[N-1][N-1]*u_old[N-1][N-1] 
	       	    - dt/2*(eta1[0][N-1]*u_old[0][N-1]*u_old[0][N-1] + g/2*eta1[0][N-1]*eta1[0][N-1]
		     	  - eta1[N-2][N-1]*u_old[N-2][N-1]*u_old[N-2][N-1] - g/2*eta1[N-2][N-1]*eta1[N-2][N-1])
		    -dt/2*(eta1[N-1][0]*u_old[N-1][0]*v_old[N-1][0]
			 - eta1[N-1][N-2]*u_old[N-1][N-2]*v_old[N-1][N-2]));

    v_new[0][0] = (1/eta2[0][0])*(eta1[0][0]*v_old[0][0] 
		- dt/2*(eta1[0][1]*v_old[0][1]*v_old[0][1] + g/2*eta1[0][1]*eta1[0][1]
	              - eta1[0][N-1]*v_old[0][N-1]*v_old[0][N-1] - g/2*eta1[0][N-1]*eta1[0][N-1])
		- dt/2*(eta1[1][0]*v_old[1][0]*u_old[1][0]
		      - eta1[N-1][0]*v_old[N-1][0]*u_old[N-1][0]));

    v_new[0][N-1] = (1/eta2[0][N-1])*(eta1[0][N-1]*v_old[0][N-1] 
		        - dt/2*(eta1[0][0]*v_old[0][0]*v_old[0][0] + g/2*eta1[0][0]*eta1[0][0]
	                      - eta1[0][N-2]*v_old[0][N-2]*v_old[0][N-2] - g/2*eta1[0][N-2]*eta1[0][N-2])
			- dt/2*(eta1[1][N-1]*v_old[1][N-1]*u_old[1][N-1]
			      - eta1[N-1][N-1]*v_old[N-1][N-1]*u_old[N-1][N-1]));

    v_new[N-1][0] = (1/eta2[N-1][0])*(eta1[N-1][0]*v_old[N-1][0] 
		        - dt/2*(eta1[N-1][1]*v_old[N-1][1]*v_old[N-1][1] + g/2*eta1[N-1][1]*eta1[N-1][1]
	                     - eta1[N-1][N-1]*v_old[N-1][N-1]*v_old[N-1][N-1] - g/2*eta1[N-1][N-1]*eta1[N-1][N-1])
			- dt/2*(eta1[0][0]*v_old[0][0]*u_old[0][0]
			      - eta1[N-2][0]*v_old[N-2][0]*u_old[N-2][0]));

    v_new[N-1][N-1] = (1/eta2[N-1][N-1])*(eta1[N-1][N-1]*v_old[N-1][N-1] 
		    - dt/2*(eta1[N-1][0]*v_old[N-1][0]*v_old[N-1][0] + g/2*eta1[N-1][0]*eta1[N-1][0]
	                  - eta1[N-1][N-2]*v_old[N-1][N-2]*v_old[N-1][N-2] - g/2*eta1[N-1][N-2]*eta1[N-1][N-2])
		    - dt/2*(eta1[0][N-1]*v_old[0][N-1]*u_old[0][N-1]
			  - eta1[N-2][N-1]*v_old[N-2][N-1]*u_old[N-2][N-1]));
}

void momentum(double** eta_old, double** eta_new, double** u, double** v,  double dt)
{
    int i,j;

    for (i = 1; i < N-1; i++)
    {
	for (j = 1; j < N-1; j++)
	{
	    eta_new[i][j] = eta_old[i][j] - dt/2*(eta_old[i+1][j]*u[i+1][j]-eta_old[i-1][j]*u[i-1][j])
		                          - dt/2*(eta_old[i][j+1]*v[i][j+1]-eta_old[i][j-1]*v[i][j-1]);
	}
    }

    eta_new[0][0] = eta_old[0][0] - dt/2*(eta_old[1][0]*u[1][0]-eta_old[N-1][0]*u[N-1][0])
	                          - dt/2*(eta_old[0][1]*v[0][1]-eta_old[0][N-1]*v[0][N-1]);

    eta_new[0][N-1] = eta_old[0][N-1] - dt/2*(eta_old[1][N-1]*u[1][N-1]-eta_old[N-1][N-1]*u[N-1][N-1])
	        		      - dt/2*(eta_old[0][0]*v[0][0]-eta_old[0][N-2]*v[0][N-2]);

    eta_new[N-1][0] = eta_old[N-1][0] - dt/2*(eta_old[0][0]*u[0][0]-eta_old[N-2][0]*u[N-2][0])
	         		      - dt/2*(eta_old[N-1][1]*v[N-1][1]-eta_old[N-1][N-1]*v[N-1][N-1]);

    eta_new[N-1][N-1] = eta_old[N-1][N-1] - dt/2*(eta_old[0][N-1]*u[0][N-1]-eta_old[N-2][N-1]*u[N-2][N-1])
	                                  - dt/2*(eta_old[N-1][0]*v[N-1][0]-eta_old[N-1][N-2]*v[N-1][N-2]);
}
