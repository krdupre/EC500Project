void a_copy(double** phi, double** phi_c)
{
    int i,j;

    for (i = 0; i < N; i++)
    {
	for (j = 0; j < N; j++)
	{
	    phi_c[i][j] = phi[i][j];
	}
    }
}

void MPI_copy(double** phi_l, double** phi_g)
{
    int i,j;

    for (i = 0; i < N; i++)
    {
	for (j = 0; j < N; j++)
	{
	    MPI_Allreduce(&phi_l[i][j], &phi_g[i][j], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}
    }
}

void MPI_l0(double** phi_l, int start, int end)
{
    int i,j;

    for (i = 0; i < start; i++)
    {
	for (j = 0; j < N; j++)
	{
	    phi_l[i][j] = 0.;
	}
    }

    for (i = end; i < N; i++)
    {
	for (j = 0; j < N; j++)
	{
	    phi_l[i][j] = 0.;
	}
    }
}
