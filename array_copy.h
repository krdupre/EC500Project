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
