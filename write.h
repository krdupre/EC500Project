#include <iostream>
#include <fstream>

using namespace std;

void writet(double** phi, double time)
{
    int i,j;

    ofstream file;
    file.open("data." + to_string ((int) time));
    file.precision(6);
    for (i = 0; i < N; i++)
    {
	for (j = 0; j < N; j++)
	{
	    file << i << "\t" << j << "\t" << phi[i][j] << "\n";
	}
    }
    file.close();
}
