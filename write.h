#include <iostream>
#include <fstream>

#ifndef h
#define h 1.0
#endif

using namespace std;

void writet(double** phi, double time)
{
    int i,j;
    string name = "data.dat";

    ofstream file;
    if (time < 1e-6)
    {file.open(name);}
    else 
    {file.open(name, ios::out | ios::app);}
    file.precision(6);
    for (i = 0; i < N; i++)
    {
	for (j = 0; j < N; j++)
	{
	    file << phi[i][j] << "\t";
	}
    }
    file << "\n";
    file.close();
}
