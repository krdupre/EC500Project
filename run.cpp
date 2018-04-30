#define g 9.81
#define N 256
#define h 1.0

#include "run.h"

int main(int argv, char** argc)
{
    int i;

    double time = 1.0;
    double dt = 0.0001;
    int N_write = 100;

    int points = 1; //number of initialized points
    double** init = new double* [3];
    for (i = 0; i < 3; i++) { init[i] = new double [points+1]; }
    init[0][0] = points;
    init[2][0] = 1.0; //default height

    init[0][1] = N/2; //x coordinate
    init[1][1] = N/2; //y coordinate
    init[2][1] = 0.5; //deviation from default

    runSerial(time,dt,N_write,init);

    return 0;
}
