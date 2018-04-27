#define g 9.81
#define N 256
#define h 1.0

#include "run.h"

int main(int argv, char** argc)
{
    double time = 1.0;
    double dt = 0.0001;
    int N_write = 100;

    runSerial(time,dt,N_write);

    return 0;
}
