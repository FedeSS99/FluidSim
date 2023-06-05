#include "Parameters.h"

int GetIndex(int x, int y, int Nx){
    return  y*Nx + x;
}

int GetPeriodicIndex(int *index, int N){
    if ((*index) == N){
        (*index) = 0;
    }
    else if ((*index) == -1){
        (*index) = N-1;
    }
}