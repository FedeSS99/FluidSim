#include "Parameters.h"

int GetIndex(int x, int y, int N){
    return  y*N + x;
}

void GetPeriodicIndex(int *index, int N){
    (*index) = (*index) % N;
}