#ifndef ARRAYS_H
#define ARRAYS_H
#include "Parameters.h"

const int NyNx = Ny*Nx;

/* Scalar arrays */
struct ScalarArrays{
    float Vx[NyNx];
    float Vy[NyNx];
    float Dens[NyNx];
    float Pres[NyNx];
};

/* Conservative arrays */
struct ConservativeArrays{
    float Mx[NyNx];
    float My[NyNx];
    float Mass[NyNx];
    float E[NyNx];
};

extern void Init2DArrayAtZero(int Ny, int Nx, float A[]);
extern void SetRandomInitialConditions(int Ny, int Nx, float P0, float D1, float D2, float V1, float V2, ScalarArrays *Scalars);

#endif