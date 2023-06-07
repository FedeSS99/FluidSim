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

struct Gradient{
    float DX[NyNx];
    float DY[NyNx];
};

struct MidSpaceStepArrays{
    float XL[NyNx];
    float XR[NyNx];
    float YB[NyNx];
    float YT[NyNx];
};

extern void SetRandomInitialConditions(ScalarArrays *Scalars);

#endif