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

struct FluxesArrays{
    float F_DensX[NyNx];
    float F_DensY[NyNx];
    float F_MomxX[NyNx];
    float F_MomxY[NyNx];
    float F_MomyX[NyNx];
    float F_MomyY[NyNx];
    float F_EneX[NyNx];
    float F_EneY[NyNx];
};

extern void SetRandomInitialConditions(ScalarArrays *Scalars);

#endif