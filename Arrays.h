#ifndef ARRAYS_H
#define ARRAYS_H
#include "Parameters.h"

#include <stdlib.h>

const int NyNx = Ny*Nx;

/* Scalar arrays */
typedef struct ScalarArrays{
    double *Vx = (double *)calloc(NyNx, sizeof(double));
    double *Vy = (double *)calloc(NyNx, sizeof(double));
    double *Dens = (double *)calloc(NyNx, sizeof(double));
    double *Pres = (double *)calloc(NyNx, sizeof(double));
} Primitives;

/* Conservative arrays */
typedef struct ConservativeArrays{
    double *Mx = (double *)calloc(NyNx, sizeof(double));
    double *My = (double *)calloc(NyNx, sizeof(double));
    double *Mass = (double *)calloc(NyNx, sizeof(double));
    double *E = (double *)calloc(NyNx, sizeof(double));
} Conservatives;

/* Gradient arrays */
typedef struct Gradient{
    double *DX = (double *)calloc(NyNx, sizeof(double));
    double *DY = (double *)calloc(NyNx, sizeof(double));
} Grad;

/* Middle Space extrapolation arrays */
typedef struct MidSpaceStepArrays{
    double *XL = (double *)calloc(NyNx, sizeof(double));
    double *XR = (double *)calloc(NyNx, sizeof(double));
    double *YB = (double *)calloc(NyNx, sizeof(double));
    double *YT = (double *)calloc(NyNx, sizeof(double));
} MidSpaceArr;

/* Flux of Density, Momentum and Energy */
typedef struct FluxesArrays{
    double *F_DensX = (double *)calloc(NyNx, sizeof(double));
    double *F_DensY = (double *)calloc(NyNx, sizeof(double));
    double *F_MomxX = (double *)calloc(NyNx, sizeof(double));
    double *F_MomxY = (double *)calloc(NyNx, sizeof(double));
    double *F_MomyX = (double *)calloc(NyNx, sizeof(double));
    double *F_MomyY = (double *)calloc(NyNx, sizeof(double));
    double *F_EneX = (double *)calloc(NyNx, sizeof(double));
    double *F_EneY = (double *)calloc(NyNx, sizeof(double));
} FluxArr;

/* Max signal speed array */
typedef struct MaxSignalSpeed{
    double *C2 = (double *)calloc(NyNx, sizeof(double));
} MaxSpeed;

extern int GetIndex(int x, int y, int N);
extern void GetPeriodicIndex(int *index, int N);

extern void ObtainConservativeValues(ConservativeArrays *Conservatives, ScalarArrays *Scalars);
extern void ObtainScalarValues(ScalarArrays *Scalars, ConservativeArrays *Conservatives);

extern void ObtainGradient(Gradient *DA, double A[]);
extern double AverageKineticEnergy_Y_OverDomain(ScalarArrays *Scalars);

extern void SetRandomInitialConditions(ScalarArrays *Scalars);

extern void ComputeMidTimeStep(double dt, ScalarArrays *PrimeScalars, ScalarArrays *Scalars, Gradient *dDens, Gradient *dVx, Gradient *dVy, Gradient *dPres);
extern void ComputeMidSpaceStep(MidSpaceStepArrays *MidSpaceArr, double A[], Gradient *dA);
extern void ComputeMidSpaceStepForEnergy(MidSpaceStepArrays *MidSpaceE, MidSpaceStepArrays *MidSpacePres, MidSpaceStepArrays *MidSpaceDens, MidSpaceStepArrays *MidSpaceVx, MidSpaceStepArrays *MidSpaceVy);

extern void GetOptimalDiffusiveTerm(MaxSignalSpeed *MaxSigV, MidSpaceStepArrays *MidDens, MidSpaceStepArrays *MidPres, MidSpaceStepArrays *MidVx, MidSpaceStepArrays *MidVy);
extern void AddDiffusiveTerms(FluxesArrays *Fluxes, MaxSignalSpeed *MaxSigV, MidSpaceStepArrays *MidDens, MidSpaceStepArrays *MidEne, MidSpaceStepArrays *MidVx, MidSpaceStepArrays *MidVy);

extern void ComputeFluxes(FluxesArrays *Fluxes, MidSpaceStepArrays *MidDens, MidSpaceStepArrays *MidPres, MidSpaceStepArrays *MidVx, MidSpaceStepArrays *MidVy, MidSpaceStepArrays *MidEn);
extern void AddFluxesToConservatives(double dt, ConservativeArrays *Cons, FluxesArrays *Fluxes);

extern void GetOptimalTimeStep(double *dt, ScalarArrays *Scalars);

#endif