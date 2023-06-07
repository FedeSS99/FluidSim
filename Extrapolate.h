#ifndef EXTRAPOLATIONS_H
#define EXTRAPOLATIONS_H
#include "Arrays.h"

extern void ComputeMidTimeStep(float dt, ScalarArrays *PrimeScalars, ScalarArrays *Scalars, Gradient *dDens, Gradient *dVx, Gradient *dVy, Gradient *dPres);
extern void ComputeMidSpaceStep(MidSpaceStepArrays *MidSpaceArr, float A[], Gradient *dA);

#endif