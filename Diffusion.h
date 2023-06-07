#ifndef DIFFUSION_H
#define DIFFUSION_H
#include "Arrays.h"

extern void GetOptimalDiffusiveTerm(MaxSignalSpeed *MaxSigV, MidSpaceStepArrays *MidDens, MidSpaceStepArrays *MidPres, MidSpaceStepArrays *MidVx, MidSpaceStepArrays *MidVy);
extern void AddDiffusiveTerms(FluxesArrays *Fluxes, MaxSignalSpeed *MaxSigV, MidSpaceStepArrays *MidDens, MidSpaceStepArrays *MidEne, MidSpaceStepArrays *MidVx, MidSpaceStepArrays *MidVy);

#endif