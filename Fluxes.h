#ifndef FLUXES_H
#define FLUXES_H
#include "Arrays.h"

extern void ComputeFluxes(FluxesArrays *Fluxes, MidSpaceStepArrays *MidDens, MidSpaceStepArrays *MidPres, MidSpaceStepArrays *MidVx, MidSpaceStepArrays *MidVy, MidSpaceStepArrays *MidEn);

#endif