#ifndef SCALARS_CONS_EQS_H
#define SCALARS_CONS_EQS_H
#include "Arrays.h"

extern void ObtainConservativeValues(ConservativeArrays *Conservatives, ScalarArrays *Scalars);
extern void ObtainScalarValues(ScalarArrays *Scalars, ConservativeArrays *Conservatives);

#endif