#ifndef SCALARSCONSEQS_H
#define SCALARSCONSEQS_H
#include "Arrays.h"

extern void ObtainConservativeValues(ConservativeArrays *Conservatives, ScalarArrays *Scalars);
extern void ObtainScalarsValues(ScalarArrays *Scalars, ConservativeArrays *Conservatives);

#endif