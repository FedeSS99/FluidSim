#include "Parameters.h"
#include "Arrays.h"
#include "Index.h"
#include <math.h> 

float AverageKineticEnergy_Y_OverDomain(ScalarArrays *Scalars){
    int i,j,index;
    float AverageKE_Y = 0.0;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);
            AverageKE_Y += 0.5*Scalars->Dens[index]*powf(Scalars->Vy[index], 2.0);
        }
    }
    AverageKE_Y /= NyNx;
    return AverageKE_Y;
}