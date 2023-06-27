#include "Parameters.h"
#include "Arrays.h"
#include <math.h> 

double AverageKineticEnergy_Y_OverDomain(Primitives *Scalars){
    int i,j,index;
    double AverageKE_Y = 0.0;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);
            AverageKE_Y += 0.5*Scalars->Dens[index]*pow(Scalars->Vy[index], 2.0);
        }
    }
    AverageKE_Y /= NyNx;
    return AverageKE_Y;
}