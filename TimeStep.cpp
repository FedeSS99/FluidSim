#include "Parameters.h"
#include "Arrays.h"
#include "Index.h"
#include <math.h>

float TimeStepEq(ScalarArrays *Scalars, float Diff, int index){
    float Pres_ind = Scalars->Pres[index];
    float Dens_ind = Scalars->Dens[index];
    float Vx2_ind = powf(Scalars->Vx[index], 2.0);
    float Vy2_ind = powf(Scalars->Vy[index], 2.0);

    return 0.5*(Diff/( sqrtf(gas_c*Pres_ind/Dens_ind) + sqrtf(Vx2_ind + Vy2_ind)));
}

void GetOptimalTimeStep(float *dt, ScalarArrays *Scalars, float Diff){
    int i,j, index;
    float dt_ref = TimeStepEq(Scalars, Diff, 0);
    float dt_index;
    
    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);

            dt_index = TimeStepEq(Scalars, Diff, index);
            if (dt_index < dt_ref){
                dt_ref = dt_index;
                (*dt) = dt_ref;
            } 
        }
    }
}