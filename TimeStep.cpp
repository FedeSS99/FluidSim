#include "Parameters.h"
#include "Arrays.h"
#include "Index.h"
#include <math.h>

float TimeStepEq(ScalarArrays *Scalars, int index){
    float Pres_ind = Scalars->Pres[index];
    float Dens_ind = Scalars->Dens[index];
    float Vx_ind = Scalars->Vx[index];
    float Vy_ind = Scalars->Vy[index];

    float dt_x = 0.5*(dx/( sqrtf(gas_c*Pres_ind/Dens_ind) + fabs(Vx_ind)));
    float dt_y = 0.5*(dy/( sqrtf(gas_c*Pres_ind/Dens_ind) + fabs(Vx_ind)));

    float minDT = (dt_x < dt_y) ? dt_x : dt_y;

    return minDT;
}

void GetOptimalTimeStep(float *dt, ScalarArrays *Scalars){
    int i,j, index;
    float dt_ref = TimeStepEq(Scalars, 0);
    float dt_index;
    
    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);

            dt_index = TimeStepEq(Scalars, index);
            if (dt_index < dt_ref){
                dt_ref = dt_index;
                (*dt) = dt_ref;
            } 
        }
    }
}