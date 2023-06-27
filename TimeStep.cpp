#include "Parameters.h"
#include "Arrays.h"
#include <math.h>

double TimeStepEq(Primitives *Scalars, int index){
    double Pres_ind = Scalars->Pres[index];
    double Dens_ind = Scalars->Dens[index];
    double Vx2_ind = pow(Scalars->Vx[index], 2.0);
    double Vy2_ind = pow(Scalars->Vy[index], 2.0);

    double Vgas = sqrt(gas_c*Pres_ind/Dens_ind);
    double Vmag = sqrt(Vx2_ind + Vy2_ind);
    double Vtotal = Vgas + Vmag;

    double dt_x = 0.5*dx/Vtotal;
    double dt_y = 0.5*dy/Vtotal;

    double minDT = (dt_x < dt_y) ? dt_x : dt_y;

    return minDT;
}

void GetOptimalTimeStep(double *dt, Primitives *Scalars){
    int i,j, index;
    double dt_ref = TimeStepEq(Scalars, 0);
    double dt_index;
    
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