#include "Arrays.h"
#include "Parameters.h"
#include <math.h>

void ObtainConservativeValues(Conservatives *Conservatives, Primitives *Scalars){
    int i,j,index;
    double mass, v_x, v_y;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);

            mass = Scalars->Dens[index]*Volume;
            Conservatives->Mass[index] = mass;

            v_x = Scalars->Vx[index];
            v_y = Scalars->Vy[index];
            Conservatives->Mx[index] = mass*v_x;
            Conservatives->My[index] = mass*v_y;

            Conservatives->E[index] = (Scalars->Pres[index]*Volume/(gas_c - 1.0)) + 0.5*mass*(pow(v_x, 2.0) + pow(v_y, 2.0));
        }
    }
}

void ObtainScalarValues(Primitives *Scalars, Conservatives *Conservatives){
    int i,j,index;
    double mass, v_x, v_y;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);

            mass = Conservatives->Mass[index];
            Scalars->Dens[index] = mass/Volume;

            Scalars->Vx[index] = Conservatives->Mx[index] / mass;
            Scalars->Vy[index] = Conservatives->My[index] / mass;

            v_x = Scalars->Vx[index];
            v_y = Scalars->Vy[index];
            Scalars->Pres[index] = (Conservatives->E[index]/Volume - 0.5*Scalars->Dens[index]*(pow(v_x, 2.0) + pow(v_y, 2.0)));
            Scalars->Pres[index] *= (gas_c - 1.0);
        }
    }
}