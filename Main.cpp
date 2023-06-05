/*MAIN FLUID SIMULATOR BY USING FLUX METHOD*/
#include <stdio.h>
#include "Parameters.h"
#include "Arrays.h"
#include "ScalarsCons.h" 

int main(){
    int i,j;
    ScalarArrays Scalars;
    ConservativeArrays Conservatives;

    /* Fill scalar arrays with zeros */
    Init2DArrayAtZero(Ny, Nx, Scalars.Vx);
    Init2DArrayAtZero(Ny, Nx, Scalars.Vy);
    Init2DArrayAtZero(Ny, Nx, Scalars.Dens);
    Init2DArrayAtZero(Ny, Nx, Scalars.Pres);

    /*
    Setting initial random conditions as seen in:
    https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
    */
    SetRandomInitialConditions(Ny, Nx, P0, D1, D2, V1, V2, &Scalars);

    /* Fill conservative arrays with their correspondt values */
    ObtainConservativeValues(&Conservatives, &Scalars);

    return 0;
}