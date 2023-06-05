/*MAIN FLUID SIMULATOR BY USING FLUX METHOD*/
#include <stdio.h>
#include "Parameters.h"
#include "Arrays.h"
#include "ScalarsCons.h"
#include "TimeStep.h"

int main(){
    float t = 0;
    int i,j;
    ScalarArrays Scalars;
    ConservativeArrays Conservatives;

    /*
    Setting initial random conditions as seen in:
    https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
    */
    SetRandomInitialConditions(Ny, Nx, P0, D1, D2, V1, V2, &Scalars);

    while (t < FinalT){
        /* Fill conservative arrays with their values dependent of Scalar arrays */
        ObtainConservativeValues(&Conservatives, &Scalars);

        /* Find optimal timestep to mantain stability in the simulation */
        GetOptimalTimeStep(&Scalars, MinDiff);
    }

    return 0;
}