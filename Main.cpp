/*MAIN FLUID SIMULATOR BY USING FLUX METHOD*/
#include <stdio.h>
#include "Parameters.h"
#include "Arrays.h"
#include "ScalarsCons.h"
#include "TimeStep.h"
#include "MathOps.h"

int main(){
    float t = 0;
    float MinDiff, dt;
    int i,j;

    /* Primitives and Conservative arrays*/
    ScalarArrays Scalars;
    ScalarArrays PrimeScalars;
    ConservativeArrays Conservatives;

    /* Gradients for each Primitive*/
    Gradient dDens;
    Gradient dVx;
    Gradient dVy;
    Gradient dPres;

    /*
    Setting initial random conditions as seen in:
    https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
    */
    SetRandomInitialConditions(Ny, Nx, P0, D1, D2, V1, V2, &Scalars);

    /* Find min difference between dx and dy*/
    if (dy > dx){
        MinDiff = dx;
    }
    else{
        MinDiff = dy;
    }

    /* Start running simulation until FinalT is reached */
    while (t < FinalT){

        /* Fill conservative arrays with their values dependent of Scalar arrays */
        ObtainConservativeValues(&Conservatives, &Scalars);

        /* Find optimal timestep to mantain stability in the simulation */
        GetOptimalTimeStep(&dt, &Scalars, MinDiff);

        /* Obtain gradients for each primitive */
        dDens = ObtainGradient(Scalars.Dens);
        dVx = ObtainGradient(Scalars.Vx);
        dVy = ObtainGradient(Scalars.Vy);
        dPres = ObtainGradient(Scalars.Pres);

        /* Extrapolate values to midsteps in time */
        

        t += dt;
        printf("At t=%.3f with dt=%.3e\n", t, dt);

        
    }

    return 0;
}