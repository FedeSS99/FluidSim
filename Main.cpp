/*MAIN FLUID SIMULATOR BY USING FLUX METHOD*/
#include <stdio.h>
#include "Parameters.h"
#include "Arrays.h"
#include "ScalarsCons.h"
#include "TimeStep.h"
#include "MathOps.h"
#include "Extrapolate.h"

int main(){
    float t = 0;
    float MinDiff, dt;

    /* Primitives and Conservative arrays*/
    ScalarArrays Scalars;
    ScalarArrays PrimeScalars;
    ConservativeArrays Conservatives;

    /* Gradients for each Primitive*/
    Gradient dDens;
    Gradient dVx;
    Gradient dVy;
    Gradient dPres;

    /* MidStep derivatives in space for each Primitive*/
    MidSpaceStepArrays MidDens;
    MidSpaceStepArrays MidVx;
    MidSpaceStepArrays MidVy;
    MidSpaceStepArrays MidPres;

    /*
    Setting initial random conditions as seen in:
    https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
    */
    SetRandomInitialConditions(&Scalars);

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
        ObtainGradient(&dDens, Scalars.Dens);
        ObtainGradient(&dVx, Scalars.Vx);
        ObtainGradient(&dVy, Scalars.Vy);
        ObtainGradient(&dPres, Scalars.Pres);

        /* Extrapolate values to midsteps in time */
        ComputeMidTimeStep(dt, &PrimeScalars, &Scalars, &dDens, &dVx, &dVy, &dPres);

        /* Extrapolate values to midsteps in space*/
        ComputeMidSpaceStep(&MidDens, PrimeScalars.Dens, &dDens);
        ComputeMidSpaceStep(&MidVx, PrimeScalars.Vx, &dVx);
        ComputeMidSpaceStep(&MidVy, PrimeScalars.Vy, &dVy);
        ComputeMidSpaceStep(&MidPres, PrimeScalars.Pres, &dPres);

        t += dt;
        printf("At t=%.3f with dt=%.3e\n", t, dt);   
    }
    return 0;
}