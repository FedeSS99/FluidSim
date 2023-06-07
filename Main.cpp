/*MAIN FLUID SIMULATOR BY USING FLUX METHOD*/
#include <stdio.h>
#include "Parameters.h"
#include "Arrays.h"
#include "ScalarsCons.h"
#include "TimeStep.h"
#include "MathOps.h"
#include "Extrapolate.h"
#include "Fluxes.h"
#include "Diffusion.h"

int main(){
    float t = 0;
    float dt, C_2, AveYKinetic;

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
    /* MidStep derivatives in space for Energy*/
    MidSpaceStepArrays MidEne;

    /* Maximum signal speed array */
    MaxSignalSpeed MaxSigV;

    /* Fluxes arrays */
    FluxesArrays Fluxes;
    

    /*
    Setting initial random conditions as seen in:
    https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
    */
    SetRandomInitialConditions(&Scalars);

    /* Fill conservative arrays with their values dependent of Scalar arrays */
    ObtainConservativeValues(&Conservatives, &Scalars);

    /* Start running simulation until FinalT is reached */
    while (t < FinalT){
        /* Find optimal timestep to mantain stability in the simulation */
        GetOptimalTimeStep(&dt, &Scalars);

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
        ComputeMidSpaceStepForEnergy(&MidEne, &MidDens, &MidPres, &MidVx, &MidVy);

        /* Compute fluxes for Density, Momentum and Energy */
        ComputeFluxes(&Fluxes, &MidDens, &MidPres, &MidVx, &MidVy, &MidEne);

        /* Obtain optimal diffusion coefficient */
        GetOptimalDiffusiveTerm(&MaxSigV, &MidDens, &MidPres, &MidVx, &MidVy);

        /* Add diffusive terms */
        AddDiffusiveTerms(&Fluxes, &MaxSigV, &MidDens, &MidEne, &MidVx, &MidVy);

        /* Apply fluxes to conservative terms*/
        AddFluxesToConservatives(dt, &Conservatives, &Fluxes);
        
        /* Compute Primitive values from Conservative values */
        ObtainScalarValues(&Scalars, &Conservatives);

        /* Compute average kinetic energy over Y coordinate */
        AveYKinetic = AverageKineticEnergy_Y_OverDomain(&Scalars);
        printf("At t=%.3f with dt=%.3e with K.E. of Vy=%.3e\n", t, dt, AveYKinetic);

        t += dt;
    }
    return 0;
}   