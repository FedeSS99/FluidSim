/*MAIN FLUID SIMULATOR BY USING FLUX METHOD*/
#include "Arrays.h"
#include <stdio.h>

int main(){
    double t = 0;
    double dt, C_2, AveYKinetic;

    /* Primitives and Conservative arrays*/
    Primitives Scalars;
    Primitives PrimeScalars;
    Conservatives Conservatives;

    /* Gradients for each Primitive*/
    Grad dDens;
    Grad dVx;
    Grad dVy;
    Grad dPres;

    /* MidStep derivatives in space for each Primitive*/
    MidSpaceArr MidDens;
    MidSpaceArr MidVx;
    MidSpaceArr MidVy;
    MidSpaceArr MidPres;
    /* MidStep derivatives in space for Energy*/
    MidSpaceArr MidEne;

    /* Maximum signal speed array */
    MaxSpeed MaxV;

    /* Fluxes arrays */
    FluxArr Fluxes;


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
        GetOptimalDiffusiveTerm(&MaxV, &MidDens, &MidPres, &MidVx, &MidVy);

        /* Add diffusive terms */
        AddDiffusiveTerms(&Fluxes, &MaxV, &MidDens, &MidEne, &MidVx, &MidVy);

        /* Apply fluxes to conservative terms*/
        AddFluxesToConservatives(dt, &Conservatives, &Fluxes);
        
        /* Compute Primitive values from Conservative values */
        ObtainScalarValues(&Scalars, &Conservatives);

        /* Compute average kinetic energy over Y coordinate */
        AveYKinetic = AverageKineticEnergy_Y_OverDomain(&Scalars);
        printf("At t=%.3f with dt=%.3e with K.E. of Vy=%.3e\n", t, dt, AveYKinetic);

        t += dt;
    }

    /*Free memory after simulation has finished*/
    free(Scalars.Dens);
    free(Scalars.Pres);
    free(Scalars.Vx);
    free(Scalars.Vy);

    free(PrimeScalars.Dens);
    free(PrimeScalars.Pres);
    free(PrimeScalars.Vx);
    free(PrimeScalars.Vy);

    free(Conservatives.E);
    free(Conservatives.Mass);
    free(Conservatives.Mx);
    free(Conservatives.My);

    free(dDens.DX);
    free(dDens.DY);
    free(dVx.DX);
    free(dVx.DY);
    free(dVy.DX);
    free(dVy.DY);
    free(dPres.DX);
    free(dPres.DY);

    free(MidDens.XL);
    free(MidDens.XR);
    free(MidDens.YB);
    free(MidDens.YT);
    free(MidVx.XL);
    free(MidVx.XR);
    free(MidVx.YB);
    free(MidVx.YT);
    free(MidVy.XL);
    free(MidVy.XR);
    free(MidVy.YB);
    free(MidVy.YT);
    free(MidPres.XL);
    free(MidPres.XR);
    free(MidPres.YB);
    free(MidPres.YT);
    free(MidEne.XL);
    free(MidEne.XR);
    free(MidEne.YB);
    free(MidEne.YT);

    free(MaxV.C2);

    /* Fluxes arrays */
    free(Fluxes.F_DensX);    
    free(Fluxes.F_DensY);    
    free(Fluxes.F_MomxX);    
    free(Fluxes.F_MomxY);    
    free(Fluxes.F_MomyX);    
    free(Fluxes.F_MomyY);    
    free(Fluxes.F_EneX);    
    free(Fluxes.F_EneY);    

    return 0;
}   