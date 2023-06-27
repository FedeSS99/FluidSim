#include "Parameters.h"
#include "Arrays.h"
#include <math.h>

double DiffusiveTerm(double Dens, double Pres, double V){
    return sqrt(gas_c * Pres / Dens) + fabs(V);
}

void GetOptimalDiffusiveTerm(MaxSpeed *MaxSigV, MidSpaceArr *MidDens, MidSpaceArr *MidPres, MidSpaceArr *MidVx, MidSpaceArr *MidVy){
    double C_2;
    int i,j,k,index;
    double C_terms[4] = {0.0, 0.0, 0.0, 0.0};

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);

            C_terms[0] = DiffusiveTerm(MidDens->XL[index],MidPres->XL[index],MidVx->XL[index]);
            C_terms[1] = DiffusiveTerm(MidDens->XR[index],MidPres->XR[index],MidVx->XR[index]);
            C_terms[2] = DiffusiveTerm(MidDens->YB[index],MidPres->YB[index],MidVy->YB[index]);
            C_terms[3] = DiffusiveTerm(MidDens->YT[index],MidPres->YT[index],MidVy->YT[index]);

            C_2 = C_terms[0];
            for (k = 1; k < 4; k++){
                if (C_terms[k] > C_2){
                    C_2 = C_terms[k];
                }
            }
            MaxSigV->C2[index] = C_2;
        }
    }
}

void AddDiffusiveTerms(FluxArr *Fluxes, MaxSpeed *MaxSigV, MidSpaceArr *MidDens, MidSpaceArr *MidEne, MidSpaceArr *MidVx, MidSpaceArr *MidVy){
    int i,j,index;
    double C_2;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j,i,Nx);

            C_2 = MaxSigV->C2[index];
            
            Fluxes->F_DensX[index] -= C_2*(MidDens->XL[index] - MidDens->XR[index]);
            Fluxes->F_DensY[index] -= C_2*(MidDens->YB[index] - MidDens->YT[index]);

            Fluxes->F_MomxX[index] -= C_2*(MidDens->XL[index] * MidVx->XL[index] - MidDens->XR[index] * MidVx->XR[index]);
            Fluxes->F_MomxY[index] -= C_2*(MidDens->YB[index] * MidVx->YB[index] - MidDens->YT[index] * MidVx->YT[index]);

            Fluxes->F_MomyX[index] -= C_2*(MidDens->XL[index] * MidVy->XL[index] - MidDens->XR[index] * MidVy->XR[index]);
            Fluxes->F_MomyY[index] -= C_2*(MidDens->YB[index] * MidVy->YB[index] - MidDens->YT[index] * MidVy->YT[index]);

            Fluxes->F_EneX[index] -= C_2*(MidEne->XL[index] - MidEne->XR[index]);
            Fluxes->F_EneY[index] -= C_2*(MidEne->YB[index] - MidEne->YT[index]);
        }
    }
}