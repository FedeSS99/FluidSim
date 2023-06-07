#include "Parameters.h"
#include "Arrays.h"
#include "Index.h"
#include <math.h> 

void ComputeFluxes(FluxesArrays *Fluxes, MidSpaceStepArrays *MidDens, MidSpaceStepArrays *MidPres, MidSpaceStepArrays *MidVx, MidSpaceStepArrays *MidVy, MidSpaceStepArrays *MidEn){
    int i,j,index;
    float DensX_Prom, DensY_Prom, MomxX_Prom, MomxY_Prom;
    float MomyX_Prom, MomyY_Prom, EneX_Prom, EneY_Prom;
    float PresX_Prom, PresY_Prom;
    
    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);

            DensX_Prom = 0.5*(MidDens->XL[index] + MidDens->XR[index]);
            DensY_Prom = 0.5*(MidDens->YB[index] + MidDens->YT[index]);
            MomxX_Prom = 0.5*(MidDens->XL[index] * MidVx->XL[index] + MidDens->XR[index] * MidVx->XR[index]);
            MomxY_Prom = 0.5*(MidDens->YB[index] * MidVx->YB[index] + MidDens->YT[index] * MidVx->YT[index]);
            MomyX_Prom = 0.5*(MidDens->XL[index] * MidVy->XL[index] + MidDens->XR[index] * MidVy->XR[index]);
            MomyY_Prom = 0.5*(MidDens->YB[index] * MidVy->YB[index] + MidDens->YT[index] * MidVy->YT[index]);
            EneX_Prom = 0.5*(MidEn->XL[index] + MidEn->XR[index]);
            EneY_Prom = 0.5*(MidEn->YB[index] + MidEn->YT[index]);

            PresX_Prom = (gas_c-1)*(EneX_Prom-0.5 * (powf(MomxX_Prom, 2.0) + powf(MomyX_Prom, 2.0)/DensX_Prom));
            PresY_Prom = (gas_c-1)*(EneY_Prom-0.5 * (powf(MomxY_Prom, 2.0) + powf(MomyY_Prom, 2.0)/DensY_Prom));

            Fluxes->F_DensX[index] = MomxX_Prom;
            Fluxes->F_DensY[index] = MomxY_Prom;
            Fluxes->F_MomxX[index] = powf(MomxX_Prom, 2.0)/DensX_Prom + PresX_Prom;
            Fluxes->F_MomxY[index] = MomyY_Prom * MomxY_Prom/DensY_Prom;
            Fluxes->F_MomyX[index] = MomxX_Prom * MomyX_Prom/DensX_Prom;
            Fluxes->F_MomyY[index] = powf(MomyY_Prom, 2.0)/DensY_Prom + PresY_Prom;
            Fluxes->F_EneX[index] = (EneX_Prom+PresX_Prom) * MomxX_Prom/DensX_Prom;
            Fluxes->F_EneY[index] = (EneY_Prom+PresY_Prom) * MomyY_Prom/DensY_Prom;
        }
    }
}

void AddFluxesToConservatives(float dt, ConservativeArrays *Cons, FluxesArrays *Fluxes){
    int i,j,index;
    int i_lower, j_lower, transX_index, transY_index;
    float dtdx = dt*dx;
    float dtdy = dt*dy;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);
            i_lower = i - 1;
            j_lower = j - 1;
            GetPeriodicIndex(&i_lower, Ny);
            GetPeriodicIndex(&j_lower, Nx);
            transX_index = GetIndex(j_lower, i, Nx);
            transY_index = GetIndex(j, i_lower, Nx);

            Cons->Mass[index] += (dtdx*Fluxes->F_DensX[transX_index] + dtdy*Fluxes->F_DensY[transY_index]); 
            Cons->Mx[index] += (dtdx*Fluxes->F_MomxX[transX_index] + dtdy*Fluxes->F_MomxY[transY_index]); 
            Cons->My[index] += (dtdx*Fluxes->F_MomyX[transX_index] + dtdy*Fluxes->F_MomyY[transY_index]); 
            Cons->E[index] += (dtdx*Fluxes->F_EneX[transX_index] + dtdy*Fluxes->F_EneY[transY_index]); 

            Cons->Mass[index] -= (dtdx*Fluxes->F_DensX[index] + dtdy*Fluxes->F_DensY[index]); 
            Cons->Mx[index] -= (dtdx*Fluxes->F_MomxX[index] + dtdy*Fluxes->F_MomxY[index]); 
            Cons->My[index] -= (dtdx*Fluxes->F_MomyX[index] + dtdy*Fluxes->F_MomyY[index]); 
            Cons->E[index] -= (dtdx*Fluxes->F_EneX[index] + dtdy*Fluxes->F_EneY[index]); 
        }
    }
}