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