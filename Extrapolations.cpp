#include "Parameters.h"
#include "Arrays.h"
#include "Index.h"
#include <math.h> 

void ComputeMidTimeStep(float dt, ScalarArrays *PrimeScalars, ScalarArrays *Scalars, Gradient *dDens, Gradient *dVx, Gradient *dVy, Gradient *dPres){
    int i,j, index;
    float half_dt = 0.5*dt;
    float Dens_ind, Vx_ind, Vy_ind, Pres_ind;
    float DensX_ind, VxX_ind, VyX_ind, PresX_ind;
    float DensY_ind, VxY_ind, VyY_ind, PresY_ind;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);

            Dens_ind = Scalars->Dens[index];
            Vx_ind = Scalars->Vx[index];
            Vy_ind = Scalars->Vy[index];
            Pres_ind = Scalars->Pres[index];

            DensX_ind = dDens->DX[index];
            DensY_ind = dDens->DY[index];
            VxX_ind = dVx->DX[index];
            VxY_ind = dVx->DY[index];
            VyX_ind = dVy->DX[index];
            VyY_ind = dVy->DY[index];
            PresX_ind = dPres->DX[index];
            PresY_ind = dPres->DY[index];

            PrimeScalars->Dens[index] = Dens_ind - half_dt*(Vx_ind*DensX_ind + Vy_ind*DensY_ind + Dens_ind*(VxX_ind + VyY_ind));
            PrimeScalars->Vx[index] = Vx_ind - half_dt*(Vx_ind*VxX_ind + Vy_ind*VxY_ind + PresX_ind/Dens_ind);
            PrimeScalars->Vy[index] = Vy_ind - half_dt*(Vx_ind*VyX_ind + Vy_ind*VyY_ind + PresY_ind/Dens_ind);
            PrimeScalars->Pres[index] = Pres_ind - half_dt*(Vx_ind*PresX_ind + Vy_ind*PresY_ind + gas_c*Pres_ind*(VxX_ind + VyY_ind));
        }
    }
}

void ComputeMidSpaceStep(MidSpaceStepArrays *MidSpaceArr, float A[], Gradient *dA){
    int i,j,index,j_upper,i_upper;
    int transX_index,transY_index;
    float half_stepX = 0.5*dx;
    float half_stepY = 0.5*dy;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);
            
            i_upper = i + 1;
            j_upper = j + 1;
            GetPeriodicIndex(&i_upper, Ny);
            GetPeriodicIndex(&j_upper, Nx);

            transX_index = GetIndex(j_upper, i, Nx);
            transY_index = GetIndex(j, i_upper, Nx);

            MidSpaceArr->XL[index] = A[transX_index] - dA->DX[transX_index]*half_stepX;
            MidSpaceArr->XR[index] = A[index] + dA->DX[index]*half_stepX;

            MidSpaceArr->YB[index] = A[transY_index] - dA->DY[transY_index]*half_stepY;
            MidSpaceArr->YT[index] = A[index] + dA->DY[index]*half_stepY;
        }
    }
}

void ComputeMidSpaceStepForEnergy(MidSpaceStepArrays *MidSpaceE, MidSpaceStepArrays *MidSpacePres, MidSpaceStepArrays *MidSpaceDens, MidSpaceStepArrays *MidSpaceVx, MidSpaceStepArrays *MidSpaceVy){
    int i,j,index;
    float half_stepX = 0.5*dx;
    float half_stepY = 0.5*dy;
    float gas_c_1 = gas_c - 1.0;
    float E_XL, E_XR, E_YB, E_YT;

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);

            E_XL = MidSpacePres->XL[index]/gas_c_1 + 0.5*MidSpaceDens->XL[index]*(powf(MidSpaceVx->XL[index], 2.0) + powf(MidSpaceVy->XL[index], 2.0));
            E_XR = MidSpacePres->XR[index]/gas_c_1 + 0.5*MidSpaceDens->XR[index]*(powf(MidSpaceVx->XR[index], 2.0) + powf(MidSpaceVy->XR[index], 2.0));
            E_YB = MidSpacePres->YB[index]/gas_c_1 + 0.5*MidSpaceDens->YB[index]*(powf(MidSpaceVx->YB[index], 2.0) + powf(MidSpaceVy->YB[index], 2.0));
            E_YT = MidSpacePres->YT[index]/gas_c_1 + 0.5*MidSpaceDens->YT[index]*(powf(MidSpaceVx->YT[index], 2.0) + powf(MidSpaceVy->YT[index], 2.0));

            MidSpaceE->XL[index] = E_XL;
            MidSpaceE->XR[index] = E_XR;

            MidSpaceE->YB[index] = E_YB;
            MidSpaceE->YT[index] = E_YT;
        }
    }
}