#include "Parameters.h"
#include "Arrays.h"
#include "Index.h"
#include <stdio.h> 

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

            printf("%d, %d", transX_index, transY_index);

            transX_index = GetIndex(j_upper, i, Nx);
            transY_index = GetIndex(j, i_upper, Nx);

            MidSpaceArr->XL[index] = A[transX_index] - dA->DX[transX_index]*half_stepX;
            MidSpaceArr->XR[index] = A[index] + dA->DX[index]*half_stepX;

            MidSpaceArr->YB[index] = A[transY_index] - dA->DY[transY_index]*half_stepY;
            MidSpaceArr->YT[index] = A[index] + dA->DY[index]*half_stepY;
        }
    }
}