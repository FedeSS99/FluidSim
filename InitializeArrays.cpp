#include "Parameters.h"
#include "Arrays.h"
#include <math.h>
#include <random>

void Init2DArrayAtZero(int Ny, int Nx, float A[]){
    int i,j, index;

    for (i = 0; i < Ny; i++){
        for (j = 0; j < Nx; j++){
            index = i*Nx + j;
            A[index] = 0.0;
            }
        }
}

void SetRandomInitialConditions(int Ny, int Nx, float P0, float D1, float D2, float V1, float V2, ScalarArrays *Scalars){
    int i,j,index;  
    int Ny_bottom = Ny/4;
    int Ny_top = 3*(Ny/4);
    float V1rand, V2rand;

    std::random_device RD;
    std::mt19937 gen(RD());
    std::uniform_real_distribution<> UniDis(-0.05, 0.05);

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = i*Nx + j;
            V1rand = UniDis(gen);
            V2rand = UniDis(gen);

            Scalars->Pres[index] = P0;
            Scalars->Vy[index] = UniDis(gen);

            if ((Ny_bottom < i) && (i < Ny_top)){
                Scalars->Vx[index] = V1;
                Scalars->Dens[index] = D1;
            }
            else{
                Scalars->Vx[index] = V2;
                Scalars->Dens[index] = D2;
            }
        }
    }
}