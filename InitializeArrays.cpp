#include "Parameters.h"
#include "Arrays.h"
#include "Index.h"
#include <math.h>
#include <random>

void SetRandomInitialConditions(ScalarArrays *Scalars){
    int i,j,index;  
    int Ny_bottom = Ny/4;
    int Ny_top = 3*(Ny/4);
    float VxRand, VyRand;

    std::random_device RD;
    std::mt19937 gen(RD());
    std::uniform_real_distribution<> UniDis(-0.05, 0.05);

    for (i=0; i<Ny; i++){
        for (j=0; j<Nx; j++){
            index = GetIndex(j, i, Nx);
            VxRand = UniDis(gen);
            VyRand = UniDis(gen);

            Scalars->Pres[index] = P0;
            Scalars->Vy[index] = VyRand;

            if ((Ny_bottom < i) && (i < Ny_top)){
                Scalars->Vx[index] = V1 + VxRand;
                Scalars->Dens[index] = D1;
            }
            else{
                Scalars->Vx[index] = V2 + VxRand;
                Scalars->Dens[index] = D2;
            }
        }
    }
}