#include "Parameters.h"
#include "Arrays.h"

void ObtainGradient(Grad *DA, double A[]){
    int i,j, index_c, index_L, index_R, index_B, index_T;
    int i_upper,j_upper, i_lower,j_lower;
    double InvDx = 0.5/dx;
    double InvDy = 0.5/dy;

    for (i=0; i<Ny; i++){
        i_upper = i+1;
        i_lower = i-1;
        GetPeriodicIndex(&i_upper, Ny);
        GetPeriodicIndex(&i_lower, Ny);

        for (j=0; j<Nx; j++){
            j_upper = j+1;
            j_lower = j-1;
            GetPeriodicIndex(&j_upper, Nx);
            GetPeriodicIndex(&j_lower, Nx);

            index_B = GetIndex(j, i_lower, Nx);
            index_T = GetIndex(j, i_upper, Nx);

            index_L = GetIndex(j_lower, i, Nx);
            index_R = GetIndex(j_upper, i, Nx);

            index_c = GetIndex(j, i, Nx);
            DA->DX[index_c] = (A[index_R] - A[index_L])*InvDx;
            DA->DY[index_c] = (A[index_T] - A[index_B])*InvDy;
        }
    }
}