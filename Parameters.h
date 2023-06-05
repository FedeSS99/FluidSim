#ifndef PARAMS_H
#define PARAMS_H
    /* PARAMETERS LIST
    Ny, Nx : Number of cells for each dimension
    Ly, Lx : Spatial dimensions for each dimension

    P0: Initial pressure of the whole system
    D_i, V_i: Density and horizontal velocities for each
    region (i=1,2)

    gas_c : Constante de gas ideal
    */
    const int Ny = 10;
    const int Nx = 10;

    const float Ly = 2.0;
    const float Lx = 1.0;

    const float dy = Ly/Ny;
    const float dx = Lx/Nx;
    const float Volume = dx*dy;

    const float P0 = 2.5;
    const float D1 = 1.0;
    const float D2 = 2.0;
    const float V1 = 0.5;
    const float V2 = -0.5;
    const float gas_c = 1.4;
#endif