#ifndef PARAMS_H
#define PARAMS_H
    /* PARAMETERS LIST
    Ny, Nx : Number of cells for each dimension
    Ly, Lx : Spatial dimensions for each dimension

    FinalT : Max value of time to end simulation

    P0: Initial pressure of the whole system
    D_i, V_i: Density and horizontal velocities for each
    region (i=1,2)

    gas_c : Ideal gas constant
    */
    const int Ny = 128;
    const int Nx = 128;

    const double Ly = 1.0;
    const double Lx = 2.0;

    const double dy = Ly/Ny;
    const double dx = Lx/Nx;
    const double Volume = dx*dy;

    const double FinalT = 5.0;

    const double P0 = 2.5;
    const double D1 = 1.0;
    const double D2 = 2.0;
    const double V1 = 0.5;
    const double V2 = -0.5;
    const double gas_c = 1.4;
#endif