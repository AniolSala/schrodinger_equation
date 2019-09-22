/*
    TODO:

    PROGRAM WRITTEN BY ANIOL SALA (github.com/AniolSala)

    Program that solves the time-dependent Scrödinger equation
    in 1D using the leapfrog algorithm.

    **************************************************************

    We will work in natural units (NU): h_bar = c = 1. Then the
    unit of time is defined as h_bar / 1eV, and the unit of length
    as h_bar · c / 1eV.

    Note that for the stability of the solutions, the relation
    between dx and dt must be:

    dx² / dt <= 2 / m

    where dt and dx are the time and space intervals.

    The value of x goes between -L/2 and L/2, and thus the space
    dimension has a length of L. The value of time goes from 0 to
    tmax.

    Note also that the function slopeLF defined in functions.c
    computes one step of the leapfrog algorithm with a potential
    equal to 0. To add a potential:

        psi2[j] = psi0[j] + slopeLF(i, j) - 2 * dt * I * V[i, j]

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./defs.h"

// To convert length and time from NU to standard units:
const double h_bar = 4.135667731e-15; // eV·s
const double c = 299792458;           // m/s

// Physical constans used in the calc:
const double m = 1.;     // eV
const double L = 1.;     // eV^-1 (NU)
const double tmax = .01; // eV^-1 (nu)

int main() {
    // Name of the output file
    char *outputFileName = "../data/schrEq1D_gaussian_wpacketpot_high.txt";

    // Space / time intervals and space / time steps
    double dx, dt;
    unsigned int nx, nt;

    // Params of the potential wall:
    // V: Potential, can be 0 or V0.
    // V0: Value of the potential wall.
    // x0: Position of the barrier
    // d: Width of the barrier
    // E: Energy of the particle
    // T: Transmission coefficient
    // j0: First index of the grid containing the barrier
    // dj: Interval of the grid corresponding to the width of the barrier
    double V, V0, x0, d, E, T;
    unsigned j0, dj;

    // Params of initial conditions
    double sigma = L / 30., mu = -L / 3.5; // Gaussian
    double k0 = 100;

    // Define space interval dx and the number of nodes in the grid nx
    dx = 0.001;
    nx = L / dx + 1;

    // Parameters of the potential wall
    j0 = nx / 2;
    dj = (1. * L / 20.) / dx;
    x0 = -L / 2. + dx * j0;
    d = dx * dj;
    E = k0 * k0 / (2. * m) + 1. / (8. * m * sigma * sigma);
    V0 = E * 1.05; // 1.05; // .3
    T = exp(-2. * (double)dj * dx * sqrt(2. * m * (V0 - E)));

    // Define the time interval dt and the number of time steps nt
    dt = (double)1 / ((double)2 / (m * dx * dx) + V0);
    // dt = m * dx * dx / (double)2;  // If there is no potential
    nt = tmax / dt + 1;

    // Allocate memory for the solution
    cnum *psi0 = malloc(nx * sizeof(cnum)), *psi1 = malloc(nx * sizeof(cnum)),
         *psi2 = malloc(nx * sizeof(cnum));
    if (!psi0 || !psi1 || !psi2) {
        printf("Error: Malloc psi\n");
        exit(1);
    }

    // Boundary conditions: psi(t, 0) = psi(t, L) = 0
    psi0[0] = 0 + 0 * I;
    psi1[0] = 0 + 0 * I;
    psi2[0] = 0 + 0 * I;
    psi0[nx - 1] = 0 + 0 * I;
    psi1[nx - 1] = 0 + 0 * I;
    psi2[nx - 1] = 0 + 0 * I;

    // Print all params used
    printf("\n----------PARAMS----------------------------\n");
    printf("L = %lf, tmax = %lf\n", L, tmax);
    printf("dx = %lf, dt = %1.10lf, dv = %lf\n", dx, dt, dx / dt);
    printf("nx = %u, nt = %u\n", nx, nt);
    printf("dx² / dt = %lf, h / 2m = %lf\n", dx * dx / dt, 2 / m);
    printf("V0 = %lf, E = %lf\n", V0, E);
    printf("Potential wall:\n");
    printf("x0 = %lf, d = %lf, j0 = %u, dj = %u, dx = %lf, T = %lf\n", x0, d,
           3 * nx / 4, dj, (double)dj * dx, T);
    printf("--------------------------------------------\n\n");

    // Initial conditions: (The leapfrog method requires to define psi0 and also psi1)
    for (unsigned j = 0; j < nx; j++) {
        double x = -L / 2. + (double)j * dx;
        psi0[j] =
            sqrt(gaussian(x, sigma, mu)) * (cos(k0 * x) + I * sin(k0 * x));
    }
    for (unsigned j = 1; j < nx - 1; j++) {
        // double x = -L / 2 + j * dx;
        if (j >= j0 && j <= j0 + dj) {
            V = V0;
        } else {
            V = 0.;
        }
        psi1[j] = psi0[j] + (double)0.5 * slopeLF(dt, dx, psi0, j, m) -
                  2. * dt * I * psi0[j] * V; // Wall
    }

    // File to save the solution
    FILE *output1 = fopen(outputFileName, "w");
    if (!output1) {
        printf("Error: Couldn't open file\n");
        exit(1);
    }

    // Write initial conditions
    printf("Writing initial conditions...\n");
    for (unsigned j = 0; j < nx; j++) {
        writeSol(output1, psi0, 0, j, dt, dx, 0, -L / 2);
    }
    for (unsigned j = 0; j < nx; j++) {
        writeSol(output1, psi1, 1, j, dx, dt, 0, -L / 2);
    }
    printf("Done\n\n");

    // Solve equations
    printf("Solving equations...\n");
    unsigned nPlots = 800; // Write this number of plots into the file
    if (nPlots >= nt) {
        printf("Error: Invalid number (%u) of nPlots\n", nt / nPlots);
    }
    for (unsigned i = 2; i < nt - 1; i++) {
        // Next value of psi
        for (unsigned j = 1; j < nx - 1; j++) {
            // double x = -L / 2 + j * dx;
            if (j >= j0 && j <= j0 + dj) {
                V = V0;
            } else {
                V = 0.;
            }
            psi2[j] = psi0[j] + slopeLF(dt, dx, psi1, j, m) -
                      2. * dt * I * psi1[j] * V; // Wall
            // dt * I * psi1[j] * x * x; //  Harmonic oscillator
        }

        // Update values
        for (unsigned j = 0; j < nx; j++) {
            psi0[j] = psi1[j];
            psi1[j] = psi2[j];
        }
        // Write the solution to the file
        if (i % (nt / nPlots) == 0) {
            for (unsigned j = 0; j < nx; j++) {
                writeSol(output1, psi2, i, j, dt, dx, 0, -L / 2.);
            }
        }
    }

    printf("Done\n\n");

    // Release memory
    printf("Releasing memory...\n");
    fclose(output1);
    free(psi0);
    free(psi1);
    free(psi2);
    psi0 = NULL;
    psi1 = NULL;
    psi2 = NULL;
    printf("Done\n\n");
    printf("Program terminated with exit\n\n");
}
