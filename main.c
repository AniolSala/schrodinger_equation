/*
    TODO:
    1) Change description
    2) Error analysis
    3) 3D plot
    4) Change X axis label from `iteration` to `time (ev^-1)`

    PROGRAM WRITTEN BY ANIOL SALA (github.com/AniolSala)

    Program that solves the time-dependent Scrödinger equation
    in 1D using the leapfrog algorithm.

    **************************************************************

    We will work in natural units (NU): h_bar = c = 1. Then the
    unit of time is defined as h_bar / 1eV, and the unit of length
    as h_bar · c / 1eV.

    Note that for the stability of the solutions, the relation
    between dx and dt must be:

    dx² / dt >= ...

    where dt and dx are the time and space intervals.

    The value of x goes between -L/2 and L/2, and thus the space
    dimension has a length of L. The value of time goes from 0 to
    tmax.

    Note also that the function slopeLF defined in functions.c
    computes one step of the leapfrog algorithm with a potential
    equal to 0. To add a potential:

        phi2[j] = phi0[j] + slopeLF(i, j) - 2 * dt * I * V[i, j]

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
const double m = 1;      // eV
const double L = 1.;     // eV^-1 (NU)
const double tmax = .02; // eV^-1 (nu)

int main() {

    double dx, dt;
    unsigned int nx, nt;
    dx = 0.001;
    // dt = (double)1 / ((double)2 / (m * dx * dx) + L * L / (double)4);
    dt = m * dx * dx / (double)2;
    nx = L / dx + 1;
    nt = tmax / dt + 1;

    printf("\n----------PARAMS----------------------------\n");
    printf("L = %lf, tmax = %lf\n", L, tmax);
    printf("dx = %lf, dt = %1.10lf, dv = %lf\n", dx, dt, dx / dt);
    printf("nx = %u, nt = %u\n", nx, nt);
    printf("dx² / dt = %lf, h / 2m = %lf\n", dx * dx / dt, 2 / m);
    printf("--------------------------------------------\n\n");

    // TODO: phi description
    cnum *phi0 = malloc(nx * sizeof(cnum)), *phi1 = malloc(nx * sizeof(cnum)),
         *phi2 = malloc(nx * sizeof(cnum));
    if (!phi0 || !phi1 || !phi2) {
        printf("Error: Malloc phi\n");
        exit(1);
    }

    // Boundary conditions: phi(t, 0) = phi(t, L) = 0
    phi0[0] = 0 + 0 * I;
    phi1[0] = 0 + 0 * I;
    phi2[0] = 0 + 0 * I;
    phi0[nx - 1] = 0 + 0 * I;
    phi1[nx - 1] = 0 + 0 * I;
    phi2[nx - 1] = 0 + 0 * I;

    // Parameters of the initial shape
    double sigma = L / 30., mu = -L / 3.5;
    double k0 = m * 20.;

    // For quantum tunneling
    double V0 = k0 * k0 / (2. * m) + 1.; // Value of the potential wall
    unsigned dj = (10. * L / 100.) / dx; // (index) width of potential wall

    double T =
        exp(-2. * (double)dj * dj * sqrt(2. * m * (V0 - k0 * k0 / (2. * m))));
    printf("dj = %u, dx = %lf, T = %lf\n", dj, (double)dj * dj, T);

    // Initial conditions: gaussian centered at the origin
    // The leapfrog method requires to define phi0 and also phi1
    for (unsigned j = 0; j < nx; j++) {
        double x = -L / 2 + j * dx;
        // weights: .5, .3, .2
        phi0[j] =
            sqrt(gaussian(x, sigma, mu)) * (cos(k0 * x) + I * sin(k0 * x));
    }
    for (unsigned j = 1; j < nx - 1; j++) {
        // double x = -L / 2 + j * dx;
        double pot = 0;
        if (j >= 3 * nx / 4 && j <= 3 * nx / 4 + dj) {
            pot = V0;
        }
        phi1[j] = phi0[j] + (double)0.5 * slopeLF(dt, dx, phi0, j, m) -
                  2. * dt * I * phi0[j] * pot;
    }

    // File to save the solution
    FILE *output1 = fopen("../data/schrEq1D_gaussian_wpacketpot.txt", "w");
    if (!output1) {
        printf("Error: Couldn't open file\n");
        exit(1);
    }

    // Write initial conditions
    printf("Writing initial conditions...\n");
    for (unsigned j = 0; j < nx; j++) {
        writeSol(output1, phi0, 0, j, dt, dx, 0, -L / 2);
    }
    for (unsigned j = 0; j < nx; j++) {
        writeSol(output1, phi1, 1, j, dx, dt, 0, -L / 2);
    }
    printf("Done\n\n");

    // Solve equations
    printf("Solving equations...\n");
    unsigned nPlots = 800; // Write this number of plots into the file
    if (nPlots >= nt) {
        printf("Error: Invalid number (%u) of nPlots\n", nt / nPlots);
    }
    for (unsigned i = 2; i < nt - 1; i++) {
        // Next value of phi
        for (unsigned j = 1; j < nx - 1; j++) {
            // double x = -L / 2 + j * dx;
            double pot = 0;
            if (j >= 3 * nx / 4 && j <= 3 * nx / 4 + dj) {
                pot = V0;
            } else {
                phi2[j] = phi0[j] + slopeLF(dt, dx, phi1, j, m) -
                          2. * I * dt * phi0[j] * pot;
            }
        }
        // Update values
        for (unsigned j = 0; j < nx; j++) {
            phi0[j] = phi1[j];
            phi1[j] = phi2[j];
        }
        // Write the solution to the file
        if (i % (nt / nPlots) == 0) {
            for (unsigned j = 0; j < nx; j++) {
                writeSol(output1, phi2, i, j, dt, dx, 0, -L / 2);
            }
        }
    }
    printf("Done\n\n");

    // Release memory
    printf("Releasing memory...\n");
    fclose(output1);
    free(phi0);
    free(phi1);
    free(phi2);
    phi0 = NULL;
    phi1 = NULL;
    phi2 = NULL;
    printf("Done\n\n");
    printf("Program terminated with exit\n\n");
}
