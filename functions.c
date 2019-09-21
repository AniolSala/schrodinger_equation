#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./defs.h"

void writeSol(FILE *input, cnum *phi, unsigned indt, unsigned indx, double dt,
              double dx, double t0, double x0) {
    double t = t0 + indt * dt, x = x0 + indx * dx;
    fprintf(input, "%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n", t, x,
            creal(phi[indx]), cimag(phi[indx]),
            pow(creal(phi[indx]), 2) + pow(cimag(phi[indx]), 2));
}

double gaussian(double x, double sigma, double mu) {
    double factor = sqrt((double)1 / ((double)2 * M_PI * sigma * sigma));
    return factor * exp(-(x - mu) * (x - mu) / ((double)2 * sigma * sigma));
}

cnum slopeLF(double dt, double dx, cnum *phi, unsigned ind, double m) {
    cnum s = I * dt / (m * dx * dx);
    return s * (phi[ind + 1] - (double)2 * phi[ind] + phi[ind - 1]);
}

static double H(double x, int ord) {
    // Hermite polynomials
    double value;

    switch (ord) {
    case 0:
        value = (double)1;
        break;
    case 1:
        value = (double)2 * x;
        break;
    case 2:
        value = (double)2 - (double)4 * x * x;
        break;
    case 3:
        value = (double)12 * x - (double)8 * x * x * x;
        break;
    }
    return value;
}

static int factorial(int n) {
    if (n == 0) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}

double hoNormModes(double x, int ord) {
    double factor =
        sqrt((double)1 / (sqrt(M_PI) * pow((double)2, ord) * factorial(ord)));
    double fx = exp(-x * x / (double)2) * H(x, ord);
    return factor * fx;
}
