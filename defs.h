#ifndef DEFS_H_
#define DEFS_H_

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double complex cnum;

void writeSol(FILE *input, cnum *phi, unsigned indt, unsigned indx, double dt,
              double dx, double t0, double x0);
double gaussian(double x, double sigma, double mu);
cnum slopeLF(double dt, double dx, cnum *phi, unsigned ind, double m);
double hoNormModes(double x, int ord);

#endif
