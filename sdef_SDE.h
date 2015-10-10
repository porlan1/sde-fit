#ifndef SDEF_SDE_H
#define SDEF_SDE_H

#include <gsl/gsl_rng.h>

double mu(double y, double parms[]);
double sigma(double y, double parms[]);
double sigma_x(double y, double parms[]);
double milstein_scheme(double y,double dt, double parms[], gsl_rng *rand);
void simulate_SDE(double y0, double tf, double dt, double Y[], double T[], double parms[],gsl_rng *rand);

#endif /* sdef_sde.h */