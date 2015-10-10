/* sdef_Fokker_Planck.h */
#ifndef SDEF_FOKKER_PLANCK_H
#define SDEF_FOKKER_PLANCK_H

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_randist.h>


typedef struct {
	unsigned int N;
	double xmax;
	double r;
	double K;
	double c;
} parameters;


void linspace(double min, double max, int points, double *c);
double trapz(double x[], double y[], int n);
void elementWise(double a[], double b[], int n, double *c);
int func(double t, const double y[], double f[], void *params);
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params);
void solve_FP(double *Sol, parameters par, double ic_mean, double ic_sd, double tf);

#endif /* sdef_Fokker_Planck.h */