#include "sdef_SDE.h"
#include <math.h>
#include <gsl/gsl_randist.h>


double milstein_scheme(double y,double dt, double parms[], gsl_rng *rand) {
	/*numerical solution of a stochastic differential equation for a single step
	 * using the milstein scheme */
	double Y_new;
	double dW;
	double z = gsl_ran_ugaussian(rand);
	dW = sqrt(dt)*z;
	Y_new = y + mu(y,parms)*dt + sigma(y,parms)*dW + 0.5*sigma(y,parms)*sigma_x(y,parms)*(dW*dW-dt);
	
	return Y_new;
}

double mu(double y, double parms[]) {
	/*computes the deterministic component mu(y) for the stochastic differential equation:
	 * dy = mu(y)dt + sigma(y)dW*/
	double x;
	double r = *(parms);
	double K = *(parms + 1);
	x = y*r*(1-y/K);
	return x;
}

double sigma(double y, double parms[]) {
	/*computes the stochastic component sigma(y) for the stochastic differential equation:
	 * dy = mu(y)dt + sigma(y)dW*/
	double x;
	double K = *(parms + 1);
	double c = *(parms + 2);
	//x = y*(1-y/K)*c;
	x = y*c;
	return x;
}

double sigma_x(double y, double parms[]) {
	/*computes derivative of sigma 
	 * as used for the milstein scheme*/
	double x;
	double K = *(parms + 1);
	double c = *(parms + 2);
	//x = c - c*2*y/K;
	x = c;
	return x;
}

void simulate_SDE(double y0, double tf, double dt, double Y[], double T[], double parms[], gsl_rng *rand) {
	/*this function simulates a single path from a stocastic differential equation
	 * with initial condition y0 and final time tf
	 * the solutions are added into the Y array.
	 * and times into the T array. */
	double t = 0.0;
	int i = 0;
	while (t < tf+dt) {
		if (i==0) {
			*(Y+i) = y0;
			*(T+i) = t;
		} else {
			*(Y+i) = milstein_scheme(*(Y+i-1),dt,parms,rand);
			*(T+i) = t;
		}
		t += dt;
		i++;
	}
}