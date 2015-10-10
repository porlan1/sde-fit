/* sdef_misc.h */
#ifndef SDEF_MISC_H
#define SDEF_MISC_H

#include <gsl/gsl_rng.h>

double interpValueAt(double x,double T[],double Y[], unsigned int length);

typedef struct {
	double *Y_obs;
	unsigned int M;
	unsigned int K;
	unsigned int Tf;
	gsl_rng *rand;
	double x[3];
} Model_data;

typedef struct {
	double tmin;
	double mu;
	double startTemp;
	unsigned int no_iter;
	double k;
	double *step_size;
	double *upper_bound;
} sa_parms;

#endif /* sdef_misc.h */