/* sdef_siman.h */
#ifndef SDEF_SIMAN_H
#define SDEF_SIMAN_H

#include "sdef_misc.h"


double Energy(void *New);
void Step(const gsl_rng *r, Model_data *New, double *step_size, double *upper_bound);
void simulated_annealing(Model_data *Best, Model_data *Current, Model_data *New, sa_parms sap);

#endif /* sdef_siman.h */