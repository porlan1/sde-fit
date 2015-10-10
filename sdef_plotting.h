/* sdef_plotting.h */
#ifndef SDEF_PLOTTING_H
#define SDEF_PLOTTING_H

#include "sdef_misc.h"

void plotFit(Model_data *Best, double parms[]);
void threeD_plot_NLL(Model_data *MD, sa_parms sap, double *parms);

#endif /* sdef_plotting.h */