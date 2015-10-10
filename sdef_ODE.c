#include "sdef_ODE.h"
#include "sdef_SDE.h"
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


int func_ode(double t, const double y[], double f[], void *params) {
	f[0] = mu(y[0],(double *)params);
    return GSL_SUCCESS;
}


void simulate_ODE(double *Y_deterministic, double *parms, double y0, unsigned int Tf) {
	
    gsl_odeiv2_system sys = {func_ode, (void *)0, 1, (void *)parms};
    gsl_odeiv2_driver *d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,1e-8, 1e-8, 0.0);
	
	double t = 0.0;
    double y[1];
	double f[1];
	y[0] = y0;
	
	func_ode(t,y,f,(void *)parms);
    *(Y_deterministic) = y[0];

    for (int i = 1; i <= Tf; i++)
    {
        double ti = i;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        
        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }
		*(Y_deterministic+i) = y[0];
    }
    gsl_odeiv2_driver_free (d);
}
