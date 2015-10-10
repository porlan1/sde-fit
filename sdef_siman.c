#include "sdef_siman.h"
#include "sdef_NLL.h"
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
#include "sdef_Fokker_Planck.h"
#include "sdef_SDE.h"


void simulated_annealing(Model_data *Best, Model_data *Current, Model_data *New, sa_parms sap) {
	/* a simulated annealing algorithm very similar to the gsl algorithm described here:
	 * https://www.gnu.org/software/gsl/manual/html_node/Simulated-Annealing.html
	 * This function is a global optimization method that randomly searches to the parameter space to 
	 * find the optimal set of parameters that minimizes the negative log likelihood function */
	double bolt = 0.0;
	double tmin = sap.tmin;
	double mu = sap.mu;
	double Temp = sap.startTemp;
	unsigned int no_iter = sap.no_iter;
	double k = sap.k;
	double *step_size = sap.step_size;
	double *upper_bound = sap.upper_bound;
	double rand = 0.0;
	double Best_E = 0.0;
	double Current_E = 0.0;
	double New_E = 0.0;
	unsigned int z;
	gsl_rng *r = gsl_rng_alloc(gsl_rng_uni);
	gsl_rng_set(r,10);
	gsl_rng *r1 = gsl_rng_alloc(gsl_rng_uni);
	gsl_rng_set(r1,22);
	
	/* start of the simualation with the current point as the best */
	Current_E = Energy((void *)Current);
	Best_E = Energy((void *)Best);
	printf("Iter Temp	r	 K	    c	    Current	    Best\n");
	while (Temp > tmin) {
		z++;
		for (int i = 0; i < no_iter; i++) {
			/* reset "new" parameters to that of the current parameters */
			memmove((void *)New,(void *)Current,sizeof(Model_data));
			/*take a random step from current parameter position */
			Step(r1, (Model_data *)New, step_size, upper_bound);
			/* evaluate energy associated with new parameters */
			New_E = Energy(New); 
			if (New_E <= Best_E) {
				/* change new to best */
				memmove((void *)Best,(void *)New,sizeof(Model_data));
				Best_E = New_E;
				/* and also to current */
				memmove((void *)Current,(void *)New,sizeof(Model_data));
				Current_E = New_E;
			} else if (New_E <= Current_E) {
				/*change new to current */
				memmove((void *)Current,(void *)New,sizeof(Model_data));
				Current_E = New_E;
			} else {
				bolt = exp(-(New_E-Best_E)/(k*Temp));
				rand = gsl_rng_uniform(r);
				if (rand < bolt) {
				/*take a step leading to higher energy (for global optimization)
				 * change new to current */
				memmove((void *)Current,(void *)New,sizeof(Model_data));
				Current_E = New_E;
				}
			}
		}
		/*temperature cooling by damping factor mu */
		Temp /= mu;
		/*print "new" parameters along with current and best energy */ 
		printf("%d	  %.5e	%.2f %.2f	%.2f	%.2f		%.2f\n",z,Temp,New->x[0],New->x[1],New->x[2],Current_E,Best_E);
	}
	/*print the optimal solution */
	printf("Parameters for best NLL:\n");
	printf("r = %f, K = %f, c = %f\n",Best->x[0],Best->x[1],Best->x[2]);
	gsl_rng_free(r);
	gsl_rng_free(r1);
}

double Energy(void *New) {
	/* this function calculates the simulated annealing "energy"
	 * the simulated annealing algorithm minimizes the energy function.
	 * In this implementation, the model parameters of the new step 
	 * are passed to the NLL function*/
	return NLLF_MBB(*(Model_data *)New);
}


void Step(const gsl_rng *r, Model_data *New, double *step_size, double *upper_bound) {
	double z;
	/* this function produces a new random set of parameters
	 * by taking a step from the current set of parameters*/
	 
	/*loop over the three parameters, each has different step size and upper bound */
	for (int i = 0; i < 3; i++) {
		/*draw a random double from uniform distribution between 0 and 1 */
		z = gsl_rng_uniform(r);
		/*take a step */
		New->x[i] = New->x[i]+(2.0*z-1.0)*step_size[i];
		/*check for lower bound and correct */
		if (New->x[i] <= 0) {
			New->x[i] = 0;
			New->x[i] = New->x[i]+z*step_size[i];
		/*check for upper bound and correct */
		} else if (New->x[i] >= upper_bound[i]) {
			New->x[i] = upper_bound[i];
			New->x[i] = New->x[i]+(1.0*z-1.0)*step_size[i];
		}
	}
}
