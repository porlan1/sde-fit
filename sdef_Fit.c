#include "sdef_Fokker_Planck.h"
#include "sdef_ODE.h"
#include "sdef_SDE.h"
#include "sdef_siman.h"
#include "sdef_misc.h"
#include "sdef_NLL.h"
#include "sdef_plotting.h"
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <string.h>

/*in terminal compile this file with:
 * gcc -lgsl -lm -lblas -o sdef sdef_Fit.c sdef_Fokker_Planck.c sdef_ODE.c sdef_SDE.c sdef_siman.c sdef_misc.c sdef_NLL.c sdef_plotting.c
 * 
 * run with:
 * ./sdef
 * */


int main() {
	FILE *obs = fopen("logistic_observations.csv","w");
	FILE *fit = fopen("logistic_fit.csv","w");
	Model_data *Current = (Model_data *)malloc(sizeof(Model_data));
	Model_data *Best = (Model_data *)malloc(sizeof(Model_data));
	Model_data *New = (Model_data *)malloc(sizeof(Model_data));
	gsl_rng *rand = gsl_rng_alloc(gsl_rng_ranlux389);
	gsl_rng_set(rand,71);
	
	double upper_bound[3] = {3.0, 6.0, 1.0};
	double step_size[3] = {1e-1,5e-1,1e-1};
	
	sa_parms sap;
	sap.tmin = 2e-5;
	sap.mu = 1.005;
	sap.startTemp = 1;
	sap.no_iter = 200;
	sap.k = 1;
	sap.step_size = step_size;
	sap.upper_bound = upper_bound;
	
	const unsigned int length = 4001;
	double parms[3] = {0.05,5.3,0.04};
	double T[length] = {0};
	double Y[length] = {0};
	double dt = 0.01;
	double y0 = 1.0;
	double tf = 40.0;
	const unsigned int Tf = 40;
	double Y_obs[Tf+1] = {0};
	unsigned int K = 8;
	unsigned int M = 8;
	double NLL_MBB = 0.0;
	double NLL_FP = 0.0;
	double NLL_ens = 0.0;
	double Y_deterministic[Tf+1];
	
	simulate_SDE(y0,tf,dt,Y,T,parms,rand);

	/*capture discrete observations using interpValueAt function*/
	for (int i = 0; i < (Tf+1); i++) {
		Y_obs[i] = interpValueAt((double)i,T,Y,length);
		fprintf(obs,"%d %f\n",i,Y_obs[i]);
		printf("Y_obs[%d] = %f\n",i,Y_obs[i]);
	}
	
	(*Best).Y_obs = Y_obs;
	(*Best).M = M;
	(*Best).K = K;
	(*Best).Tf = Tf;
	(*Best).rand = rand;
	(*Best).x[0] = 1.0;
	(*Best).x[1] = 5.0;
	(*Best).x[2] = 0.2;
	memcpy((void *)Current, (void *)Best, sizeof(Model_data));
	
	/*run simulated annealing algorithm to find optimal parameter fit */
	simulated_annealing(Best, Current, New, sap);
	
	/*for verification a comparison of the negative log likelihoods generated by different methods */
	parameters par;
	par.N = 201;
	par.xmax = 8.0;
	par.r = Best->x[0];
	par.K = Best->x[1];
	par.c = Best->x[2];
	NLL_FP = NLLF_FP(*Best, par);
	
	for (int i = 0; i < Tf+1; i++) {
		printf("%d %f\n",i,Best->Y_obs[i]);
	}
	
	NLL_MBB = NLLF_MBB(*Best);
	
	printf("nLL_FP = %f\nNLL_MBB = %f\n",NLL_FP,NLL_MBB);
	
	/*solution of the logistic ODE at best parameters */
	simulate_ODE(Y_deterministic, Best->x, y0, Tf);
	
	/*see what the negative log likelihood surface looks like */
	threeD_plot_NLL(Best, sap, parms);
	plotFit(Best, parms);
	
	fclose(fit);
	fclose(obs);
	gsl_rng_free(rand);
	free(New);
	free(Current);
	free(Best);
	return 0;

}


