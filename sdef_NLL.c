#include "sdef_NLL.h"
#include <math.h>
#include "sdef_SDE.h"
#include "sdef_Fokker_Planck.h"
#include "sdef_misc.h"


double NLLF_FP(Model_data MD, parameters par) {
	/* this function calculates the negative log likelihood 
	 * of a data set Y_obs, 
	 * using numerical solutions of the Fokker Planck PDE */
	double *Y_obs = MD.Y_obs;
	unsigned int Tf = MD.Tf;
	double NLL = 0.0;
	double like = 1.0;
	unsigned int N = par.N;
	double Sol[N];
	double xmax = par.xmax;
	//double h = xmax/((double)(N-1));
	double x[N];  
	linspace(0.0, xmax, N, x);
	 
	for (int i = 0; i < Tf; i++) {
		/*for observation i solve FP equation up to time of next observation*/
		solve_FP(Sol, par, Y_obs[i], 0.001, 1.0);
		/*interpolate the results to estimate the likelihood of the next observation 
		 * and add to the product of the likelihood of all other observtions */
		like *= interpValueAt(Y_obs[i+1],x,Sol,N);
		/*this solution can be slow, so uncomment the printf function
		 * if you want feedback during the solution */
		printf("dn_fp = %.3e\n",interpValueAt(Y_obs[i+1],x,Sol,N));
	}
	return -log(like);
}


double NLLF_MBB(Model_data MD){
	/*This function calculates the negative log likelihood based on a
	 * modified brownian bridge importance sampler combined with
	 * an Euler subdensity, as described on page 305 of the following article:
	 * Durham, G.B., Gallant A.R. 2002. Numerical techniques for maximum likelihood estimation of
	 * continuous time diffusion processes. Journal of Business and Economic Statistics. 20:297-338;
	 * DOI: 10.1198/073500102288618397
	 * */
	unsigned int M = MD.M;
	double *parms = MD.x;
	double delta = 1/(double)M;
	double NLL = 0.0;
	double um_old, um_new, sigma2_p, sigma2_q, ps, qs, mu_p, mu_q, L, t, Like, P, uM;
	double *Y_obs = MD.Y_obs;
	double Tf = MD.Tf;
	unsigned int K = MD.K;
	gsl_rng *rand = MD.rand;
	double time = 0.0;
	/*For each observation */
	Like = 1.0;
	for (int z = 0; z < Tf; z++){
		/*sample q times to average over a distribution */
		L = 0.0;
		for (int j = 0; j < K; j++){
			ps = 1.0;
			qs = 1.0;
			P = 1.0;
			um_old = Y_obs[z];
			uM = Y_obs[z+1];
			time = 0.0;
			/*break each sample into M time chunks */
			for (int m = 0; m < M; m++){
				t = delta*m;
				time += delta;
				/*Euler mean estimate */
				mu_p = um_old+mu(um_old,parms)*delta;
				/* Euler variance estimate */
				sigma2_p = sigma(um_old,parms)*sigma(um_old,parms)*delta;
				/* mean of importance sampler */
				mu_q = (uM-um_old)/(1-t);
				mu_q = um_old+mu_q*delta;
				/* variance of importance sampler */
				sigma2_q = ((M-m-1)/((double)(M-m)))*sigma2_p;
				/*generate a new random sample with the importance sampler */
				if (m == M-1) {
					um_new = uM;
					qs = 1;
				} else {
					um_new = mu_q + gsl_ran_gaussian(rand, sqrt(sigma2_q));
					/*probability density of this random sample with respect to importance sampler distribution */
					qs = gsl_ran_gaussian_pdf(um_new - mu_q, sqrt(sigma2_q));
				}
				/*probability density of the random sample */
				ps = gsl_ran_gaussian_pdf(um_new - mu_p, sqrt(sigma2_p));
				um_old = um_new;
				P *= ps/qs;
			}
			L += P/(double)K;
		}
		Like *=L;
	}
	NLL = -log(Like);
	return NLL;
}
