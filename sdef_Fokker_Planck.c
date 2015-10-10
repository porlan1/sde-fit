#include "sdef_Fokker_Planck.h"
#include "sdef_SDE.h"

int func(double t, const double y[], double f[], void *params) {
    parameters par = *(parameters *)params;
	unsigned int N = par.N;
	double xmax = par.xmax;
	double r = par.r;
	double c = par.c;
	double K = par.K;
	double h = xmax/((double)(N-1));
	double x[N];  
	linspace(0.0, xmax, N, x);
	double parms[3] = {r, K, c};
	double d2Dp_dx2;
	double dmup_dx;
	
	for (int i = 0; i < N; i++) {
		d2Dp_dx2 = 0.5*(sigma(x[i+1],parms)*sigma(x[i+1],parms)*y[i+1]-2*sigma(x[i],parms)*sigma(x[i],parms)*y[i]+sigma(x[i-1],parms)*sigma(x[i-1],parms)*y[i-1])/(h*h);
		if (x[i] > K) {
			dmup_dx = -(mu(x[i+1],parms)*y[i+1]-mu(x[i],parms)*y[i])/h;
		} else {
			dmup_dx = -(mu(x[i],parms)*y[i]-mu(x[i-1],parms)*y[i-1])/h;
		}
		if (i==0) {
			dmup_dx = -mu(x[i],parms)*y[i]/h;
			d2Dp_dx2 = 0.5*(sigma(x[i+1],parms)*sigma(x[i+1],parms)*y[i+1]-sigma(x[i],parms)*sigma(x[i],parms)*y[i])/(h*h);
		}
		if (i==(N-1)) {
			dmup_dx = mu(x[i],parms)*y[i]/h;
			d2Dp_dx2 = 0.5*(-sigma(x[i],parms)*sigma(x[i],parms)*y[i]+sigma(x[i-1],parms)*sigma(x[i-1],parms)*y[i-1])/(h*h);
		}
		*(f+i) = dmup_dx + d2Dp_dx2;
	}
	*(f+N-1) = 0.0;
    return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
	parameters par = *(parameters *)params;
	unsigned int N = par.N;
	double xmax = par.xmax;
	double r = par.r;
	double c = par.c;
	double K = par.K;
	double h = xmax/((double)(N-1));
	double x[N];  
	linspace(0.0, xmax, N, x);
	double parms[3] = {r, K, c};
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, N, N);
	gsl_matrix *m = &dfdy_mat.matrix;
	gsl_matrix_set_zero(m);
	//only tridiagonal elements are non-zero:
	for (int i = 0; i < N; i++) {
		if (x[i] > K) {
			if (i == N-1) {
				gsl_matrix_set(m,i,i, (mu(x[i],parms)/h)-(sigma(x[i],parms)*sigma(x[i],parms)/(h*h)));
				gsl_matrix_set(m,i,i-1, 0.5*(sigma(x[i-1],parms)*sigma(x[i-1],parms)/(h*h)));
			} else {
				gsl_matrix_set(m,i,i, (mu(x[i],parms)/h)-(sigma(x[i],parms)*sigma(x[i],parms)/(h*h)));
				gsl_matrix_set(m,i,i-1, 0.5*(sigma(x[i-1],parms)*sigma(x[i-1],parms)/(h*h)));
				gsl_matrix_set(m,i,i+1, -(mu(x[i+1],parms)/h)+0.5*(sigma(x[i+1],parms)*sigma(x[i+1],parms)/(h*h)));
			}
		} else {
			if (i == 0) {
				gsl_matrix_set(m,i,i,-(mu(x[i],parms)/h)-(sigma(x[i],parms)*sigma(x[i],parms)/(h*h)));
				gsl_matrix_set(m,i,i+1,0.5*(sigma(x[i+1],parms)*sigma(x[i+1],parms)/(h*h)));
			} else {
				gsl_matrix_set(m,i,i,-(mu(x[i],parms)/h)-(sigma(x[i],parms)*sigma(x[i],parms)/(h*h)));
				gsl_matrix_set(m,i,i-1,(mu(x[i-1],parms)/h)+0.5*(sigma(x[i-1],parms)*sigma(x[i-1],parms)/(h*h)));
				gsl_matrix_set(m,i,i+1,0.5*(sigma(x[i+1],parms)*sigma(x[i+1],parms)/(h*h)));
			}
		}
	}
	for (int i = 0; i < N; i++) {
		dfdt[i] = 0;
	}
	return GSL_SUCCESS;
}

void linspace(double min, double max, int points, double *c) {
	double h = (max-min)/(points-1);
	for (int i = 0; i < points; i++) {
		*(c+i) = min + h*i;
	}
}

double trapz(double x[], double y[], int n) {
	double z = 0.0;
	double triangle = 0.0;
	double rectangle = 0.0;
	for (int i = 1; i < n; i++) {
		rectangle = y[i-1]*(x[i]-x[i-1]);
		triangle = (1/(float)2)*(x[i]-x[i-1])*(y[i]-y[i-1]);
		z = z + triangle + rectangle;
	}
	return z;
}

void elementWise(double a[], double b[], int n, double *c) {
	for (int i = 0; i < n; i++) {
		*(c+i) = a[i]*b[i];
	}
}

void solve_FP(double *Sol, parameters par, double ic_mean, double ic_sd, double tf) {
	
    gsl_odeiv2_system sys = {func, jac, (int)(par.N), &par};
    gsl_odeiv2_driver *d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf,1e-8, 1e-8, 0.0);
	
	double t = 0.0, t1 = tf;
	unsigned int N = par.N;
    double y[N];
	double f[N];
	double x[N];  
	linspace(0.0, par.xmax, N, x);
	for (int j = 0; j < N; j++) {
		y[j] = gsl_ran_gaussian_pdf(x[j]-ic_mean,ic_sd);
	}
	double z = trapz(x, y, N);
	for (int j = 0; j < N; j++) {
		y[j] /= z;
	}
	
	func(t,y,f,(void *)&par);
    
    for (int i = 1; i <= (int)tf; i++)
    {
        double ti = i * t1 / tf;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        
        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }
    }
	//put the output into Solution array:
	for (int i = 0; i < N; i++) {
		*(Sol+i) = y[i];
	}
    gsl_odeiv2_driver_free (d);
}
