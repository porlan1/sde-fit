#include "sdef_plotting.h"
#include <stdio.h>
#include "sdef_NLL.h"
#include "sdef_SDE.h"
#include "sdef_Fokker_Planck.h"

void threeD_plot_NLL(Model_data *MD, sa_parms sap, double *parms) {
	unsigned int n = 100;
	double *x = (*MD).x;
	double NLL = 0.0;
	double increment_r = (sap.upper_bound[0] - 0.0)/((double)n - 1.0);
	double increment_K = (sap.upper_bound[1] - 0.0)/((double)n - 1.0);
	double NLL_lowest = 20000.0;
	double r,K;
	FILE *file = fopen("3d_NLL.csv","w");
	FILE *file1 = fopen("NLL_contour_plot.gnu","w");
	fprintf(file1,"set xrange[0:%f]\nset yrange [0:%f]\nunset key\nunset surface\n",sap.upper_bound[0],sap.upper_bound[1]);
	fprintf(file1,"set contour both\nset table 'cont.dat'\nsplot '3d_NLL.csv'\n");
	fprintf(file1,"unset table\nset object circle at %f,%f fillstyle solid fc rgb 'black' front size 0.05\n", *x,*(x+1));
	fprintf(file1,"set label \"point found with \\nsimulated annealing\" at 0.1,%f front font 'Arial, 15'\n", *(x+1)-0.3);
	fprintf(file1,"set palette rgbformulae -10,0,10\nset xlabel 'r' font 'Arial,20' offset 0,-2,0\nset ylabel 'K' font 'Arial, 20' offset -2,0,0\n");
	fprintf(file1,"set cblabel 'NLL' font 'Arial, 20' offset 5,0,0\nset xtic font 'Arial, 20'\nset ytic font 'Arial, 20'\n");
	fprintf(file1,"set cbtic font 'Arial, 20'\nset size square\nset term png\n");
	fprintf(file1,"set output \"Actual_r%.2f_K%.2f_c%.2f__fit_r%.2f_K%.2f_c%.2f_contours.png\"\n",*parms,*(parms+1),*(parms+2),*x,*(x+1),*(x+2));
	fprintf(file1,"plot '3d_NLL.csv' with image, 'cont.dat' w l lt -1 lw 2");
	*(x) = 0.0;
	*(x+1) = 0.0;
	for (int i = 0; i < n; i++) {
		*(x) += increment_r;
		*(x+1) = 0.0;
		for (int j = 0; j < n; j++) {
			*(x+1) += increment_K;
			NLL = NLLF_MBB(*MD);
			if(NLL<NLL_lowest) {
				NLL_lowest = NLL;
				r = *(x);
				K = *(x+1);
			}
			fprintf(file,"%f %f %f\n", *x, *(x+1), NLL);
		}
		fprintf(file,"\n");
	}
	printf("Optimal Soln: r = %f, K = %f, NLL = %f\n",r,K,NLL_lowest);
	fclose(file);
	fclose(file1);
}

void plotFit(Model_data *Best, double parms[]) {
	FILE *file = fopen("population_plot.gnu","w");
	fprintf(file,"set size square\nset xlabel 'time' font 'Arial, 20'\n");
	fprintf(file,"set ylabel 'population size' font 'Arial, 20'\nset term png\n");
	fprintf(file,"set output \"Actual_r%.2f_K%.2f_c%.2f__fit_r%.2f_K%.2f_c%.2f_population.png\"\n",parms[0],parms[1],parms[2],Best->x[0],Best->x[1],Best->x[2]);
	fprintf(file,"set key right bottom\n");
	fprintf(file,"plot 'logistic_observations.csv' with linespoints pt 7 lc rgb 'black' title 'observations', \\\n");
	fprintf(file,"'logistic_fit.csv' smooth mcsplines lw 2 lc rgb 'red' title 'fit'");
}

