#include "sdef_misc.h"
#include <math.h>

double interpValueAt(double x,double T[],double Y[], unsigned int length) {
	/*find closest index c of where x is located in T */
	unsigned long int d = 0;
	for (unsigned long int i = 1; i < length; i++) {
		if (fabs(T[i]-x) < fabs(T[i-1]-x)) {
			d = i;
		}
	}
	/*Use this index to perform linear interpolation in Y */
	return Y[d]+((Y[d+1]-Y[d])/(T[d+1]-T[d]))*(x-T[d]);
}

