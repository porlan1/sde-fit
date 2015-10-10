//
//  sdef_ODE.h
//
//  Created by Paul Orlando on 8/19/15.
//  Copyright (c) 2015 Paul Orlando. All rights reserved.
//
#ifndef SDEF_ODE_H
#define SDEF_ODE_H

int func_ode(double t, const double y[], double f[], void *params);
void simulate_ODE(double *Y_deterministic, double *parms, double y0, unsigned int Tf);

#endif /* sdef_ODE.h */