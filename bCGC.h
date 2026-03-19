#ifndef BCGC_H
#define BCGC_H

#include <iostream>

double x_0 = 1.84e-6; // valor típico para x_0 em modelos bCGC

double QS_bCGC( double x, double b, double x0 = x_0);

double N_IIM(double r,  double x,  double x0 = x_0);

double prof_bCGC(double r, double x, double b, double x0 = x_0);

double sigma_qq_bCGC(double r, double x, double b, double x0 = x_0);






#endif// BCGC_H