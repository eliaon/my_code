#include "utils.h"
#include "GBW.h"
#ifndef CORRECS_H
#define CORRECS_H


double lnA( double y, double Q2, const Meson& M);

double calculate_lambda(double x, double Q2, const Meson& M);

double RG(double x, double Q2, double lambda_e, const Meson& M);

double beta(double x, double Q2, double lambda_e, const Meson& M);

void debug_correc(void);

double f_c(double r,  double B, double omega, double R = 6.8);

#endif //CORRECS_H