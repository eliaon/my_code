#include "utils.h"
#include "GBW.h"
#include "dipoleamplitude.hpp"
#ifndef CORRECS_H
#define CORRECS_H

using namespace MZ_ipsat;

double lnA( double y, double Q2, const Meson& M, std::string dipolemodel, DipoleAmplitude& dipole);

double calculate_lambda(double x, double Q2, const Meson& M, std::string dipolemodel);

double RG(double x, double Q2, double lambda_e, const Meson& M);

double beta(double x, double Q2, double lambda_e, const Meson& M);

void debug_correc(void);

double f_c(double r,  double B, double omega, double R = 6.8);

#endif //CORRECS_H