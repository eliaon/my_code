#ifndef MINUIT_HPP
#define MINUIT_HPP

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>
#include <iostream>


void minimizar_chi2();

double calc_chi2(double omega, double B);

#endif // MINUIT_HPP