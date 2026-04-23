#ifndef DGLAP_HPP
#define DGLAP_HPP

#include "dipoleamplitude.hpp"
#include "utils.h"
#include "ctes.h"

using namespace MZ_ipsat;

namespace DGLAP{
double amplitude(double x, double Q2, const Meson& M,
                      DipoleAmplitude& dipole,
                       int Nr , 
                     double rmin , double rmax );

void sigma_p_csv();

void N_p_csv();

string N_csv(double x);
}
#endif // DGLAP_HPP