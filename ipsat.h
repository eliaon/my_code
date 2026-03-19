#ifndef IPSAT_H
#define IPSAT_H

#include "utils.h"
#include "ctes.h"
#include "wavefunctions.h"
#include "correcs.h"
#include "integration.hpp"
#include "dglap_cpp/EvolutionLO_nocoupling.h"
#include <cmath>

namespace MZ_ipsat {
    class DipoleAmplitude;
}

double sigma_dipolo_ipsat(double r, double x, double Delta, MZ_ipsat::DipoleAmplitude& dipole);

double amplitude_ipsat(double x, double Delta, double Q2, const Meson& M,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);

double dsigma_dt(double x, double Q2, const Meson& M, double t,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);


double sigma_ipsat_slope(double x, double Q2, const Meson& M,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);

double sigma_ipsat_integrado(double x, double Q2, const Meson& M,
                 int Nr, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);

void dsigma_dt_csv(double W, const Meson& M_GLC);

void dsigma_dump(void);

void N_ipsat_csv(void);

#endif //IPSAT_H