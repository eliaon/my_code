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

namespace IPSAT{

double sigma_dipolo(double r, double x, double Delta, MZ_ipsat::DipoleAmplitude& dipole);

double amplitude(double x, double Delta, double Q2, const Meson& M,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);

double dsigma_dt(double x, double Q2, const Meson& M, double t,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);


double sigma_slope(double x, double Q2, const Meson& M,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);

double sigma_integrado(double x, double Q2, const Meson& M,
                 int Nr, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);

void dsigma_dt_csv(double W, const Meson& M_GLC);

void dsigma_dump(void);

std::string N_csv(double x);

void sigma_integrado_csv(void);
};

#endif //IPSAT_H