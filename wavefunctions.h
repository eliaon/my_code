
#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H
#include "utils.h"
#include "ctes.h"



double laplacian_phi_L(double r, double z, const Meson& M);

double psi_Vpsi_T(double z, double r, double Q2, const Meson& M);

double psi_Vpsi_L(double z, double r, double Q2, const Meson& M);

double phi_T(double r, double z, const Meson& M);

double phi_L(double r, double z, const Meson& M);

double overlap_r(double r, double Q2, const Meson& M, int Nz = 200);

void overlap_csv(void);

#endif// WAVEFUNCTIONS_H