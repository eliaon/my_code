#include <cmath>
#include "utils.h"

#ifndef GBW_H
#define GBW_H

// ---------------- Classe de parâmetros GBW --------------
class parametros_GBW {
    public:
        double sigma0; // GeV^-2
        double x0;
        double lambda;

        parametros_GBW(double sigma0_, double x0_, double lambda_)
            : sigma0(sigma0_), x0(x0_), lambda(lambda_) {}
};

extern parametros_GBW gbw;

double QS2_GBW( double x, parametros_GBW params);

double N_GBW(double r,  double x, parametros_GBW params);

double sigma_qq_GBW(double r, double x, parametros_GBW params);

void run_sigma_csv_GBW();

double amplitude_GBW(double x, double Q2, const Meson& M,
                const parametros_GBW& gbw,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);

double sigma_x_GBW(double x, double Q2 , const Meson& M,
               int Nr = 600, int Nz = 200,
               double rmin = 1e-4, double rmax = 10.0);





#endif// GBW_H