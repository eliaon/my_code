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

extern parametros_GBW GBW_10; // Parâmetros do fit GBW com 10% de correção

struct TA_Table {
    std::vector<double> b_vals; // valores de b
    std::vector<double> TA_vals; // valores de T_A(b)
};

double interpolate_TA(double b, const TA_Table& table);

TA_Table precompute_TA(int Nb, double bmax, double R, double a, double rho0);

namespace GBW {

double QS2( double x, parametros_GBW params);

double Np(double r,  double x, parametros_GBW params);

double sigma_qq(double r, double x, parametros_GBW params);

void sigma_p_csv();

double amplitude(double x, double Q2, const Meson& M,
                const parametros_GBW& gbw,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0);

double sigma_x(double x, double Q2 , const Meson& M,
               int Nr = 600, int Nz = 200,
               double rmin = 1e-4, double rmax = 10.0);

double amplitude_model(double x,  double B, double omega,
                        const Meson& M, double Q2 = 0,  int Nr = 600,
                        int Nz = 200, double rmin = 1e-4, double rmax = 10.0);

double sigma_model(double x, double B, double omega, const Meson& M, double Q2 = 0,
                     int Nr = 600, int Nz = 200, double rmin = 1e-4, double rmax = 10.0);


double compute_rho0(int A, double R, double a);

double rho_WS(double r, double R, double a, double rho0);

double sigma_A(double x, double Q2, const Meson& M, const TA_Table& table,
                    int Nr = 600, int Nz = 200, double rmin = 1e-4, double rmax = 10.0);

void sigma_A_csv();




void rapidez_PbPb_csv();

std::string N_csv(double x);

}
#endif// GBW_H