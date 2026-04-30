#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <iomanip>

#include "ctes.h"
#include "utils.h"

namespace DGLAP{

    // -------------- CONSTANTES DGLAP ----------------
    const double lambda_QCD = 0.2; // GeV
    const double lambda2 = 0.04; // GeV^2, para evitar singularidades em mu^2 muito pequenas

    // parâmetros para a distribuição inicial de gluons (golec-biernat's Saturation model of DIS: an update (2018))
    const double A_g = 1.18; // normalização da distribuição inicial de gluons
    const double lambda_g = 0.11; // controle do comportamento em x pequeno
    const double n = 5.6; // expoente do ( 1 - x ) ^ n na distribuição inicial de gluons

    // -------------- FUNÇÕES DGLAP ----------------

    // -------------- ALPHA_S (1-LOOP) ----------------
    double alpha_s(double mu2)
    {
        double beta0 = (11.0 * Nc - 2.0 * Nf) / 3.0;
        return 4.0 * M_PI / (beta0 * std::log(mu2 / lambda2));
    }

    // -------------- SPLITTING FUNCTION P_{gg} ----------------

    double P_gg(double z)
    {
        double eps = 1e-6; // para evitar singularidades
        if (z < eps) z = eps;
        if (z > 1.0 - eps) z = 1.0 - eps;

        return 6.0 *( (1.0 / (1.0 - z)) + (1.0 / z) - 2.0 + z * (1.0 - z) );
    }

    // -------------- DISTRIBUIÇÃO INICIAL DE GLUONS ----------------

    double xg0(double x)
    {
        // parametrização : A_g * x^{-lambda_g} * (1-x)^{n_g}
        return A_g * std::pow(x, -lambda_g) * std::pow(1.0 - x, n);
    }
    // -------------- CONVOLUÇÃO PARA EVOLUÇÃO DGLAP ----------------

    double convolution(double x, double mu2,
                   const std::function<double(double)> &xg)
{
    int Nz = 200;
    double zmin = x;
    double zmax = 1.0;

    double dz = (zmax - zmin) / Nz;

    double sum = 0.0;
    double xg_x = xg(x);

    double eps = 1e-6;

    for (int i = 0; i < Nz; ++i)
    {
        double z = zmin + (i + 0.5) * dz;

        if (z > 1.0 - eps) continue;

        double xp = x / z;
        if (xp > 1.0) continue;

        double fz = xg(xp) / z;

        double denom = (1.0 - z);
        if (denom < 1e-8) continue;

        double plus_term = (fz - xg_x) / denom;

        double regular =
            (1.0 / z) * xg(xp)
            - 2.0 * fz
            + z * (1.0 - z) * fz;

        sum += 6.0 * (plus_term + regular);
    }

    sum *= dz;

    // log term protegido
    if (x < 1.0 - 1e-12)
        sum += 6.0 * xg_x * std::log(1.0 - x);

    // delta term
    double beta0 = (11.0 * Nc - 2.0 * Nf) / 3.0;
    sum += (beta0 / 2.0) * xg_x;

    return sum;
}

    struct DGLAPTable
    {
        std::vector<double> xgrid;
        std::vector<double> mu2grid;
        std::vector<std::vector<double>> xg; // xg[imu][ix]
    };

    DGLAPTable build_table(int Nx, int Nmu,
                       double xmin, double xmax,
                       double mu2_min, double mu2_max)
{
    DGLAPTable table;

    table.xgrid.resize(Nx);
    table.mu2grid.resize(Nmu);
    table.xg.resize(Nmu, std::vector<double>(Nx));

    // grid em x (log)
    for (int i = 0; i < Nx; ++i)
    {
        double t = (double)i / (Nx - 1);
        table.xgrid[i] = xmin * std::pow(xmax / xmin, t);
    }

    // grid em μ² (log)
    for (int i = 0; i < Nmu; ++i)
    {
        double t = (double)i / (Nmu - 1);
        table.mu2grid[i] = mu2_min * std::pow(mu2_max / mu2_min, t);
    }

    // condição inicial
    for (int ix = 0; ix < Nx; ++ix)
        table.xg[0][ix] = xg0(table.xgrid[ix]);

    // evolução passo a passo
    for (int im = 1; im < Nmu; ++im)
    {
        double mu2_prev = table.mu2grid[im-1];
        double mu2_curr = table.mu2grid[im];

        double t0 = std::log(mu2_prev);
        double t1 = std::log(mu2_curr);
        double dt = t1 - t0;

        std::vector<double>& prev = table.xg[im-1];
        std::vector<double>& curr = table.xg[im];

        curr = prev;

        for (size_t ix = 0; ix < table.xgrid.size(); ++ix)
        {
            double x = table.xgrid[ix];

            auto interp = [&](double xp) {
                for (size_t j = 0; j < table.xgrid.size() - 1; ++j)
                {
                    if (xp >= table.xgrid[j] && xp <= table.xgrid[j+1])
                    {
                        double t = (xp - table.xgrid[j]) /
                                   (table.xgrid[j+1] - table.xgrid[j]);
                        return prev[j]*(1-t) + prev[j+1]*t;
                    }
                }
                return 0.0;
            };

            double conv = convolution(x, mu2_prev, interp);

            double deriv = (alpha_s(mu2_prev)/(2.0*M_PI)) * conv;

            curr[ix] += deriv * dt;
        }
    }

    return table;
}


    double xg(double x, double mu2)
{
    // ---- encontra índices em x ----
    DGLAPTable tab = build_table(100, 50,
                                 1e-6, 1.0,
                                  1.0, 100.0);

    int ix = 0;
    while (ix < tab.xgrid.size()-2 && tab.xgrid[ix+1] < x) ix++;

    // ---- encontra índices em μ² ----
    int im = 0;
    while (im < tab.mu2grid.size()-2 && tab.mu2grid[im+1] < mu2) im++;

    double x1 = tab.xgrid[ix];
    double x2 = tab.xgrid[ix+1];

    double m1 = tab.mu2grid[im];
    double m2 = tab.mu2grid[im+1];

    double tx = (x - x1)/(x2 - x1);
    double tm = (mu2 - m1)/(m2 - m1);

    double f00 = tab.xg[im][ix];
    double f10 = tab.xg[im][ix+1];
    double f01 = tab.xg[im+1][ix];
    double f11 = tab.xg[im+1][ix+1];

    return (1-tx)*(1-tm)*f00
         + tx*(1-tm)*f10
         + (1-tx)*tm*f01
         + tx*tm*f11;
}

}

