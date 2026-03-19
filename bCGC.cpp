#include <cmath>
#include "ctes.h"
#include "utils.h"
#include "integration.hpp"


// ------------ calculo do N -----------
 

 
double QS_bCGC( double x, double b, double x0)
{
    double B_CGC = 7.5; // GeV^-2
    double Qs2 = std::pow(x0/x, lambda/2.0)*
                std::pow(std::exp(-b*b/(2*B_CGC)), 1.0/(2.0*gamma_s));
    return Qs2;
}
 

 double N_IIM(double r,  double x,  double x0)
{
    const double lambda_iim = 0.119;
    const double gamma_s = 0.46;     // anomalous dimension
    const double kappa = 9.9;
    const double N_0 = 0.558;

    double Qs  = std::pow(x0/x, lambda_iim/2.0);
    double rQs = r * Qs;

    rQs = std::max(rQs, 1e-12);      // proteção numérica

    if (rQs < 2.0) {

        double Lx = std::log(1.0/x);

        double exp =
            2.0 * (gamma_s +
            std::log(2.0/rQs) /
            (kappa * lambda_iim * Lx));

        return N_0 * std::pow(rQs/2.0, exp)*std::pow(1.0-x, 5.26);

    } else {

        double a = -N_0*N_0*gamma_s*gamma_s /
                   ((1-N_0)*(1-N_0)*std::log(1.0-N_0));

        double b = 0.5 * std::pow(1.0-N_0,
                   -(1.0-N_0)/(N_0*gamma_s));

        double ln = std::log(b * rQs);

        return 1.0 - std::exp(-a * ln * ln)*std::pow(1.0-x, 5.26);
    }
}

double prof_bCGC(double r, double x, double b, double x0)
{
    double Qs = QS_bCGC(x, b, x0);
    const double lambda_iim = 0.119;
    const double gamma_s = 0.46;     // anomalous dimension
    const double kappa = 9.9;
    const double N_0 = 0.558;

    double rQs = r * Qs;

    rQs = std::max(rQs, 1e-12);      // proteção numérica

    if (rQs < 2.0) {

        double Lx = std::log(1.0/x);

        double expnt =
            2.0 * (gamma_s +
            std::log(2.0/rQs) /
            (kappa * lambda_iim * Lx));

        double N_b = N_0 * std::pow(rQs/2.0, expnt);
        return N_b;

    } else {

        double a = -N_0*N_0*gamma_s*gamma_s /
                   ((1-N_0)*(1-N_0)*std::log(1.0-N_0));

        double b = 0.5 * std::pow(1.0-N_0,
                   -(1.0-N_0)/(N_0*gamma_s));

        double ln = std::log(b * rQs);

        double N_b = 1.0 - std::exp(-a * ln * ln);
        return N_b * std::pow(1.0-x, 5.26);
    }
}

double sigma_qq_bCGC(double r, double x, double x_0)
{
    auto integrand = [&](double b) {
        return b * prof_bCGC(r, x, b, x_0);
    };
    double bmax = 10.0; // Limite superior para a integração em b
    return 4.0 * M_PI * integrate_simpson(integrand, 0.0, bmax, 200);
}
