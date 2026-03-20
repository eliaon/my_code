#include "utils.h"
#include "ctes.h"
#include "wavefunctions.h"
#include "GBW.h"
#include "integration.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>


// Fatores de correção Rg e deltinha
 
 double lnA( double y, double Q2, double amp, const Meson& M)
{
     double x = std::exp(-y);
    if (amp <= 0) {
        std::cerr << "Atenção: amplitude não positiva em y=" << y << std::endl;
        return 0.0;
    }

    return std::log(amp);
}


// -------------- lamnda_e --------------
// pipeline é a seguinte: calcula sua amplitude, usa ela para o lambda e depois o lambda para o Rg e beta, que são os fatores de correção
// fica mais leve calculando o lambda uma vez só.
double calculate_lambda(double x, double Q2, double amp, const Meson& M)
{
    const double h = 1e-4;
    double err;
    double y = -std::log(x);

    if (x > 1e-2) {
        return 0.2; // valor aproximado para x muito grande
    }
    else{
         auto f_lnA = [&](double y) {
        return lnA(y, Q2, amp, M);
    };

    double dlnA_dy = dfridr(f_lnA, y, h, err);
    return dlnA_dy;   
    }
}

// ----------------- fator de correção Rg skeweness ----------------

double RG(double x, double Q2, double lambda_e, const Meson& M)
{
    return std::pow(2.0, 2.0*lambda_e + 3.0) * tgamma(lambda_e + 2.5) 
                         / (sqrt(M_PI) * tgamma(lambda_e + 4.0));
}


double beta(double x, double Q2, double lambda_e, const Meson& M){
    return std::tan(M_PI * lambda_e / 2.0);
}

// ----------------- função pra printar os valores de lambda_e e Rg ----------------
void debug_correc(void)
{
    double x = 1e-4;
    double Q2 = 0.0;
    const Meson& M = Jpsi_GLC;

    for (int i = 0; i < 100; ++i) {
        double xi = x * std::pow(10.0, i * 0.1);
        double amp = amplitude_GBW(xi, Q2, M, gbw);
        double lambda_e = calculate_lambda(xi, Q2, amp, M);
        double Rg = RG(xi, Q2, lambda_e, M);
        double beta_val = beta(xi, Q2, lambda_e, M);

        std::cout << "x: " << xi << "  lambda_e: " << lambda_e << "  Rg: " << Rg << "  beta: " << beta_val << std::endl;
    }
}


// --------------- correção não perturbativa para o phi --------------

double f_c(double r)
{
    double R = 6.8; //GeV^-1
    double B = -0.9;;
    double omega = 0.15; //GeV
    double omega2 = omega * omega;
    double fc_num = 1.0 + B * std::exp(-omega2*(r - R)*(r - R));
    double fc_den = 1.0 + B * std::exp(-omega2*R*R);
    return fc_num / fc_den;
}