#include "utils.h"
#include "ctes.h"
#include "wavefunctions.h"
#include "GBW.h"
#include "DGLAP.hpp"
#include "ipsat.h"
#include "bCGC.h"
#include "integration.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>


// Fatores de correção Rg e deltinha
 
double lnA(double y, double Q2, const Meson& M,
           std::string dipolemodel,
           DipoleAmplitude& dipole)
{
    double x = std::exp(-y);

    double amp = 0.0;

    if (dipolemodel == "GBW") {
        amp = GBW::amplitude(x, Q2, M, gbw);
    }
    else if (dipolemodel == "DGLAP") {
        amp = dipole_DGLAP::amplitude(x, Q2, M, dipole, 200,  1e-4, 10.0);
    }
    else if (dipolemodel == "ipsat") {
        amp = IPSAT::amplitude(x, 0.0, Q2, M, 200, 1e-4, 10.0);
    }
    else if (dipolemodel == "bCGC") {
        amp = bCGC::amplitude_p(x, Q2, M);
    }
    else {
        throw std::invalid_argument("Modelo de dipolo desconhecido: " + dipolemodel);
    }

    if (amp == 0.0) {
        std::cerr << "Amp = 0 em x=" << x << std::endl;
        return -1e30;
    }

    return std::log(std::abs(amp));
}


// -------------- lamnda_e --------------
// pipeline é a seguinte: calcula sua amplitude, usa ela para o lambda e depois o lambda para o Rg e beta, que são os fatores de correção
// fica mais leve calculando o lambda uma vez só.
double calculate_lambda(double x, double Q2, const Meson& M, std::string dipolemodel)
{
    const double h = 1e-4;
    double err;
    double y = -std::log(x);

    DipoleAmplitude dipole(MZ_IPSAT);
    dipole.EnableLookupTable();

    auto f_lnA = [&](double y_var) {
        return lnA(y_var, Q2, M, dipolemodel, dipole);
    };

    double dlnA_dy = dfridr(f_lnA, y, h, err);

    if (x > 1e-2) {
        double y_b = -log(1e-2);
        double dlnA_dy_b = dfridr(f_lnA, y_b, h, err);
        
        return dlnA_dy_b; // para x > 1e-2, usa o valor de lambda em x=1e-2
    }

    return dlnA_dy;
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
    const Meson& M = Jpsi_GLC_GBW;

    for (int i = 0; i < 120; ++i) {
        double xi = x + i * 1e-4; // varre x de 1e-4 a 1e-2
        double amp = GBW::amplitude(xi, Q2, M, gbw);
        double lambda_e = calculate_lambda(xi, Q2, M, "GBW");
        double Rg = RG(xi, Q2, lambda_e, M);
        double beta_val = beta(xi, Q2, lambda_e, M);

        std::cout << "x: " << xi << "  lambda_e: " << lambda_e << "  Rg: " << Rg << "  beta: " << beta_val << std::endl;
    }
}


// --------------- correção não perturbativa para o phi --------------

double f_c(double r,  double B, double omega, double R)
{
    //[R]  = GeV^-1
    double omega2 = omega * omega;
    double fc_num = 1.0 + B * std::exp(-omega2*(r - R)*(r - R));
    double fc_den = 1.0 + B * std::exp(-omega2*R*R);
    return fc_num / fc_den;
}