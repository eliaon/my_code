#ifndef IIM_HPP
#define IIM_HPP

#include "ctes.h"
#include <iostream>
#include "utils.h"

class param_IIM {
    public:
    string nome;
    double sigma0, gamma_s, x_0, lambda, kappa, N_0;
    param_IIM(string nome_, double sigma0_, double gamma_s_, double x_0_, double lambda_, double kappa_, double N_0_)
        : nome(nome_), sigma0(sigma0_), gamma_s(gamma_s_), x_0(x_0_), lambda(lambda_), kappa(kappa_), N_0(N_0_) {}
};

extern param_IIM IIM_S;
extern param_IIM IIM_RS;



namespace IIM {

    double N(double r,  double x, const param_IIM& param);

    double sigma_qq(double r, double x, const param_IIM& param);

    double amplitude_p(double x, double Q2, const Meson& M);

    void N_csv();

    void sigma_p_csv();
}





#endif