#ifndef BCGC_H
#define BCGC_H

#include <iostream>

namespace bCGC {
    inline double lambda = 0.2063; // valor típico para lambda em modelos bCGC
    inline double gamma_s = 0.6599; // valor típico para gamma_s em modelos bCGC
    inline double x_0 = 0.00105; // valor típico para x0 em modelos bCGC
    inline double kappa = 9.9;
    inline double N_0 = 0.3358;
    inline double B_CGC = 5.5; // GeV^-2


    double QS( double x, double b);

    double N_IIM(double r,  double x);

    double prof(double r, double x, double b);

    double sigma_qq(double r, double x);

    double amplitude_p(double x, double Q2, const Meson& M);

    void N_csv();

    void sigma_p_csv();
}






#endif// BCGC_H