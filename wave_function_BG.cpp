#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>  // Funções de Bessel K0 e K1

// Constantes globais
const double pi = 4.0 * atan(1.0);
const double alfem = 1.0 / 137.0;
double rprime; // variável global como no COMMON do Fortran

// Protótipos
double overlap(double r);
double zint(double z);
double wfgt(double rprime, double zzt);

// Integração numérica simples (Simpson)
double integrate(double a, double b, int n, double (*func)(double)) {
    if (n % 2 == 1) n++; // garantir que n seja par
    double h = (b - a) / n;
    double s = func(a) + func(b);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        s += func(x) * (i % 2 == 0 ? 2.0 : 4.0);
    }
    return s * h / 3.0;
}

// Função principal
int main() {
    double mJpsi = 3.097; // GeV
    double mv = mJpsi;

    std::ofstream fout("overlap_bg_Jpsi.dat");

    double r_min = 0.001;
    double r_max = 100.0;
    int k_points = 4000;

    // Aqui nesse loop fiz o r variar em escala logarítmica
    // o motivo é que pra r pequendo a curva não estava suave
    // então o passo teria que ser MUITO pequendo pra suavizar a curva
    // o que estava fazendo o código demorar

    for (int i = 0; i < k_points; i++) {
    double frac = (double)i / (k_points - 1); // varia de 0 a 1
    double r = r_min * pow(r_max / r_min, frac); // escala logarítmica
    double resultado = overlap(r);
    double rfm = r / 5.07; // conversão para fm

        std::cout << rfm << " " << resultado << " " << i << std::endl;
        fout << rfm << " " << resultado << std::endl;
    }

    fout.close();
    return 0;
}

// Função overlap
double overlap(double r) {
    rprime = r;
    double liz = 1e-6;
    double lsz = 1.0 - 1e-6;


    // integração em z de [0,1]
    double integral = integrate(liz, lsz, 1000, zint);

    return (rprime / 2.0) * integral;
}

// Função integrando em z
double zint(double zz) {
    return wfgt(rprime, zz);
}

// Função de onda transversa Boosted Gaussian
double wfgt(double rprime, double zzt) {
    double Mf2 = pow(1.4, 2);   // massa do quark charm^2
    double ehat = 2.0 / 3.0;
    double NT = 0.578;
    double R2 = 2.3;            // GeV^-2
    double q2 = 0.0;

    double Nc = 3.0;
    double ANORM = ehat * sqrt(4.0 * pi * alfem) * Nc / (pi * zzt * (1.0 - zzt));

    double EPS2 = zzt * (1.0 - zzt) * q2 + Mf2;
    double EPS = sqrt(EPS2);

    double mK0 = gsl_sf_bessel_K0(EPS * rprime);
    double mK1 = gsl_sf_bessel_K1(EPS * rprime);

    double PHITpart1 = NT * zzt * (1.0 - zzt);
    double ARGpart1 = -(Mf2 * R2) / (8.0 * zzt * (1.0 - zzt));
    double ARGpart2 = -(2.0 * zzt * (1.0 - zzt) * rprime * rprime) / R2;
    double ARGpart3 = (Mf2 * R2) / 2.0;
    double ARGexp   = ARGpart1 + ARGpart2 + ARGpart3;
    double PHIT = PHITpart1 * exp(ARGexp);

    double ZZZ = pow(zzt, 2) + pow(1.0 - zzt, 2);
    double DPHIT = -((4.0 * rprime * zzt * (1.0 - zzt)) / R2) * PHIT;

    return ANORM * ((Mf2 * mK0 * PHIT) - (ZZZ * EPS * mK1 * DPHIT));
}
