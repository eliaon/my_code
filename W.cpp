#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <vector>
#include <boost/math/special_functions/bessel.hpp>

// ---------------- Constantes globais ----------------
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double Q2    = 0.0;
const double alfem = 1.0 / 137.0;
const double Nc    = 3.0;
const double CFAC  = 5.07;     // conversão GeV^-1 -> fm
const double CFAC2 = CFAC*CFAC;

// ---------------- Classe Meson ----------------
class Meson {
public:
    std::string nome;
    double MV;   // massa do méson
    double mf;   // massa do quark
    double ef;   // carga efetiva
    double NT;   // normalização transversa
    double NL;   // normalização longitudinal
    double R2T;  // parâmetro transverso GLC
    double R2L;  // parâmetro longitudinal GLC
    double R2;   // parâmetro BG
    double B;    // parâmetro B
    bool isGLC;  // se é GLC ou BG

    // Construtor Gaus-LC
    Meson(std::string n, double MV_, double mf_, double ef_,
          double NT_, double R2T_, double NL_, double R2L_)
        : nome(n), MV(MV_), mf(mf_), ef(ef_), NT(NT_), NL(NL_),
          R2T(R2T_), R2L(R2L_), R2(0.0), isGLC(true)
    {
        B =  0.60 * (14.0 / std::pow((Q2 + MV*MV) / CFAC2, 0.26) + 1.0);
    }

    // Construtor Boosted Gaussian
    Meson(std::string n, double MV_, double mf_, double ef_,
          double NT_, double NL_, double R2_)
        : nome(n), MV(MV_), mf(mf_), ef(ef_), NT(NT_), NL(NL_),
          R2T(0.0), R2L(0.0), R2(R2_), isGLC(false)
    {
        B =  0.60 * (14.0 / std::pow((Q2 + MV*MV) / CFAC2, 0.26) + 1.0);
    }
};

// ---------------- phi_T e derivada ----------------
inline double phi_T(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        double zz = z * (1.0 - z);
        return M.NT * (zz * zz) * std::exp(-(r * r) / (2.0 * M.R2T));
    } else {
        double zz = z * (1.0 - z);
        double part1 = M.NT * zz;
        double arg1  = -(M.mf*M.mf * M.R2) / (8.0 * zz);
        double arg2  = -(2.0 * zz * r * r) / M.R2;
        double arg3  = (M.mf*M.mf * M.R2) / 2.0;
        return part1 * std::exp(arg1 + arg2 + arg3);
    }
}

inline double dphi_dr(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        return -(r / M.R2T) * phi_T(r, z, M);
    } else {
        double zz = z * (1.0 - z);
        return -((4.0 * r * zz) / M.R2) * phi_T(r, z, M);
    }
}

// ---------------- integrando psi_V psi_T ----------------
inline double psi_Vpsi_T(double z, double r, const Meson& M)
{
    double Mf2 = M.mf * M.mf;
    double EPS2 = z * (1.0 - z) * Q2 + Mf2;
    double EPS  = std::sqrt(EPS2);

    double K0 = boost::math::cyl_bessel_k(0, EPS * r);
    double K1 = boost::math::cyl_bessel_k(1, EPS * r);

    double PHIT  = phi_T(r, z, M);
    double DPHIT = dphi_dr(r, z, M);

    double ZZZ = z*z + (1.0 - z)*(1.0 - z);
    double ANORM = M.ef * std::sqrt(4.0 * M_PI * alfem) * Nc / (M_PI * z * (1.0 - z));

    return ANORM * ((Mf2 * K0 * PHIT) - (ZZZ * EPS * K1 * DPHIT));
}

// ---------------- integração numérica (Simpson) ----------------
double integrate_simpson(const std::function<double(double)>& f,
                         double a, double b, int n = 1000)
{
    if (n < 2) n = 2;
    if (n % 2 != 0) ++n;
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        const double x = a + i * h;
        s += f(x) * ((i % 2 == 0) ? 2.0 : 4.0);
    }
    return s * h / 3.0;
}

// ---------------- W(r) ----------------
inline double W_r(double r, const Meson& M)
{
    auto integrand = [r, &M](double z) { return psi_Vpsi_T(z, r, M); };
    double integral = integrate_simpson(integrand, 1e-6, 1.0 - 1e-6);
    return r * integral / 2.0;
}

// ---------------- principal ----------------
int main()
{
    const double q = 2.0 / 3.0; // carga do quark estranho
    // J/psi nos dois modelos
    // Meson(std::string n, double MV_, double mf_, double ef_,
         // double NT_, double R2T_, double NL_, double R2L_)
        //: nome(n), MV(MV_), mf(mf_), ef(ef_), NT(NT_), NL(NL_),
        //  R2T(R2T_), R2L(R2L_), R2(0.0), isGLC(true)
    {
    Meson Jpsi_GLC("Jpsi_GLC", 3.097, 1.4, 2.0/3.0, 1.23, 6.5, 0.0, 0.0);
    Meson Jpsi_BG ("Jpsi_BG",  3.097, 1.4, 2.0/3.0, 0.578, 0.575, 2.3);
    Meson phi_GLC("phi_GLC", 1.019, 0.14, 1.0/3.0, 4.75, 16.0, 0.0, 0.0);
    Meson phi_BG("phi_BG", 1.019, 0.14, 1.0/3.0, 0.919, 0.825, 11.2);

    std::ofstream fout("W_results_phi.csv");
    fout << "r_fm,W_GLC,W_BG\n";

    const int N = 4000;
    const double rmin = 0.001, rmax = 100.0;

    for (int i = 0; i < N; ++i) {
        double frac = static_cast<double>(i) / (N - 1);
        double r = rmin * std::pow(rmax / rmin, frac); // escala log
        fout << r/CFAC << "," 
             << W_r(r, phi_GLC) << "," 
             << W_r(r, phi_BG) << "\n";
    }

    fout.close();
    std::cout << "Arquivo 'W_results_phi.csv' gerado." << std::endl;
    return 0;
}}
