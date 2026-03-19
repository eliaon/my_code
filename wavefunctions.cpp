#include "utils.h"
#include "ctes.h"
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include "integration.hpp"



// ---------------- phi_T e derivada ----------------
double phi_T(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        double zz = z * (1.0 - z);
        return M.NT * (zz * zz) * exp(-(r * r) / (2.0 * M.R2T));
    } else {
        double zz = z * (1.0 - z);
        double part1 = M.NT * zz;
        double arg1  = -(M.mf*M.mf * M.R2) / (8.0 * zz);
        double arg2  = -(2.0 * zz * r * r) / M.R2;
        double arg3  = (M.mf*M.mf * M.R2) / 2.0;
        return part1 * exp(arg1 + arg2 + arg3);
    }
}
double phi_L(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        double zz = z * (1.0 - z);
        return M.NL * zz * exp(-(r * r) / (2.0 * M.R2L));
    } else {
        double zz = z * (1.0 - z);
        double part1 = M.NL * zz;
        double arg1  = -(M.mf*M.mf * M.R2) / (8.0 * zz);
        double arg2  = -(2.0 * zz * r * r) / M.R2;
        double arg3  = (M.mf*M.mf * M.R2) / 2.0;
        return part1 * exp(arg1 + arg2 + arg3);
    }
}

double dphiT_dr(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        return -(r / M.R2T) * phi_T(r, z, M);
    } else {
        double zz = z * (1.0 - z);
        return -((4.0 * r * zz) / M.R2) * phi_T(r, z, M);
    }
}
double dphiL_dr(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        return -(r / M.R2L) * phi_L(r, z, M);
    } else {
        double zz = z * (1.0 - z);
        return -((4.0 * r * zz) / M.R2) * phi_L(r, z, M);
    }
}

double laplacian_phi_L(double r, double z, const Meson& M)
{
    const double dr = 1e-4;

    if (r < 1e-6)
    {
        double phi_p = phi_L(dr, z, M);
        double phi_0 = phi_L(0.0, z, M);

        double d2phi = (phi_p - phi_0) / (dr * dr);
        return 2.0 * d2phi;
    }

    double phi_p = phi_L(r + dr, z, M);
    double phi_m = phi_L(r - dr, z, M);
    double phi_0 = phi_L(r,       z, M);

    double dphi_dr   = (phi_p - phi_m) / (2.0 * dr);
    double d2phi_dr2 = (phi_p - 2.0 * phi_0 + phi_m) / (dr * dr);

    return d2phi_dr2 + (1.0 / r) * dphi_dr;
}



// ----------------psi_V psi_T ----------------
double psi_Vpsi_T(double z, double r,double Q2,  const Meson& M)
{
    double Mf2 = M.mf * M.mf;
    double EPS2 = z * (1.0 - z) * Q2 + Mf2;
    double EPS  = sqrt(EPS2);

    double K0 = boost::math::cyl_bessel_k(0, EPS * r);
    double K1 = boost::math::cyl_bessel_k(1, EPS * r);

    double PHIT  = phi_T(r, z, M);
    double DPHIT = dphiT_dr(r, z, M);

    double ZZZ = z*z + (1.0 - z)*(1.0 - z);
    double ANORM = M.ef * sqrt(4.0 * M_PI * alfem) * Nc / (M_PI * z * (1.0 - z));

    return ANORM * ((Mf2 * K0 * PHIT) - (ZZZ * EPS * K1 * DPHIT));
}

// ----------------psi_V*psi_L ----------------
double psi_Vpsi_L(double z, double r, double Q2, const Meson& M)
{
    const double mf2 = M.mf * M.mf;
    const double eps2 = z * (1.0 - z) * Q2 + mf2;
    const double eps  = std::sqrt(eps2);

    const double K0 = boost::math::cyl_bessel_k(0, eps * r);

    const double phiL = phi_L(r, z, M);
    const double lap_phiL = laplacian_phi_L(r, z, M);

    const double pref =
        M.ef * std::sqrt(4.0 * M_PI * alfem) * Nc / M_PI;

    const double bracket =
        M.MV * phiL
        +
        (mf2 * phiL - lap_phiL) / (M.MV * z * (1.0 - z));

    return pref
           * 2.0 * std::sqrt(Q2)
           * z * (1.0 - z)
           * K0
           * bracket;
}


// ---------------- overlap integrado em z ----------------
double overlap_r(double r, double Q2, const Meson& M, int Nz) {
    auto fz = [r, Q2, &M](double z) {
        return psi_Vpsi_T(z, r, Q2, M); //+ psi_Vpsi_L(z, r, Q2, M);
    };
    return integrate_simpson( fz, 1e-6, 1.0 - 1e-6, Nz );
}


void overlap_csv(void)
{
    std::string meson_input;
    std::cout << "Insira o meson (Jpsi, phi): ";
    std::cin >> meson_input;

    // normalização simples
    if (meson_input == "jpsi") meson_input = "Jpsi";
    if (meson_input == "Phi")  meson_input = "phi";

    auto it = meson_models.find(meson_input);
    if (it == meson_models.end()) {
        std::cerr << "Meson invalido. Usando Jpsi por padrao.\n";
        it = meson_models.find("Jpsi");
    }

    const Meson& M_GLC = it->second.M_GLC;
    const Meson& M_BG  = it->second.M_BG;

    int Nz = 200;
    std::string filename = M_GLC.meson + "_overlap_r.csv";
    std::ofstream fout(filename);
    fout << "r,overlap_GLC,overlap_BG\n";
    const int Npoints = 1000;
    double rmin = 1e-4, rmax = 20.0;
    double Q2 = 0.0;

    for (int i = 0; i < Npoints; ++i) {
        double frac = static_cast<double>(i) / (Npoints - 1);
        double r = rmin * pow(rmax / rmin, frac);

        double overlap_glc = 0.5*r*overlap_r(r, Q2, M_GLC, Nz);
        double overlap_bg  = 0.5*r*overlap_r(r, Q2, M_BG, Nz);

        fout << r/CFAC << "," << overlap_glc << "," << overlap_bg << "\n";
    }
    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado." << std::endl;
}
