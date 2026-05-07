#include "LHAPDF.hpp"
#include "LHAPDF/LHAPDF.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include "ctes.h"
#include "utils.h"
#include "plot.h"

DipoleModel::DipoleModel(const std::string& setname, int member)
{
    pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(setname, member));

    // valores default (ajuste conforme seu fit)
    sigma0 = 23.0;   // mb
    C = 0.29;
    mu0_2 = 1.85; // GeV²
}

DipoleModel::~DipoleModel()  = default; 

// --- PDF ---
double DipoleModel::gluon(double x, double Q2) const {
    double Q = std::sqrt(Q2);
    return pdf->xfxQ(21, x, Q) / x;  // CORRETO
}

double DipoleModel::alpha_s(double Q2) const {
    return pdf->alphasQ(std::sqrt(Q2));
}

// --- escala ---
double DipoleModel::mu2(double r) const {
    return C/(r*r) + mu0_2;
}

// --- dipolo ---
double DipoleModel::sigma_dip(double x, double r) const {
    double mu2_val = mu2(r);

    double xg = gluon(x, mu2_val);
    double as = alpha_s(mu2_val);

    double arg = (M_PI*M_PI/6.0) * r*r * as * xg;

    return sigma0 * (1.0 - std::exp(-arg));
}

double DipoleModel::N(double x, double r) const {
    return sigma_dip(x, r) / sigma0;
}

void DipoleModel::N_csv()
{

    const int Npoints = 5000;
    double rmin = 1e-4, rmax = 10.0;

    std::vector<double> x_values = {1e-5, 1e-4, 1e-3, 1e-2};

    std::vector<std::string> filenames;

    for (double x : x_values)
    {
        std::string filename = "csv/N_LHAPDF_x_" + doubleParaString(x) + ".csv";
        std::ofstream fout(filename);

        for (int i = 0; i < Npoints; ++i)
        {
            double frac = (double)i / (Npoints - 1);
            double r = rmin * std::pow(rmax / rmin, frac); // log spacing

            double N_val = N(x, r);

            fout << r * CFAC << "," << N_val << "\n"; // r em fm
        }
        filenames.push_back(filename);
        std::cout << "Gerado csv para N(r) em x = " << x << ": " << filename << std::endl;
        fout.close();
    }
    plot_N_multi(filenames, x_values, "LHAPDF - CT14lo");

}

std::string DipoleModel::xf_vs_x(double Q2, int flavor)
{
    int nflavor = flavor; // 21 para xg, 1: xd, 2: xu, 3: xs, 4: xc, 5: xb, 6: xt
    std::string fstring = flavorName(nflavor);
    std::string xf_string = flavorToString(flavor);

    std::string filename = "csv/"+fstring+"_LHAPDF_Q2=" + doubleParaString(Q2) + ".csv";
    std::ofstream fout(filename);
    fout << "x," << xf_string << "\n";
    for(double x = 1e-5; x < 1.0; x *= 1.2) {
    double Q = sqrt(Q2);
    double xf = pdf->xfxQ(nflavor, x, Q); 
    fout << x << "," << xf << "\n";
    std::cout << x << " " << xf << std::endl;
}
    return filename;
}

// --- setters ---
void DipoleModel::set_sigma0(double val) { sigma0 = val; }
void DipoleModel::set_C(double val) { C = val; }
void DipoleModel::set_mu0_2(double val) { mu0_2 = val; }