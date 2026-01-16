#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <vector>
#include <iomanip>
#include <sstream>
#include <format>
#include <boost/math/special_functions/bessel.hpp>
#include "integration.hpp"

// ---------------- Constantes globais ----------------
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double alfem = 1.0 / 137.0;   // constante de estrutura fina
const double Nc    = 3.0;           // número de cores
const double lambda = 0.29;       // parâmetro do modelo GBW
const double gamma_s = 0.46; // parâmetro do modelo bCGC
const double CFAC = 5.07; // conversão fm para GeV^-1
const double sigma0 = 23.0/ 0.3894; // GBW old gev^-2
//const double sigma0 = 27.43 / 0.3894; // GBW new gev^-2
 // Cargas dos quarks
double qJ = 2.0/3.0;
double qS = 1.0/3.0;




// Conversão de unidades
const double GeV2_to_nb = 3.89379e5; // 1 GeV^-2 = 3.89379×10^5 nb
// Parâmetros de seção de choque
//const double sigma0_GeV2 = 0.23; // Constante sigma0 em GeV^-2. (23.0 mb = 0.23 GeV^-2)
// ---------------- Classe Meson ----------------
class Meson {
public:
    std::string meson;
    std::string nome;
    double MV;   // massa do méson (GeV)
    double mf;   // massa do quark (GeV)
    double ef;   // carga efetiva
    double NT;   // normalização transversa
    double NL;   // normalização longitudinal
    double R2T;  // parâmetro transverso GLC
    double R2L;  // parâmetro longitudinal GLC
    double R2;   // parâmetro BG
    bool isGLC;  // se é GLC ou BG

    // Construtor Gaus-LC
    Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
          double NT_, double R2T_, double NL_, double R2L_)
        : meson(m), nome(n), MV(MV_), mf(mf_), ef(ef_), NT(NT_), NL(NL_),
          R2T(R2T_), R2L(R2L_), R2(0.0), isGLC(true) {
            //if m = Jpsi {
            //    B = 4.5 + 2* 0.164 *log(W2/95);
           //     else if m = phi {
            //        B = 0.60 * (14.0 / pow((Q2 + MV*MV), 0.26) + 1.0);
           // }
          }

    // Construtor Boosted Gaussian
    Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
          double NT_, double NL_, double R2_)
        : meson(m), nome(n), MV(MV_), mf(mf_), ef(ef_), NT(NT_), NL(NL_),
          R2T(0.0), R2L(0.0), R2(R2_), isGLC(false) {}
};

// Definição dos mésons
Meson Jpsi_GLC("Jpsi", "Jpsi_GLC", 3.097, 1.4, qJ, 1.23, 6.5, 0.0, 0.0);
Meson Jpsi_BG ("Jpsi", "Jpsi_BG",  3.097, 1.4, qJ, 0.578, 0.575, 2.3);
Meson phi_GLC("phi", "phi_GLC", 1.019, 0.14, qS, 4.75, 16.0, 0.0, 0.0);
Meson phi_BG("phi", "phi_BG", 1.019, 0.14, qS, 0.919, 0.825, 11.2);



// ---------------- função para calcular B(Q2) ----------------
inline double B(double x, double Q2, const Meson& M) {
    double W = std::sqrt(M.MV*M.MV/x);
    if (M.meson == "Jpsi"){
        double B1 = 4.99 + 4.0* 0.25 *log(W/90.0);
        return B1;
    }
        else if (M.meson == "phi"){
            double B2 = 0.60 * (14.0 / pow((Q2 + M.MV*M.MV), 0.26) + 1.0);
            return B2;}
            else {
                std::cerr << "Méson desconhecido para cálculo de B: " << M.meson << std::endl;
                return 0.0;
            }
}
// ---------------- phi_T e derivada ----------------
inline double phi_T(double r, double z, const Meson& M)
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
inline double phi_L(double r, double z, const Meson& M)
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

inline double dphiT_dr(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        return -(r / M.R2T) * phi_T(r, z, M);
    } else {
        double zz = z * (1.0 - z);
        return -((4.0 * r * zz) / M.R2) * phi_T(r, z, M);
    }
}
inline double dphiL_dr(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        return -(r / M.R2L) * phi_L(r, z, M);
    } else {
        double zz = z * (1.0 - z);
        return -((4.0 * r * zz) / M.R2) * phi_L(r, z, M);
    }
}

// ----------------psi_V psi_T ----------------
inline double psi_Vpsi_T(double z, double r,double Q2,  const Meson& M)
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
inline double psi_Vpsi_L(double z, double r,double Q2,  const Meson& M)
{
    double Mf2 = M.mf * M.mf;
    double EPS2 = z * (1.0 - z) * Q2 + Mf2;
    double EPS  = sqrt(EPS2);

    double K0 = boost::math::cyl_bessel_k(0, EPS * r);

    double PHIL  = phi_L(r, z, M);

    double ANORM = M.ef * sqrt(4.0 * M_PI * alfem) * Nc / (M_PI * z * (1.0 - z));

    return ANORM * (2.0 * sqrt(Q2) * z * (1.0 - z) * K0 * PHIL);
}

// ---------------- overlap integrado em z ----------------
double overlap_r(double r, double Q2, const Meson& M, int Nz = 200) {
    auto fz = [r, Q2, &M](double z) {
        return psi_Vpsi_T(z, r, Q2, M); //+ psi_Vpsi_L(z, r, Q2, M);
    };
    return integrate_simpson( fz, 1e-6, 1.0 - 1e-6, Nz );
}


void plot_overlap(const Meson& M_GLC, const Meson& M_BG, std::string N_method, int Nz = 200)
{
    std::string filename = M_GLC.meson + "_overlap_r.csv";
    std::ofstream fout(filename);
    fout << "r,overlap_GLC,overlap_BG\n";
    const int Npoints = 1000;
    double rmin = 1e-4, rmax = 10.0;
    double Q2 = 0.0;

    for (int i = 0; i < Npoints; ++i) {
        double frac = static_cast<double>(i) / (Npoints - 1);
        double r = rmin * pow(rmax / rmin, frac);

        double overlap_glc = overlap_r(r, Q2, M_GLC, Nz);
        double overlap_bg  = overlap_r(r, Q2, M_BG, Nz);

        fout << r/CFAC << "," << overlap_glc << "," << overlap_bg << "\n";
    }
    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado." << std::endl;
}

// ------------ calculo do N -----------
double QS2_GBW(double x)
{
    const double x_0 = 3e-4;
    return pow(x_0 / x, lambda);
}   

double N_GBW(double r, double x)
{
    double Qs2 = QS2_GBW(x);
    double arg = (r * r) * Qs2 / 4.0;
    return (1.0 - exp(-arg))*std::pow(1-x, 5.26);
}

void N_plot()
{
    std::ofstream fout("csv/N_GBW.csv");
    fout << "r²,N_GBW\n";
    const int Npoints = 5000;
    double rmin = 1e-4, rmax = 10.0;
    double x = 1e-4;

    for (int i = 0; i < Npoints; ++i) {
        //double frac = static_cast<double>(i) / (Npoints - 1);
        //double r = rmin * pow(rmax / rmin, frac);

        double r = rmin + (rmax/(Npoints -1))*(i-1);
        double Nval = N_GBW(r, x);

        fout << r*r << "," << Nval << "\n";
    }

    fout.close();
}

double sigma_dipolo(double r, double x)
{
    double Nval = N_GBW(r, x);
    return sigma0 * Nval;
}

//----------- amplitude -----------

double amplitude(double x, double Q2, const Meson& M,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0)
{
    auto fr = [&](double r) {
        double Ov = overlap_r(r, Q2, M, Nz);
        double sigma_qq = sigma_dipolo(r, x);
        return 0.5 * r * Ov * sigma_qq; // r de d²r = 2π r dr
    };
    double Ir = integrate_simpson(fr, rmin, rmax, Nr);
    return Ir;
}


// ---------------- cálculo da seção de choque ----------------

double sigma_x(double x, double Q2 , const Meson& M,
               int Nr = 600, int Nz = 200, 
               double rmin = 1e-4, double rmax = 10.0)
{
    // A seção de choque é proporcional ao quadrado da amplitude.
    double amp = amplitude(x, Q2, M, Nr, Nz, rmin, rmax);
    double B_val = B(x, Q2,  M); // Unidades de GeV^-2
    
    return ((amp * amp) / (16.0 * M_PI * B_val));
            //;
};

void run_sigma_plot(double Q2, const Meson& M_GLC, const Meson& M_BG)
{
    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string filename = "csv/" + M_GLC.meson + "_sigma_Q2=" + Q2_str + ".csv";
    std::ofstream fout(filename);
    fout << "W,sigma_GLC,sigma_BG\n";

    const int Nw = 150;
    double Wmin = M_GLC.MV;
    double Wmax = 2e4;

    for (int i = 0; i< Nw; ++i) {
        double frac = static_cast<double>(i) / (Nw-1);
        double W = Wmin * pow(Wmax / Wmin, frac);
        double x = (M_GLC.MV * M_GLC.MV) / (W * W);
        double sigma_GLC = sigma_x(x, Q2, M_GLC) * GeV2_to_nb;
        double sigma_BG  = sigma_x(x, Q2, M_BG) * GeV2_to_nb;
        std::cout << "W = " << W << ", sigma_GLC = " << sigma_GLC << ", sigma_BG = " << sigma_BG << std::endl;
        fout << W << "," << sigma_GLC << "," << sigma_BG << "\n";
    




    //const int Nx = 40;
    //double xmin = 1e-4, xmax = 1e-1;
    //for (int i = 0; i < Nx; ++i) {
    //    double frac = static_cast<double>(i) / (Nx - 1);
    //   double x = xmin * pow(xmax / xmin, frac);
    //    double W2 = M_GLC.MV * M_GLC.MV / x;
    //    double W = sqrt(W2);

    //    double sigma_GLC = sigma_x(x, Q2, M_GLC);
    //    double sigma_BG  = sigma_x(x, Q2, M_BG);
    //   std::cout << "W = " << W << ", sigma_GLC = " << sigma_GLC * GeV2_to_nb << ", sigma_BG = " << sigma_BG * GeV2_to_nb << std::endl;
    //    fout << W << "," << sigma_GLC * GeV2_to_nb << "," << sigma_BG * GeV2_to_nb << "\n";
    }
    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado." << std::endl;
}

void perfil(const Meson& meson){
    if (meson.isGLC){
        std::cout << "Perfil do méson Gaus-LC:\n";
        std::cout << "Méson: " << meson.meson << "\n";
        std::cout << "Massa do méson (GeV): " << meson.MV << "\n";
        std::cout << "Massa do quark (GeV): " << meson.mf << "\n";
        std::cout << "Carga efetiva: " << meson.ef << "\n";
        std::cout << "Normalização transversa NT: " << meson.NT << "\n";
        std::cout << "Parâmetro transverso R2T (GeV^-2): " << meson.R2T << "\n";
        std::cout << "Normalização longitudinal NL: " << meson.NL << "\n";
        std::cout << "Parâmetro longitudinal R2L (GeV^-2): " << meson.R2L << "\n";
    } else {
        std::cout << "Perfil do méson Boosted Gaussian:\n";
        std::cout << "Méson: " << meson.meson << "\n";
        std::cout << "Massa do méson (GeV): " << meson.MV << "\n";
        std::cout << "Massa do quark (GeV): " << meson.mf << "\n";
        std::cout << "Carga efetiva: " << meson.ef << "\n";
        std::cout << "Normalização transversa NT: " << meson.NT << "\n";
        std::cout << "Normalização longitudinal NL: " << meson.NL << "\n";
        std::cout << "Parâmetro R2 (GeV^-2): " << meson.R2 << "\n";
    }
}


int main(){
    double Q2 = 0.0; // Fotoprodução

    
    //N_plot();
    //run_sigma_plot(Q2_phi, phi_GLC, phi_BG);
    run_sigma_plot(Q2, Jpsi_GLC, Jpsi_BG);
    return 0;
}
