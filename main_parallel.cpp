#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <format>
#include <boost/math/special_functions/bessel.hpp>
#include "integration.hpp"
#include "dipoleamplitude.hpp"
#include "dglap_cpp/AlphaStrong.h"
#include "dglap_cpp/EvolutionLO_nocoupling.h"
#include <chrono>
#include <omp.h>
#include "plot.h"
#include "utils.h"
#include "ctes.h"

using namespace MZ_ipsat;


// ---------------- Classe de parâmetros GBW --------------
class parametros_GBW {
    public:
        double sigma0; // GeV^-2
        double x0;
        double lambda;

        parametros_GBW(double sigma0_, double x0_, double lambda_)
            : sigma0(sigma0_), x0(x0_), lambda(lambda_) {}
};

parametros_GBW gbw_kowalski(23/mb_to_gev, 3e-4, 0.29);


// ---------------- Classe Meson ----------------





// ---------------- slope B(Q2) ----------------
inline double B(double x, double Q2, const Meson& M) {
    double W = std::sqrt(M.MV*M.MV/x);
    if (M.meson == "Jpsi"){
        double B1 = 4.80 + 4.0* 0.133 *log(W/90.0); //valores do lhcb dados pelo haimon xdxd
        return B1;
    }
        else if (M.meson == "phi"){
            double B2 = 0.55 * (14.0 / pow((Q2 + M.MV*M.MV), 0.2) + 1.0);
            std::cout << "B = " << B2 << "para x = " << x << "%\n";
            return B2;}
            else {
                std::cerr << "Méson desconhecido para cálculo de B: " << M.meson << std::endl;
                return 0.0;
            }
}

void print_B_values(const Meson& M, double Q2, double Wmin=30.0, double Wmax=10000.0, int Npoints = 30)
{
    std::cout << "Valores de B para o méson " << M.meson << " com Q2 = " << Q2 << " GeV²:\n";
    for (int i = 0; i < Npoints; ++i) {
        double frac = static_cast<double>(i) / (Npoints - 1);
        double W = Wmin * pow(Wmax / Wmin, frac);
        double x = (M.MV * M.MV) / (W * W);
        double B_val = B(x, Q2, M);
        std::cout << "W: " << W << " GeV, B: " << B_val << " GeV^-2\n";
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

inline double laplacian_phi_L(double r, double z, const Meson& M)
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
inline double psi_Vpsi_L(double z, double r, double Q2, const Meson& M)
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
double overlap_r(double r, double Q2, const Meson& M, int Nz = 200) {
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

// ------------ calculo do N -----------
 
 double QS2_GBW( double x, parametros_GBW params)
{
    return pow(params.x0 / x, params.lambda);
}   
 
double QS_bCGC( double x, double b, double x0 = 1.84e-6)
{
    double B_CGC = 7.5; // GeV^-2
    double Qs2 = std::pow(x0/x, lambda/2.0)*
                std::pow(std::exp(-b*b/(2*B_CGC)), 1.0/(2.0*gamma_s));
    return Qs2;
}
 
double N_GBW(double r,  double x, parametros_GBW params)
{
    double Qs2 = QS2_GBW(x, params);
    double arg = (r * r) * Qs2 / 4.0;
    return (1.0 - exp(-arg))*std::pow(1.0-x, 5.26);
}
 double N_IIM(double r,  double x,  double x0 =  1.84e-6)
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

double prof_bCGC(double r, double x, double b, double x0 = 1.84e-6)
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

double N_bCGC(double r, double x, double x_0 = 1.84e-6)
{
    auto integrand = [&](double b) {
        return b * prof_bCGC(r, x, b, x_0);
    };
    double bmax = 10.0; // Limite superior para a integração em b
    return 4.0 * M_PI * integrate_simpson(integrand, 0.0, bmax, 200);
}

void dump_curve_N(const std::string& fname, double x)
{
    DipoleAmplitude dipole(MZ_IPSAT);
    dipole.EnableLookupTable();
    std::ofstream fout(fname);
    const int Npoints = 5000;

    fout << "r,N_ipsat\n";

    double rmin=1e-4, rmax=10.0;

    for(int i=0;i<Npoints;++i)
    {
        double frac = (double)i/(Npoints-1);
        double r = rmin + frac*(rmax-rmin);

        double Nval = dipole.N(r, x, 0.0);
        fout << r/CFAC << "," << Nval << "\n";
    }
}

void N_plot(void){
    for(double x : {1e-4, 1e-3, 1e-2})
    {
        std::string fname = "csv/N_ipsat_x=" + std::to_string(x) + ".csv";
        dump_curve_N(fname, x);
        std::cout << "Arquivo '" << fname << "' gerado." << std::endl;
    }
}

double sigma_dipolo_ipsat(double r, double x, double Delta, DipoleAmplitude& dipole)
{
    auto N_b = [&](double b) {
        return 2.0 * M_PI * b * dipole.N(r, x, b) * boost::math::cyl_bessel_j(0, Delta * b);
    };
    double bmax = 10.0;
    double Nval = integrate_simpson(N_b, 0.0, bmax, 200);
    return 2 * Nval;
}

void csv_sigmaqq(void)
{
    for (double x : {1e-4, 1e-2})
    {
        std::ofstream fout("csv/sigma_dipolo_ipsat_x=" + std::to_string(x) + ".csv");
        fout << "r,sigma_dipolo\n";

        const int Npoints = 5000;
        double rmin = 1e-4, rmax = 10.0;

        std::vector<std::pair<double,double>> results(Npoints);

        #pragma omp parallel
        {
            DipoleAmplitude dipole(MZ_IPSAT);
            dipole.EnableLookupTable();
            for (int i = 0; i < Npoints; ++i)
            {
                double frac = (double)i / (Npoints - 1);
                double r    = rmin * std::pow(rmax / rmin, frac); // escala log
                double sigma = sigma_dipolo_ipsat(r, x, 0.0, dipole);
                results[i]  = {r * 0.1973, sigma * 0.3894}; // converte r para fm e sigma para mb
            }
        }

        for (auto& [r, s] : results)
            fout << r << "," << s << "\n";
    }
}

//----------- amplitude -----------
 
double amplitude( double x, double Delta, double Q2, const Meson& M,
                 int Nr = 600, int Nz = 200,
                 double rmin = 1e-4, double rmax = 10.0)
{
    DipoleAmplitude dipole(MZ_IPSAT);
    dipole.EnableLookupTable();
    auto fr = [&](double r) {
        double Ov = overlap_r(r, Q2, M, Nz);
        double sigma_qq = sigma_dipolo_ipsat(r, x, Delta, dipole);
        return 0.5 * r * Ov * sigma_qq; // r de d²r = 2π r dr
    };
    double Ir = integrate_simpson(fr, rmin, rmax, Nr);
    return Ir;
}

// Fatores de correção Rg e deltinha
 
 double lnA( double y, double Q2, const Meson& M)
{
     double x = std::exp(-y);
    double amp = amplitude(x, 0.0, Q2, M);

    if (amp <= 0) {
        std::cerr << "Atenção: amplitude não positiva em y=" << y << std::endl;
        return 0.0;
    }

    return std::log(amp);
}

// -------------- lamnda_e --------------

double calculate_lambda(double x, double Q2, const Meson& M)
{
    const double h = 1e-4;
    double err;
    double y = -std::log(x);

    if (x > 1e-2) {
        return 0.2; // valor aproximado para x muito grande
    }
    else{
         auto f_lnA = [&](double y) {
        return lnA(y, Q2, M);
    };

    double dlnA_dy = dfridr(f_lnA, y, h, err);
    return dlnA_dy;   
    }
}

// ----------------- fator de correção Rg skeweness ----------------

double calculate_RG(double x, double Q2, double lambda_e, const Meson& M)
{
    return std::pow(2.0,2.0*lambda_e +3.0) * tgamma(lambda_e + 2.5) 
                         /(sqrt(M_PI) * tgamma(lambda_e + 4.0));
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
        double lambda_e = calculate_lambda(xi, Q2, M);
        double Rg = calculate_RG(xi, Q2, lambda_e, M);

        std::cout << "x: " << xi << "  lambda_e: " << lambda_e << "  Rg: " << Rg << std::endl;
    }
}


// ---------------- cálculo da seção de choque ----------------
// ---------------- primeiro em ipsat, depois com os fatores de correção ----------------
double dsigma_dt(double x, double Q2, const Meson& M, double t,
                 int Nr = 600, int Nz = 200, 
                 double rmin = 1e-4, double rmax = 10.0)
{
    double Delta = std::sqrt(t);
    double amp = amplitude(x, Delta, Q2, M, Nr, Nz, rmin, rmax);
    double lambda_e = calculate_lambda(x, Q2, M);
    double Rg = calculate_RG(x, Q2, lambda_e, M);
    double beta_val = beta(x, Q2, lambda_e, M);
    return ((amp * amp) / (16.0 * M_PI )) * Rg * Rg * (1.0 + beta_val * beta_val);
}
// ----------------- seção de choque total sigma(x) a partir de dsigma/dt em t=0 e B -----------------
// ----------------- mas para usar a integração em t, basta descomentar o código e comentar a linha de baixo -----------------    
double sigma_ipsat(double x, double Q2, const Meson& M,
                 int Nr = 600, int Nz = 200, 
                 double rmin = 1e-4, double rmax = 10.0)
{
    double B_val = B(x, Q2, M);
    double sigma = dsigma_dt(x, Q2, M, 0)/(B_val);
    //double N_t = 100;
    //auto ft = [&](double t) {
    //    return dsigma_dt(x, Q2, M, t, Nr, Nz, rmin, rmax);
    //};
    //double sigma = integrate_simpson(ft, 0.0, 2.5, N_t);
    return sigma;
}


// ----------------- função para gerar csv da seção de choque dsigma/dt para um W específico ---------- 
void dsigma_dt_csv(double W, const Meson& M_GLC)
{
    double x = W_to_x(W, M_GLC);

    const Meson& M_BG  = meson_models.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    double Q2 = 0.0;
    string W_str = doubleParaString(W, 0);
    std::string filename ="csv/" + M_GLC.meson + "_dsigma_dt_W=" + W_str + "GeV.csv";
    std::ofstream fout(filename);
    fout << "t,dsigma_dt_GLC,dsigma_dt_BG\n";

    const int Npoints = 150;
    double tmin = 0.0, tmax = 4.0; // GeV^2

    //esse for escolhe os valores de t e calcula a seção de choque para cada modelo, salvando no arquivo csv
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Npoints; ++i) {
        double frac = static_cast<double>(i) / (Npoints - 1);
        double t = tmin + frac * (tmax - tmin);

        double dsdt_glc = dsigma_dt(x, Q2, M_GLC, t) * GeV2_to_nb; // converte para nb/GeV^2
        double dsdt_bg  = dsigma_dt(x, Q2, M_BG, t) * GeV2_to_nb;  // converte para nb/GeV^2
        cout << t << "," << dsdt_glc << "," << dsdt_bg << "\n";
        fout << t << "," << dsdt_glc << "," << dsdt_bg << "\n";
    }
    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado." << std::endl;

    plot_dsigma_dt(M_GLC.meson);
}
// ----------------- função para gerar os csvs de dsigma/dt para os W escolhidos ---------- 
void dsigma_dump(void)
{
    const Meson M_GLC = input_meson(); //importa o glc do meson escolhido

    for(double W : {70.0}){
        dsigma_dt_csv(W, M_GLC);
    }
  
}
// ----------------- função de calcular a seção de choque total sigma(x) -----------------
double sigma_x(double x, double Q2 , const Meson& M,
               int Nr = 600, int Nz = 200, 
               double rmin = 1e-4, double rmax = 10.0)
{
    // A seção de choque é proporcional ao quadrado da amplitude.
    double amp = amplitude(x, 0.0, Q2, M, Nr, Nz, rmin, rmax);
    double B_val = B(x, Q2,  M); // Unidades de GeV^-2
    double lambda_e = calculate_lambda(x, Q2, M);
    double Rg = calculate_RG(x, Q2, lambda_e, M);
    double beta_val = beta(x, Q2, lambda_e, M);

    return ((amp * amp) / (16.0 * M_PI * B_val))
            * Rg * Rg * 
            (1.0+beta_val*beta_val);
};

// ---------------- distribuição de rapidez ----------------

double fluxo_fotons(double omega, double sqrt_s)
{
    if (omega <= 0.0 || omega >= sqrt_s / 2.0)
        return 0.0;

    double gamma_L = sqrt_s / (2.0 * 0.938);
    double Q2min = std::pow(omega / gamma_L, 2);

    if (Q2min <= 0.0) return 0.0;

    double Ohm = 1.0 + 0.71 / Q2min;

    double pref = alfem / (2.0 * M_PI * omega);
    double kin  = 1.0 + std::pow(1.0 - 2.0 * omega / sqrt_s, 2);
    double logt = std::log(Ohm) - 11.0/6.0
                  + 3.0/Ohm - 3.0/(2.0*Ohm*Ohm)
                  + 1.0/(3.0*Ohm*Ohm*Ohm);

    return pref * kin * logt;
}


double d_sigma_dy(double y, double sqrt_s, double Q2, const Meson& M)
{
    double MV = M.MV;

    double omega_plus  = (MV / 2.0) * std::exp(+y);
    double omega_minus = (MV / 2.0) * std::exp(-y);

    double W_plus  = std::sqrt(2.0 * omega_plus  * sqrt_s);
    double W_minus = std::sqrt(2.0 * omega_minus * sqrt_s);

    double x_plus  = (MV * MV) / (W_plus * W_plus);
    double x_minus = (MV * MV) / (W_minus * W_minus);

    if (x_plus <= 0.0 || x_plus >= 1.0)  return 0.0;
    if (x_minus <= 0.0 || x_minus >= 1.0) return 0.0;

    double n_plus  = omega_plus  * fluxo_fotons(omega_plus,  sqrt_s);
    double n_minus = omega_minus * fluxo_fotons(omega_minus, sqrt_s);

    double sigma_plus  = sigma_x(x_plus,  Q2, M);
    double sigma_minus = sigma_x(x_minus, Q2, M);

    return n_plus * sigma_plus + n_minus * sigma_minus;
}


void dump_curve_rap(const std::string& meson_input, double sqrt_s)
{
    
    double Q2=0.0;

    auto it = meson_models.find(meson_input);
    if (it == meson_models.end()) {
        std::cerr << "Meson invalido. Usando Jpsi por padrao.\n";
        it = meson_models.find("Jpsi");
    }

    const Meson& M_GLC = it->second.M_GLC;
    const Meson& M_BG  = it->second.M_BG;

    std::string sqrt_s_str = std::format("{:.3g}", sqrt_s);
    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string filename = "csv/" + M_GLC.meson + "_rapidez_" + sqrt_s_str + "GeV.csv";
    std::ofstream fout(filename);
    fout << "y,d_sigma_dy_GLC,d_sigma_dy_BG\n";

    const int Ny = 100;
    double ymax = std::log(sqrt_s / M_GLC.MV);
    #pragma omp parallel for schedule(dynamic)
    for (int i = -Ny; i <= Ny; ++i) {
        double y = (static_cast<double>(i) / Ny) * ymax;
        double dsdy_GLC = d_sigma_dy(y, sqrt_s, Q2, M_GLC) * GeV2_to_nb;
        double dsdy_BG  = d_sigma_dy(y, sqrt_s, Q2, M_BG) * GeV2_to_nb;
        fout << y << "," << dsdy_GLC << "," << dsdy_BG << "\n";
    }
    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado." << std::endl;
}

void run_rapidez_csv(void){

    double sqrt_s;
    std::cout << "Insira o valor de sqrt(s) (GeV): ";
    std::cin >> sqrt_s;


    std::string Jpsi="Jpsi";
    std::string phi = "phi";
    using clock = std::chrono::steady_clock;
    auto start = clock::now();
    dump_curve_rap(Jpsi, sqrt_s);
    auto end = clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução para Jpsi: " << duration << " ms" << std::endl;

     start = clock::now();
    dump_curve_rap(phi, sqrt_s);
    end = clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução para phi: " << duration << " ms" << std::endl;
}


void run_sigma_csv()
{
    double Q2;
    std::cout << "Insira o valor de Q2 (GeV^2): ";
    std::cin >> Q2;

    const Meson& M_GLC = input_meson(); //importa o glc do meson escolhido
    const Meson& M_BG  = meson_models.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string filename = "csv/" + M_GLC.meson + "_sigma_Q2=" + Q2_str + ".csv";
    std::ofstream fout(filename);
    fout << "W,sigma_GLC,sigma_BG\n";

    const int Nw = 100;
    double Wmin = 25.0;
    double Wmax = 3e3;

    using clock = std::chrono::steady_clock;
    auto start = clock::now();

    std::vector<double> Wv(Nw), sigma_GLC_v(Nw), sigma_BG_v(Nw);


    std::cout << "W (Gev), sigma_GLC (nb), sigma_BG (nb)\n";
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Nw; ++i) {
        double perc = (double)i / (Nw - 1) * 100.0;
        double frac = static_cast<double>(i) / (Nw - 1);
        double W = Wmin * pow(Wmax / Wmin, frac);
        double x = (M_GLC.MV * M_GLC.MV) / (W * W);

        Wv[i] = W;
        sigma_GLC_v[i] = sigma_ipsat(x, Q2, M_GLC) * GeV2_to_nb;
        sigma_BG_v[i]  = sigma_ipsat(x, Q2, M_BG) * GeV2_to_nb;
        std::cout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << ", " << perc << "%\n";

    }

    for (int i = 0; i < Nw; ++i) {
    fout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
}

    fout.close();
    auto end = clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução: " << duration << " ms" << std::endl;
    std::cout << "Arquivo '" << filename << "' gerado em " << duration << " ms" << std::endl;

    plot_sigma(M_GLC.meson);
}





int main(){
    int threads;
    std::cout << "Insira o número de threads para OpenMP: ";
    std::cin >> threads;
    omp_set_num_threads(threads);// Ajuste conforme o número de núcleos disponíveis
    //draw_wavefunctions(Jpsi_BG);
    //print_B_values(Jpsi_GLC, 0.0);
    //N_plot();
    //plot_sigmaqq();
    //run_rapidez_plot();
    //run_sigma_csv();
    //plot_sigma(std::string("phi"));    //debug_correc();
    //plot_overlap();
    //printf("Qs = %g\n", QS_bCGC(1e-4,0));
    //plot_dsigma_dt(std::string("phi"));
    dsigma_dump();
    return 0;
}