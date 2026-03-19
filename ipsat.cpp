#include "integration.hpp"
#include "dipoleamplitude.hpp"
#include "dglap_cpp/AlphaStrong.h"
#include "dglap_cpp/EvolutionLO_nocoupling.h"
#include <chrono>
#include <omp.h>
#include "plot.h"
#include "utils.h"
#include "ctes.h"
#include "wavefunctions.h"
#include "correcs.h"
#include "boost/math/special_functions/bessel.hpp"
#include <fstream>
#include <sstream>

using namespace MZ_ipsat;


// ----------------- função de calcular a seção de choque de dipolo -----------------
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

void N_ipsat_csv(void){
    for(double x : {1e-4, 1e-3, 1e-2})
    {
        std::string fname = "csv/N_ipsat_x=" + std::to_string(x) + ".csv";
        dump_curve_N(fname, x);
        std::cout << "Arquivo '" << fname << "' gerado." << std::endl;
    }
}



//----------- amplitude -----------
double amplitude_ipsat(double x, double Delta, double Q2, const Meson& M,
                 int Nr, int Nz,
                 double rmin, double rmax)
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

// ---------------- cálculo da seção de choque ----------------
// ---------------- primeiro em ipsat, depois com os fatores de correção ----------------
double dsigma_dt(double x, double Q2, const Meson& M, double t,
                 int Nr, int Nz, 
                 double rmin, double rmax)
{
    double Delta = std::sqrt(t);
    double amp = amplitude_ipsat(x, Delta, Q2, M, Nr, Nz, rmin, rmax);
    double lambda_e = calculate_lambda(x, Q2, amp, M);
    double Rg_val = RG(x, Q2, lambda_e, M);
    double beta_val = beta(x, Q2, lambda_e, M);
    return ((amp * amp) / (16.0 * M_PI )) * Rg_val * Rg_val * (1.0 + beta_val * beta_val);
}

// ----------------- seção de choque total sigma(x) a partir de dsigma/dt em t=0 e B -----------------

double sigma_ipsat_slope(double x, double Q2, const Meson& M,
                 int Nr, int Nz, 
                 double rmin, double rmax)
{
    double B_val = B(x, Q2, M);
    double sigma = dsigma_dt(x, Q2, M, 0, Nr, Nz, rmin, rmax)/(B_val);
    return sigma;
}

// ----------------- integrando em t -----------------    
double sigma_ipsat_integrado(double x, double Q2, const Meson& M,
                 int Nr, int Nz, 
                 double rmin, double rmax)
{
    double N_t = 100;
    auto ft = [&](double t) {
        return dsigma_dt(x, Q2, M, t, Nr, Nz, rmin, rmax);
    };
    double sigma = integrate_simpson(ft, 0.0, 2.5, N_t);
    return sigma;
}


// ----------------- função para gerar csv da seção de choque dsigma/dt para um W específico ---------- 
void dsigma_dt_csv(double W, const Meson& M_GLC)
{
    double x = W_to_x(W, M_GLC);

    const Meson& M_BG  = meson_models.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    double Q2 = 0.0;
    std::string W_str = doubleParaString(W, 0);
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

        double dsdt_glc = dsigma_dt(x, Q2, M_GLC, t,  600, 200, 1e-4, 10.0) * GeV2_to_nb; // converte para nb/GeV^2
        double dsdt_bg  = dsigma_dt(x, Q2, M_BG, t, 600, 200, 1e-4, 10.0) * GeV2_to_nb;  // converte para nb/GeV^2
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



