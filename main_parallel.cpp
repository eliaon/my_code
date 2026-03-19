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
#include "wavefunctions.h"
#include "GBW.h"
#include "correcs.h"
#include "ipsat.h"

using namespace MZ_ipsat;






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



// ----------------- função de calcular a seção de choque total sigma(x) -----------------
//double sigma_x(double x, double Q2 , const Meson& M,
 //              int Nr = 600, int Nz = 200, 
  //             double rmin = 1e-4, double rmax = 10.0)
//{
    // A seção de choque é proporcional ao quadrado da amplitude.
  //  double amp = amplitude(x, 0.0, Q2, M, Nr, Nz, rmin, rmax);
  //  double B_val = B(x, Q2,  M); // Unidades de GeV^-2
  //  double lambda_e = calculate_lambda(x, Q2, amp, M);
  //  double Rg_val = RG(x, Q2, lambda_e, M);
  //  double beta_val = beta(x, Q2, lambda_e, M);

   // return ((amp * amp) / (16.0 * M_PI * B_val))
    //        * Rg_val * Rg_val * 
      //      (1.0+beta_val*beta_val);
//};

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

    double sigma_plus  = sigma_x_GBW(x_plus,  Q2, M);
    double sigma_minus = sigma_x_GBW(x_minus, Q2, M);

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
    run_sigma_csv_GBW();
    //plot_sigma(std::string("phi"));    //debug_correc();
    //plot_overlap();
    //printf("Qs = %g\n", QS_bCGC(1e-4,0));
    //plot_dsigma_dt(std::string("phi"));
    //dsigma_dump();
    return 0;
}