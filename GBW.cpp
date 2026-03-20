#include "utils.h"
#include "ctes.h"
#include <cmath>
#include <iostream>
#include "GBW.h"
#include "fstream"
#include <sstream>
#include "wavefunctions.h"
#include "integration.hpp"
#include "correcs.h"
#include <format>
#include <chrono>
#include "plot.h"

parametros_GBW gbw(23/mb_to_gev, 3e-4, 0.29);


// ----------- parametros GBW NEW 2018 Saturation model of DIS: an update https://doi.org/10.1007/JHEP03(2018)102
// sequência: sigma_0 , x_0, lambda
parametros_GBW gbw_5(28.18/mb_to_gev, 0.31e-4, 0.237);

parametros_GBW gbw_10(27.43/mb_to_gev, 0.40e-4, 0.248);

parametros_GBW gbw_20(26.60/mb_to_gev, 0.53e-4, 0.259);

parametros_GBW gbw_50(25.21/mb_to_gev, 0.80e-4, 0.281);


double QS2_GBW( double x, parametros_GBW params)
{
    return pow(params.x0 / x, params.lambda);
}   

double N_GBW(double r,  double x, parametros_GBW params)
{
    double Qs2 = QS2_GBW(x, params);
    double arg = (r * r) * Qs2 / 4.0;
    return (1.0 - exp(-arg))*std::pow(1.0-x, 5.26);
}

double sigma_qq_GBW(double r, double x, parametros_GBW params)
{
    return params.sigma0 * N_GBW(r, x, params);
}

void dipolo_csv(double x)
{
    std::ofstream fout("csv/sigma_qq_GBW_x=" + std::to_string(x) + ".csv");
    fout << "r,sigma_qq\n";

    const int Npoints = 5000;
    double rmin = 1e-4, rmax = 10.0;

    for (int i = 0; i < Npoints; ++i)
    {
        double frac = (double)i / (Npoints - 1);
        double r = rmin * std::pow(rmax / rmin, frac); // escala log
        double sigma = sigma_qq_GBW(r, x, gbw);
        fout << r * CFAC << "," << sigma * mb_to_gev << "\n"; // converte r para fm e sigma para mb
    }
}

double amplitude_GBW(double x, double Q2, const Meson& M,
                const parametros_GBW& gbw,
                 int Nr , int Nz,
                 double rmin , double rmax)
{
    auto amp_r = [x, Q2, M, gbw, Nz](double r) {
        double Ov = overlap_r(r, Q2, M, Nz);
        double sigma_qq = sigma_qq_GBW(r, x, gbw);
        //double sqrt_fc = std::sqrt(f_c(r));
        return 0.5 * r * Ov * sigma_qq; //* sqrt_fc; // r de d²r = 2π r dr
    };
    double amp = integrate_simpson(amp_r, rmin, rmax, Nr);
    return amp;
}

double sigma_x_GBW(double x, double Q2 , const Meson& M,
               int Nr, int Nz, 
               double rmin, double rmax)
{
    parametros_GBW params = gbw_50; // ou escolha outro conjunto de parâmetros se desejar
    double amp = amplitude_GBW(x, Q2, M, params, Nr, Nz, rmin, rmax);
    double B_val = B(x, Q2, M);
    double lambda_e = calculate_lambda(x, Q2, amp, M);
    double RG_val = RG(x, Q2, lambda_e, M);
    double beta_val = beta(x, Q2, lambda_e, M);
    double correction_factor = RG_val * RG_val * (1.0 + beta_val * beta_val);
    return correction_factor * (amp * amp) / (16.0 * M_PI * B_val);
}

void run_sigma_csv_GBW()
{
    double Q2;
    std::cout << "Insira o valor de Q2 (GeV^2): ";
    std::cin >> Q2;

    const Meson& M_GLC = input_meson(); //importa o glc do meson escolhido
    const Meson& M_BG  = meson_models.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string filename = "csv/" + M_GLC.meson + "_sigma_GBW(50)_Q2=" + Q2_str + ".csv";
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
                                                                                
        double frac = static_cast<double>(i) / (Nw - 1);
        double W = Wmin * pow(Wmax / Wmin, frac);
        double x = (M_GLC.MV * M_GLC.MV) / (W * W);

        Wv[i] = W;
        sigma_GLC_v[i] = sigma_x_GBW(x, Q2, M_GLC) * GeV2_to_nb;
        sigma_BG_v[i]  = sigma_x_GBW(x, Q2, M_BG) * GeV2_to_nb;
        std::cout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "%\n";

    }

    for (int i = 0; i < Nw; ++i) {
    fout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
}

    fout.close();
    auto end = clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução: " << duration << " ms" << std::endl;
    std::cout << "Arquivo '" << filename << "' gerado em " << duration << " ms" << std::endl;

    plot_sigma(M_GLC.meson, filename);
}
