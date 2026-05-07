#include <iostream>
#include <fstream>
#include <chrono>

#include "ctes.h"
#include "utils.h"
#include "plot.h"
#include "correcs.h"
#include "IIM.hpp"
#include "wavefunctions.h"
#include "integration.hpp"


param_IIM IIM_S("IIM_S", 26.25 * mb_to_gev2, 0.6194, 0.2131e-4, 0.2545, 9.9, 0.7);
param_IIM IIM_RS("IIM_RS", 21.85 * mb_to_gev2, 0.762, 6.226e-5, 0.2319, 9.9, 0.7);


namespace IIM{

double N(double r,  double x, const param_IIM& param)
{

    double Qs  = std::pow(param.x_0/x, param.lambda/2.0);
    double rQs = r * Qs;

    rQs = std::max(rQs, 1e-12);      // proteção numérica

    if (rQs < 2.0) {

        double Lx = std::log(1.0/x);

        double exp =
            2.0 * (param.gamma_s +
            std::log(2.0/rQs) /
            (param.kappa * param.lambda * Lx));

        return param.N_0 * std::pow(rQs/2.0, exp)*std::pow(1.0-x, 6);

    } else {

        double a = -param.N_0*param.N_0*param.gamma_s*param.gamma_s /
                   ((1-param.N_0)*(1-param.N_0)*std::log(1.0-param.N_0));

        double b = 0.5 * std::pow(1.0-param.N_0,
                   -(1.0-param.N_0)/(param.N_0*param.gamma_s));

        double ln = std::log(b * rQs);

        return 1.0 - std::exp(-a * ln * ln)*std::pow(1.0-x, 6);
    }
}

double sigma_qq(double r, double x, const param_IIM& param){return param.sigma0 * N(r, x, param);}

double amplitude_p(double x, double Q2, const Meson& M)
{
    param_IIM param = IIM_S; // escolha do parâmetro IIM (pode ser IIM_S ou IIM_RS)
    auto amp_r = [x, Q2, &M, param](double r) {
        double Ov = sigma_qq(r, x, param) * overlap_r(r, Q2, M);
        return r * Ov;
    };
    double amp = 0.5 * integrate_simpson(amp_r, 1e-4, 10.0, 300);
    return amp;
}

double sigma_p(double x, double Q2, const Meson& M)
{
    param_IIM param = IIM_S; // escolha do parâmetro IIM (pode ser IIM_S ou IIM_RS)
    double amp = amplitude_p(x, Q2, M);
    double slope = B_slope(x, Q2, M);
    double lambda_e = calculate_lambda(x, Q2, M, "IIM");
    double R_g = RG(x, Q2, lambda_e, M);
    double beta_corr = beta(x, Q2, lambda_e, M);
    double A2 = (amp * amp) / (16.0 * M_PI * slope);
    double correction_factor = R_g * R_g * (1.0 + beta_corr * beta_corr);
    return  A2*correction_factor; // B=4.0 GeV^-2 é um valor típico para o slope
}

void N_csv()
{
    vector<double> x_values = {1e-4, 1e-3, 1e-2};
    vector<string> filenames;
    
    for (double x : x_values) {
        std::string filename = "csv/N_IIM_x=" + doubleParaString(x) + ".csv";
    std::ofstream fout(filename);
    fout << "r,N\n";

    const int Npoints = 5000;
    double rmin = 1e-4, rmax = 10.0;

    for (int i = 0; i < Npoints; ++i)
    {
        double frac = (double)i / (Npoints - 1);
        double r = rmin * std::pow(rmax / rmin, frac)*CFAC; // escala log
        double N_val = N(r, x, IIM_S); // exemplo para x=1e-4 e x0=1e-2
        fout << r / CFAC << "," << N_val << "\n"; // converte r para fm
    }
    filenames.push_back(filename);
    fout.close();
    }
    plot_N_multi(filenames, x_values, "IIM");

}

void sigma_p_csv()
{
    double Q2 = 0.0;
    const Meson& M_GLC = input_meson("GBW");
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG;

    std::string filename = "csv/" + M_GLC.meson + "_sigma_gammap_IIM.csv";
    std::ofstream fout(filename);

    double Wmin = x_to_W(0.01, M_GLC);
    double Wmax = 3e+3;

    const int Nw = 100;

    std::vector<double> W_vals(Nw), sigma_GLC_vals(Nw), sigma_BG_vals(Nw);

    using clock = std::chrono::steady_clock;
    auto start = clock::now();

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nw; ++i) {
        double frac = static_cast<double>(i) / (Nw - 1);
        double W = Wmin * pow(Wmax / Wmin, frac);
        double x = (M_GLC.MV * M_GLC.MV) / (W * W);

        double s_GLC = sigma_p(x, Q2, M_GLC);
        double s_BG  = sigma_p(x, Q2, M_BG);

        W_vals[i]         = W;
        sigma_GLC_vals[i] = s_GLC * GeV2_to_nb;
        sigma_BG_vals[i]  = s_BG  * GeV2_to_nb;
        cout << "W = " << W << " GeV, sigma_GLC = " << sigma_GLC_vals[i] << " nb, sigma_BG = " << sigma_BG_vals[i] << " nb\n";
    }

    auto end = clock::now();

    fout << "W,sigma_GLC,sigma_BG\n";

    for (int i = 0; i < Nw; ++i) {
        fout << W_vals[i] << "," 
             << sigma_GLC_vals[i] << "," 
             << sigma_BG_vals[i] << "\n";
    }

    fout.close();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução: " << duration << " ms\n";
    std::cout << "Arquivo '" << filename << "' gerado com sucesso.\n";

    plot_sigma(M_GLC.meson, filename, "IIM");
}










}