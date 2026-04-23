


#include <chrono>
#include <omp.h>

#include "plot.h"
#include "utils.h"
#include "ctes.h"
#include "wavefunctions.h"
#include "fstream"
#include "correcs.h"
#include "integration.hpp"

#include "dipoleamplitude.hpp"
#include "dglap_cpp/AlphaStrong.h"
#include "dglap_cpp/EvolutionLO_nocoupling.h"

using namespace MZ_ipsat;


namespace DGLAP{
// --------------------- DIPOLO DGLAP ---------------------

double N_p(double r, double x, DipoleAmplitude& dipole)
{
    double musqrd = dipole.C/(r*r) + dipole.GetMu0() * dipole.GetMu0();
    

    double alpha_S = dipole.Alphas(std::sqrt(musqrd));
    double xg = dipole.xg(x, musqrd);

    double alphasxg = alpha_S * xg;

    double arg = std::pow(M_PI, 2) * r * r * alphasxg / (3.0 * sigma0_dglap);

    cout <<"x=" << x << " | mu^2: " << musqrd << " | alphasxg: " << alphasxg << " | arg: " << arg << std::endl;
    return 1.0 - exp(-arg);
}

double amplitude(double x, double Q2, const Meson& M,
                      DipoleAmplitude& dipole,
                      int Nr = 200, double rmin = 1e-4, double rmax = 10.0)
{
    auto amp_r = [x, Q2, &M, &dipole](double r) {

        //double fc = f_c(r,-0.979599, 0.403569, 6.8); // B=0.979599, omega=0.403569, R=6.8 fm
        //double sqrt_fc = std::sqrt(fc); // fator de correção da função de onda do fóton
        double Ov = sigma0_dglap * overlap_r(r, Q2, M) * N_p(r, x, dipole);


        return  r * Ov;
    };
    double amp = 0.5 * integrate_simpson(amp_r, rmin, rmax, Nr);
    return amp;
}



double sigma_p_DGLAP(double x, double Q2, const Meson& M, DipoleAmplitude& dipole)
{
    double amp = amplitude(x, Q2, M, dipole);
    

    double lambda_e = calculate_lambda(x, Q2, M, "DGLAP");
    double Rg = RG(x, Q2, lambda_e, M);

    double beta_corr = beta(x, Q2, lambda_e, M);
    double B_val = B_slope(x, Q2, M);

    std::cout << "x = " << x << " | Amplitude DGLAP: " << amp << std::endl;
    std::cout << "lambda_e: " << lambda_e << " | Rg: " << Rg << " | beta_corr: " << beta_corr << std::endl;
    std::cout << "B_slope: " << B_val << std::endl;
    double correction_factor = Rg * Rg * (1.0 + beta_corr * beta_corr);
    return  correction_factor * (amp * amp) / (16.0 * M_PI * B_val);
}

void sigma_p_csv()
{
    DipoleAmplitude dipole(MZ_IPSAT);
    dipole.EnableLookupTable();
    const Meson& M_GLC = input_meson("ipsat");
    const Meson& M_BG = meson_modelsipsat.find(M_GLC.meson)->second.M_BG;


    std::string filename = "csv/" + M_GLC.meson + "_sigma_DGLAP.csv";
    std::ofstream fout(filename);

    double Wmin = x_to_W(0.01, M_GLC);
    double Wmax = 3e+3;

    const int Nw = 10;
    // Criamos vetores para armazenar os resultados temporariamente
    std::vector<double> W_vals(Nw), sigma_GLC_vals(Nw), sigma_BG_vals(Nw);

    using clock = std::chrono::steady_clock;
    auto start = clock::now();

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Nw; ++i) {
        double frac = static_cast<double>(i) / (Nw - 1);
        double W = Wmin * pow(Wmax / Wmin, frac);
        double x = (M_GLC.MV * M_GLC.MV) / (W * W);
        double Q2 = 0.0; // Para fotoprodução, Q² é zero

        
        double s_GLC = sigma_p_DGLAP(x, Q2, M_GLC, dipole);
        double s_BG = sigma_p_DGLAP(x, Q2, M_BG, dipole);

        W_vals[i] = W;
        sigma_GLC_vals[i] = s_GLC * GeV2_to_nb;
        sigma_BG_vals[i] = s_BG * GeV2_to_nb;
        // Opcional: imprimir progresso parcial
        double i_float = static_cast<double>(i + 1) / Nw * 100;
        std::cout << "Progresso: " << i_float << "%" << std::endl;
    }

    auto end = clock::now();

    fout << "W,sigma_GLC,sigma_BG\n"; // Cabeçalho (ajuste se preferir outro)
    for (int i = 0; i < Nw; ++i) {
        fout << W_vals[i] << "," << sigma_GLC_vals[i] << "," << sigma_BG_vals[i] << "\n";
        
        // Opcional: imprimir no console para acompanhar o progresso final
        std::cout << "W: " << W_vals[i] << " | GLC: " << sigma_GLC_vals[i] << " | BG: " << sigma_BG_vals[i] << std::endl;
    }
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução: " << duration << " ms" << std::endl;

    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado com sucesso.\n";

    plot_sigma(M_GLC.meson, filename);

}

void N_p_csv()
{
    DipoleAmplitude dipole(MZ_IPSAT);
    dipole.EnableLookupTable();

    std::vector<double> x_vals = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6};

    double rmin = 1e-4;
    double rmax = 10.0;
    int Nr = 300;

    // agora já no formato esperado pelo plot
    std::vector<std::pair<std::string,std::string>> files_for_plot;

    for (const auto& x_val : x_vals)
    {
        std::string x_str = std::to_string(x_val);
        std::string filename = "csv/N_p_DGLAP_x=" + x_str + ".csv";

        std::ofstream fout(filename);
        fout << "r2,N_p_DGLAP\n";

        for (int i = 0; i < Nr; ++i) 
        {
            double frac = static_cast<double>(i) / (Nr - 1);
            double r = rmin * pow(rmax / rmin, frac);

            double N_val = N_p(r, x_val, dipole);

            fout << r*r << "," << N_val << "\n";
        }

        fout.close();
        std::cout << "Arquivo '" << filename << "' gerado com sucesso.\n";

        // label mais informativo
        std::string label = "DGLAP x=" + x_str;

        files_for_plot.push_back({filename, label});
    }

    // plota todos juntos
    plot_N_models(files_for_plot, x_vals[0]); // plota usando o primeiro x só pra passar o valor para o título
}

string N_csv(double x)
{
    std::string filename = "csv/N_DGLAP_x=" + std::to_string(x) + ".csv";
    std::ofstream fout(filename);
    fout << "r2,N_p_DGLAP\n";

    DipoleAmplitude dipole(MZ_IPSAT);
    dipole.EnableLookupTable();

    double rmin = 1e-4;
    double rmax = 10.0;
    int Nr = 300;

    for (int i = 0; i < Nr; ++i) 
    {
        double frac = static_cast<double>(i) / (Nr - 1);
        double r = rmin * pow(rmax / rmin, frac)* CFAC;

        double N_val = N_p(r, x, dipole);

        fout << r*r << "," << N_val << "\n";
    }

    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado com sucesso.\n";
    return filename;
}
}