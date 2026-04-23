#include "utils.h"
#include "ctes.h"
#include "GBW.h"
#include "fstream"
#include "wavefunctions.h"
#include "integration.hpp"
#include "correcs.h"
#include "plot.h"

#include <sstream>
#include <cmath>
#include <iostream>
#include <format>
#include <chrono>
#include <boost/math/special_functions/bessel.hpp>


parametros_GBW gbw(23*mb_to_gev2, 3e-4, 0.29);


// ----------- parametros GBW NEW 2018 Saturation model of DIS: an update https://doi.org/10.1007/JHEP03(2018)102
// sequência: sigma_0 , x_0, lambda
parametros_GBW gbw_5(28.18*mb_to_gev2, 0.31e-4, 0.237);

parametros_GBW gbw_10(27.43*mb_to_gev2, 0.40e-4, 0.248);

parametros_GBW gbw_20(26.60*mb_to_gev2, 0.53e-4, 0.259);

parametros_GBW gbw_50(25.21*mb_to_gev2, 0.80e-4, 0.281);

TA_Table precompute_TA(int Nb, double bmax, double R, double a, double rho0)
{
    TA_Table table;
    table.b_vals.resize(Nb);
    table.TA_vals.resize(Nb);

    double zmax = 20.0 * 5.07; // fm para GeV^-1, suficiente para Pb-208
    int Nz = 200;

    for (int i = 0; i < Nb; ++i)
    {
        double b = bmax * i / (Nb - 1);
        table.b_vals[i] = b;

        auto integrand = [b, R, a, rho0](double z) {
            double r = sqrt(b*b + z*z);
            return GBW::rho_WS(r, R, a, rho0);
        };

        table.TA_vals[i] = integrate_simpson(integrand, -zmax, zmax, Nz);
    }

    return table;
}



double interpolate_TA(double b, const TA_Table& table)
{
    int N = table.b_vals.size();

    if (b <= table.b_vals[0]) return table.TA_vals[0];
    if (b >= table.b_vals[N-1]) return table.TA_vals[N-1];
    // busca binária
    int low = 0, high = N - 1;
    while (high - low > 1)
    {
        int mid = (low + high) / 2;
        if (table.b_vals[mid] > b)
            high = mid;
        else
            low = mid;
    }

    double b1 = table.b_vals[low];
    double b2 = table.b_vals[high];
    double T1 = table.TA_vals[low];
    double T2 = table.TA_vals[high];

    return T1 + (T2 - T1) * (b - b1) / (b2 - b1);
}



namespace GBW {
double QS2( double x, parametros_GBW params)
{
    return pow(params.x0 / x, params.lambda);
}   

double Np(double r,  double x, parametros_GBW params)
{
    double Qs2 = QS2(x, params);
    double arg = (r * r) * Qs2 / 4.0;
    return (1.0 - exp(-arg))*std::pow(1.0-x, 5.26);
}

double sigma_qq(double r, double x, parametros_GBW params)
{
    return params.sigma0 * Np(r, x, params);
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
        double sigma = sigma_qq(r, x, gbw);
        fout << r * CFAC << "," << sigma * mb_to_gev2 << "\n"; // converte r para fm e sigma para mb
    }
}

double amplitude(double x, double Q2, const Meson& M,
                const parametros_GBW& gbw,
                 int Nr , int Nz,
                 double rmin , double rmax)
{
    auto amp_r = [x, Q2, M, gbw, Nz](double r) {
        double Ov = overlap_r(r, Q2, M, Nz);
        double sigma_dip = sigma_qq(r, x, gbw);
        double sqrt_fc = std::sqrt(f_c(r, -0.979599, 0.403569));
        return 0.5 * r * Ov * sigma_dip * sqrt_fc; // r de d²r = 2π r dr
    };
    double amp = integrate_simpson(amp_r, rmin, rmax, Nr);
    return amp;
}

double sigma_x(double x, double Q2 , const Meson& M,
               int Nr, int Nz, 
               double rmin, double rmax)
{
    parametros_GBW params = gbw_10; // ou escolha outro conjunto de parâmetros se desejar
    double amp = amplitude(x, Q2, M, params, Nr, Nz, rmin, rmax);
    double B_val = B_slope(x, Q2, M);
    double lambda_e = calculate_lambda(x, Q2, M, "GBW");
    double RG_val = RG(x, Q2, lambda_e, M);
    double beta_val = beta(x, Q2, lambda_e, M);
    std::cout<< "x =" << x << " | Amplitude GBW: " << amp << std::endl;
    std::cout << "lambda_e: " << lambda_e << " | Rg: " << RG_val << " | beta_corr: " << beta_val << std::endl;
    std::cout << "B_slope: " << B_val << std::endl;

    double correction_factor = RG_val * RG_val * (1.0 + beta_val * beta_val);
    return correction_factor * (amp * amp) / (16.0 * M_PI * B_val);
}

// --------------- modelos para parametrizar B e omega ------------

double amplitude_model(double x,  double B, double omega,
                        const Meson& M, double Q2,  int Nr, 
                        int Nz, double rmin, double rmax)
{
    auto amp_r = [x, Q2, M, B, omega, Nz](double r) {
        double Ov = overlap_r(r, Q2, M, Nz);
        double sigma_qq = GBW::sigma_qq(r, x, gbw);
        double sqrt_fc = std::sqrt(f_c(r, B, omega));
        return 0.5 * r * Ov * sigma_qq * sqrt_fc; // r de d²r = 2π r dr
    };
    double amp = integrate_simpson(amp_r, rmin, rmax, Nr);
    return amp;
}

double sigma_model(double x, double B, double omega, const Meson& M, double Q2,
                     int Nr, int Nz, double rmin, double rmax)
    {
     parametros_GBW params = gbw_10; // ou escolha outro conjunto de parâmetros se desejar
     double amp = GBW::amplitude_model(x, B, omega, M, Q2, Nr, Nz, rmin, rmax);
     double B_slope_val = B_slope(x, Q2, M);
     double lambda_e = calculate_lambda(x, Q2, M, "GBW");
     double RG_val = RG(x, Q2, lambda_e, M);
     double beta_val = beta(x, Q2, lambda_e, M);
    
     double correction_factor = RG_val * RG_val * (1.0 + beta_val * beta_val);
     double sigma_gev = correction_factor * (amp * amp) / (16.0 * M_PI * B_slope_val);
     double sigma_nb = sigma_gev * GeV2_to_nb;
     return sigma_nb;
    }





void sigma_p_csv()
{
    double Q2;
    std::cout << "Insira o valor de Q2 (GeV^2): ";
    std::cin >> Q2;

    const Meson& M_GLC = input_meson("GBW"); //importa o glc do meson escolhido
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string filename = "csv/" + M_GLC.meson + "_sigma_GBW(10)_Q2=" + Q2_str + ".csv";
    std::ofstream fout(filename);
    fout << "W,sigma_GLC,sigma_BG\n";

    const int Nw = 10;
    double Wmin = x_to_W(1e-2, M_GLC); // W mínimo correspondente a x=1e-5 para o meson escolhido
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
        sigma_GLC_v[i] = sigma_x(x, Q2, M_GLC) * GeV2_to_nb;
        sigma_BG_v[i]  = sigma_x(x, Q2, M_BG) * GeV2_to_nb;
        std::cout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";

    }

    for (int i = 0; i < Nw; ++i) {
    fout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
}

    fout.close();
    auto end = clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução: " << duration << " ms" << std::endl;
    std::cout << "Arquivo '" << filename << "' gerado em " << duration << " ms" << std::endl;

    plot_sigma(M_GLC.meson, filename, "GBW");
}


//------------------ caso nuclear -----------------------

double rho_WS(double r, double R, double a, double rho0)
{
    return rho0 / (1.0 + exp((r - R)/a));
}






double N_A(double r, double x, double b, parametros_GBW params, const TA_Table& table)
{
    double sigmaqq_p = params.sigma0 * Np(r, x, params);
    double TA = interpolate_TA(b, table); 
    double arg = 0.5 * TA * sigmaqq_p;
    return (1.0 - exp(-arg));   
}

double sigma_qq_A(double r, double x, double b, parametros_GBW params, const TA_Table& table)
{
    auto b_integrand = [r, x, params, &table](double b) {
        return 2.0 * M_PI * b * N_A(r, x, b, params, table);
    };
    double bmax = table.b_vals.back(); // limite de integração em b
    int nb = table.b_vals.size(); // número de pontos para integração em b
    return integrate_simpson(b_integrand, 0.0, bmax, nb);
}

double sigma_A(double x, double Q2, const Meson& M, const TA_Table& table,
                    int Nr, int Nz, double rmin, double rmax)
{
    parametros_GBW params = gbw_10; // ou escolha outro conjunto de parâmetros se desejar

    auto amp_b2 = [x, Q2, M, params, Nz, &table, rmin, rmax, Nr](double b){
        auto amp_r = [x, Q2, b, M, params, Nz, &table](double r) {
        double sqrt_fc = std::sqrt(f_c(r, -0.979599, 0.403569)); // valores de B e omega fitados no minuit para o proton
        double Ov = overlap_r(r, Q2, M, Nz);
        double sigma_qq = N_A(r, x, b, params, table); 
        return 0.5 * r * Ov * sigma_qq * sqrt_fc; // r de d²r = 2π r dr
    };
    double amp = integrate_simpson(amp_r, rmin, rmax, Nr);
    //double B_val = B_slope(x, Q2, M);
    //double lambda_e = calculate_lambda(x, Q2, M);
    //double RG_val = RG(x, Q2, lambda_e, M);
    //double beta_val = beta(x, Q2, lambda_e, M);

    double result = (amp * amp); // (16.0 * M_PI)); //* RG_val * RG_val * (1.0 + beta_val * beta_val);

    return 2 * M_PI * b * result; // quadrado da amplitude

    };
    double bmax = table.b_vals.back(); // limite de integração em b
    int nb = table.b_vals.size(); // número de pontos para integração em b
    return integrate_simpson(amp_b2, 0.0, bmax, nb);
    
};


double integral_rho(double R, double a)
{
    double rmax = 20.0; // fm (suficiente para Pb)
    int Nr = 500;

    auto integrand = [R, a](double r) {
        double rho = 1.0 / (1.0 + exp((r - R)/a)); // rho0 = 1
        return r*r * rho;
    };

    double integral = integrate_simpson(integrand, 0.0, rmax, Nr);

    return 4.0 * M_PI * integral;
}

double compute_rho0(int A, double R, double a)
{
    double I = integral_rho(R, a);
    return A / I;
}


void sigma_A_csv()
{
    double Q2=0.0;

    const Meson& M_GLC = input_meson("GBW"); //importa o glc do meson escolhido
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG; //pega o bg correspondente ao meson escolhido

    std::string Q2_str = std::format("{:.3g}", Q2);
    std::string filename = "csv/" + M_GLC.meson + "_sigma-gammaPb_GBW(10)_Q2=" + Q2_str + ".csv";
    std::ofstream fout(filename);
    fout << "W,sigma_GLC,sigma_BG\n";

    const int Nw = 100;
    double Wmin = x_to_W(1e-2, M_GLC);
    double Wmax = 3e3;

    using clock = std::chrono::steady_clock;
    auto start = clock::now();

    double R = 6.62*5.07; // fm para GeV^-1
    double a = 0.546*5.07; // fm para GeV^-1
    double rho0 = 0.1603/std::pow(5.07, 3); // fm^-3 para GeV^3 

    TA_Table table = precompute_TA(200, 50.0, R, a, rho0); // pré-calcula T_A(b) para Pb-208

    std::vector<double> Wv(Nw), sigma_GLC_v(Nw), sigma_BG_v(Nw);


    std::cout << "W (Gev), sigma_GLC (nb), sigma_BG (nb)\n";
    #pragma omp parallel for schedule(dynamic)
for (int i = 0; i < Nw; ++i) {

    double frac = static_cast<double>(i) / (Nw - 1);
    double W = Wmin * pow(Wmax / Wmin, frac);
    double x = (M_GLC.MV * M_GLC.MV) / (W * W);

    Wv[i] = W;
    sigma_GLC_v[i] = GBW::sigma_A(x, Q2, M_GLC, table) * GeV2_to_nb;
    sigma_BG_v[i]  = GBW::sigma_A(x, Q2, M_BG, table) * GeV2_to_nb;
    std::cout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
}
//for (int i = 0; i < Nw; ++i) {
//    std::cout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
//}
    for (int i = 0; i < Nw; ++i) {
    fout << Wv[i] << "," << sigma_GLC_v[i] << "," << sigma_BG_v[i] << "\n";
}

    fout.close();
    auto end = clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Tempo de execução: " << duration << " ms" << std::endl;
    std::cout << "Arquivo '" << filename << "' gerado em " << duration << " ms" << std::endl;
if (M_GLC.meson == "Jpsi")
    plot_sigma_Jpsi_gammaPb_GBW(filename);
else if (M_GLC.meson == "phi")
    plot_sigma_phi_gammaPb_GBW(filename);
}

double dN_domega(double omega, double sqrt_s)
{
    double Z = 82.0; // número atômico do Pb
    double Z2 = Z * Z;
    double alpha = 1.0 / 137.0; // constante de estrutura fina
    double gamma = sqrt_s / (2.0 * 0.938); // fator de Lorentz para Pb (massa ~ 0.938 GeV)
    double bmin = 2.0 * 6.62 * 5.07; // raio do Pb em fm, convertido para GeV
    double ksi = (omega * bmin) / gamma;
    double K0 = boost::math::cyl_bessel_k(0, ksi);
    double K1 = boost::math::cyl_bessel_k(1, ksi);
    double coeff = (2 * Z2 * alpha) / (M_PI * omega);
    double flux = coeff * (ksi * K0 * K1 - (ksi * ksi / 2.0) * (K1 * K1 - K0 * K0));
    //std::cout << "gamma = " << gamma << std::endl;
    //std::cout << "ksi = " << ksi << std::endl;
    return flux;
    
}

double d_sigma_dy_PbPb(double y, double sqrt_s, double Q2, const Meson& M, TA_Table& table)
{

    double omega_plus = (M.MV / 2.0) * exp(y);
    double omega_minus = (M.MV / 2.0) * exp(-y);

    double W_plus = std::sqrt(2 * omega_plus * sqrt_s);
    double W_minus = std::sqrt(2 * omega_minus * sqrt_s);

    double x_plus = (M.MV * M.MV) / (W_plus * W_plus);
    double x_minus = (M.MV * M.MV) / (W_minus * W_minus);


    double n_plus  = omega_plus  * dN_domega(omega_plus, sqrt_s);
    double n_minus = omega_minus * dN_domega(omega_minus, sqrt_s);

    //std::cout<<"omega = " << omega_plus << "n(omega)=" << n_plus << std::endl;

    double sigma_plus  = sigma_A(x_plus, Q2, M, table);
    double sigma_minus = sigma_A(x_minus, Q2, M, table);
    

    return n_plus * sigma_plus + n_minus * sigma_minus;
}

void rapidez_PbPb_csv()
{
    double sqrt_s = 2.76e3; // GeV

    const Meson& M_GLC = input_meson("GBW");
    const Meson& M_BG  = meson_modelsGBW.find(M_GLC.meson)->second.M_BG;

    std::string sqrt_s_str = std::format("{:.3g}", sqrt_s);
    std::string filename = "csv/" + M_GLC.meson +
        "_rapidez_PbPb_" + sqrt_s_str + "GeV_fc.csv";

    double Q2 = 0.0;

    double R = 6.62*5.07;
    double a = 0.546*5.07;
    double rho0 = 0.1603/std::pow(5.07, 3);

    TA_Table table = precompute_TA(200, 50.0, R, a, rho0);

    const int Ny = 100;
    double ymax = std::log(sqrt_s / M_GLC.MV);

    int Ntot = 2*Ny + 1;

    std::vector<double> Y(Ntot), GLC(Ntot), BG(Ntot);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i <= Ny; ++i)
    {
        double y = (static_cast<double>(i) / Ny) * ymax;

        double dsdy_GLC = d_sigma_dy_PbPb(y, sqrt_s, Q2, M_GLC, table) * gev2_to_mb;
        double dsdy_BG  = d_sigma_dy_PbPb(y, sqrt_s, Q2, M_BG, table) * gev2_to_mb;

        int idx_pos = Ny + i;
        int idx_neg = Ny - i;

        std::cout << "y: " << y << ", dσ/dy GLC: " << dsdy_GLC << " mb, dσ/dy BG: " << dsdy_BG << " mb\n";  

        Y[idx_pos] =  y;
        GLC[idx_pos] = dsdy_GLC;
        BG[idx_pos]  = dsdy_BG;

        // simetria
        Y[idx_neg] = -y;
        GLC[idx_neg] = dsdy_GLC;
        BG[idx_neg]  = dsdy_BG;
    }

    // ✔ escrita segura
    std::ofstream fout(filename);
    fout << "y,d_sigma_dy_GLC,d_sigma_dy_BG\n";

    for (int i = 0; i < Ntot; ++i)
    {
        fout << Y[i] << "," << GLC[i] << "," << BG[i] << "\n";
    }

    fout.close();

    std::cout << "Arquivo '" << filename << "' gerado.\n";

    if (M_GLC.meson == "Jpsi")
        plot_rapidez_PbPb_Jpsi(filename, sqrt_s);
    else if (M_GLC.meson == "phi")
        plot_rapidez_PbPb_phi(filename, sqrt_s);

    std::cout << "Plot gerado para " << M_GLC.meson
              << " em Pb-Pb a sqrt(s) = "
              << sqrt_s_str << " GeV.\n";
}

string N_csv(double x)
{
    std::string filename = "csv/N_GBW_x=" + doubleParaString(x) + ".csv";
    std::ofstream fout(filename);
    fout << "r,N\n";

    parametros_GBW params = gbw_10; // ou escolha outro conjunto de parâmetros se desejar

    const int Npoints = 5000;
    double rmin = 1e-4, rmax = 10.0;

    for (int i = 0; i < Npoints; ++i)
    {
        double frac = (double)i / (Npoints - 1);
        double r = rmin * std::pow(rmax / rmin, frac)*CFAC; // escala log
        double N_val = Np(r, x, params); // exemplo para x=1e-4 e x0=1e-2
        fout << r * r << "," << N_val << "\n"; // converte r para fm
    }
    return filename;
}







}