#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <functional>
#include <map>
#include "ctes.h"
#include "utils.h"
#include <filesystem>

std::string extrair_nome_base(const std::string& caminho)
{
    return std::filesystem::path(caminho).stem().string();
}

Meson::Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
             double NT_, double R2T_, double NL_, double R2L_)
    : meson(m), nome(n), MV(MV_), mf(mf_), ef(ef_),
      NT(NT_), NL(NL_), R2T(R2T_), R2L(R2L_), R2(0.0), isGLC(true) {}

Meson::Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
             double NT_, double NL_, double R2_)
    : meson(m), nome(n), MV(MV_), mf(mf_), ef(ef_),
      NT(NT_), NL(NL_), R2T(0.0), R2L(0.0), R2(R2_), isGLC(false) {}


Meson Jpsi_GLC("Jpsi", "Jpsi_GLC", massa_psi, mc, qJ, 1.23, 6.5, 0.83, 3.0);
Meson phi_GLC ("phi",  "phi_GLC",  massa_phi, ms, qS, 4.75, 16.0, 1.41, 9.7);

Meson Jpsi_BG ("Jpsi", "Jpsi_BG",  massa_psi, mc, qJ, 0.578, 0.575, 2.3);//R2_psi);
Meson phi_BG  ("phi",  "phi_BG",   massa_phi, ms, qS, 0.919, 0.825, 11.2);//R2_phi);

std::string doubleParaString(double valor, int casas) {
    std::ostringstream stream;
    // fixed: garante o número de casas após o ponto
    // setprecision(n): define o número de casas decimais
    stream << std::fixed << std::setprecision(casas) << valor;
    return stream.str();
}

// Derivada por extrapolação de Richardson
double dfridr(
    const std::function<double(double)>& func,
    double x,
    double h,
    double& err
) {
    constexpr int NTAB = 10;
    constexpr double CON  = 1.4;
    constexpr double CON2 = CON * CON;
    constexpr double BIG  = 1.0e30;
    constexpr double SAFE = 2.0;

    if (h == 0.0) {
        throw std::runtime_error("h must be nonzero in dfridr");
    }

    double a[NTAB][NTAB];

    double hh = h;
    a[0][0] = (func(x + hh) - func(x - hh)) / (2.0 * hh);

    err = BIG;
    double dfridr_val = a[0][0];

    for (int i = 1; i < NTAB; ++i) {
        hh /= CON;
        a[0][i] = (func(x + hh) - func(x - hh)) / (2.0 * hh);

        double fac = CON2;

        for (int j = 1; j <= i; ++j) {
            a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
            fac *= CON2;

            double errt = std::max(
                std::abs(a[j][i] - a[j - 1][i]),
                std::abs(a[j][i] - a[j - 1][i - 1])
            );

            if (errt <= err) {
                err = errt;
                dfridr_val = a[j][i];
            }
        }

        if (std::abs(a[i][i] - a[i - 1][i - 1]) >= SAFE * err) {
            return dfridr_val;
        }
    }

    return dfridr_val;
}
double derivative_richardson(const std::function<double(double)>& f,
                             double x,
                             double h)
{
    double D1 = (f(x + h) - f(x - h)) / (2.0*h);
    double D2 = (f(x + h/2) - f(x - h/2)) / h;

    return (4.0*D2 - D1) / 3.0;
}

double derivative_poly5(
    const std::function<double(double)>& f,
    double x)
{
    double h = 1e-10;

    double f1=f(x-2*h);
    double f2=f(x-h);
    double f3=f(x+h);
    double f4=f(x+2*h);

    return (-f4 + 8*f3 - 8*f2 + f1)/(12*h);
}

// ----------------- conversão W <-> x ----------------

double x_to_W(double x, const Meson& M) {
    return M.MV / std::sqrt(x);
}

double W_to_x(double W, const Meson& M) {
    return (M.MV * M.MV) / (W * W);
}



std::map<std::string, MesonModels> meson_models = {
    {"Jpsi", {Jpsi_GLC, Jpsi_BG}},
    {"phi", {phi_GLC, phi_BG}},
};

// ----------------- função escolhe meson ----------
Meson input_meson()
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

    return it->second.M_GLC; // Retorna o modelo GLC por padrão
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

// ---------------- slope B(Q2) ----------------
double B(double x, double Q2, const Meson& M) {
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