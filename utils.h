#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <functional>
#include <map>
#include "ctes.h"
#include <string>

std::string extrair_nome_base(const std::string& caminho);

class Meson {
public:
    std::string meson;
    std::string nome;
    double MV;
    double mf;
    double ef;
    double NT;
    double NL;
    double R2T;
    double R2L;
    double R2;
    bool isGLC;

    Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
          double NT_, double R2T_, double NL_, double R2L_);

    Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
          double NT_, double NL_, double R2_);
};

// 🔹 Só declaração (extern)
extern Meson Jpsi_GLC;
extern Meson phi_GLC;
extern Meson Jpsi_BG;
extern Meson phi_BG;

struct MesonModels{
    const Meson& M_GLC;
    const Meson& M_BG;
};

extern std::map<std::string, MesonModels> meson_models;

// funções
std::string doubleParaString(double valor, int casas = 2);

double dfridr(const std::function<double(double)>& func,
              double x,
              double h,
              double& err);

double derivative_richardson(const std::function<double(double)>& f,
                             double x,
                             double h = 1e-4);

double derivative_poly5(const std::function<double(double)>& f,
                        double x);

double x_to_W(double x, const Meson& M);
double W_to_x(double W, const Meson& M);

Meson input_meson();

void perfil(const Meson& M);

double B(double x, double Q2, const Meson& M);
#endif