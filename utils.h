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
extern Meson Jpsi_GLC_GBW;
extern Meson phi_GLC_GBW;
extern Meson Jpsi_BG_GBW;
extern Meson phi_BG_GBW;

extern Meson Jpsi_GLC_ipsat;
extern Meson phi_GLC_ipsat;
extern Meson Jpsi_BG_ipsat;
extern Meson phi_BG_ipsat;

struct MesonModelsGBW{
    const Meson& M_GLC;
    const Meson& M_BG;
};

struct MesonModelsipsat{
    const Meson& M_GLC;
    const Meson& M_BG;
};

extern std::map<std::string, MesonModelsGBW> meson_modelsGBW;
extern std::map<std::string, MesonModelsipsat> meson_modelsipsat;

// funções
std::string doubleParaString(double x);

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

Meson input_meson(const std::string model);

void perfil(const Meson& M);

double B_slope(double x, double Q2, const Meson& M);

double calc_chi2(double B, double omega);

void minimizar_chi2();

void read_sigma_exp(const std::string& filename,
                    std::vector<int>& dataset,
                    std::vector<double>& W,
                    std::vector<double>& sigma,
                    std::vector<double>& error);

std::string timestamp();

void read_two_columns(const std::string& filename,
                      std::vector<double>& col1,
                      std::vector<double>& col2);

void read_csv(const std::string& filename,
              std::vector<double>& col1,
              std::vector<double>& col2,
              std::vector<double>& col3);


void read_xyscan_csv(const std::string& filename,
                    std::vector<double>& x,
                    std::vector<double>& y,
                    std::vector<double>& error);

bool load_set(
    const std::string& xstr,
    std::vector<double>& r2_ipsat,
    std::vector<double>& N_ipsat,
    std::vector<double>& r2_gbw,
    std::vector<double>& N_gbw,
    std::vector<double>& r2_iim,
    std::vector<double>& N_iim);

std::string get_meson();

std::string vec_to_pylist(const std::vector<double>& vec);
std::string flavorToString(int flavor);
std::string flavorName(int flavor);

void read_rapidity_hepdata(
    const std::string& filename,
    std::vector<double>& y,
    std::vector<double>& dsdy,
    std::vector<double>& err);

#endif