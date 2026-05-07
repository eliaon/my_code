#ifndef LHAPDF_WRAPPER_HPP
#define LHAPDF_WRAPPER_HPP

#include <memory>
#include <string>

// Forward declaration (evita incluir LHAPDF no header)
namespace LHAPDF {
    class PDF;
}

class DipoleModel {
public:
    // Construtor
    DipoleModel(const std::string& setname, int member = 0);
    ~DipoleModel();

    // --- PDF access ---
    double gluon(double x, double Q2) const;     // g(x,Q²)
    double alpha_s(double Q2) const;

    // --- Escala ---
    double mu2(double r) const;

    // --- Dipole ---
    double sigma_dip(double x, double r) const;
    double N(double x, double r) const;  // amplitude normalizada
    void N_csv(); // gera CSV para N(r) em dado x
    std::string xf_vs_x(double Q2, int flavor); // gera CSV para xf(x) em dado Q²

    // --- Parâmetros físicos ---
    void set_sigma0(double val);
    void set_C(double val);
    void set_mu0_2(double val);

private:
    std::unique_ptr<LHAPDF::PDF> pdf;

    // parâmetros do modelo
    double sigma0;   // mb
    double C;
    double mu0_2;
};

#endif