
#include <Minuit2/Minuit2Minimizer.h>
#include "Math/Functor.h"
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
#include <filesystem>

#include "utils.h"
#include "correcs.h"
#include "ctes.h"
#include "utils.h"
#include "GBW.h"
#include "plot.h"
#include "dipoleamplitude.hpp"

// ---------------- minimização ----------------

struct Chi2 {
    const std::vector<double>& W_exp;
    const std::vector<double>& sigma_exp;
    const std::vector<double>& error;

    double operator()(const double* par) const {
        double omega = par[0];
        double B     = par[1];

        double chi2 = 0.0;
std::cout << "N pontos = " << W_exp.size() << std::endl;

    std::vector<double> x(W_exp.size());
    for(size_t i=0; i<W_exp.size(); ++i){
        x[i] = W_to_x(W_exp[i], phi_BG_GBW);
    }

    for(size_t i=0; i<W_exp.size(); ++i){
        if (!std::isfinite(x[i]) || x[i] <= 0 || x[i] >= 1) {
            std::cerr << "x inválido: " << x[i] << " para W=" << W_exp[i] << std::endl;
        }
    }
        for(size_t i=0; i<W_exp.size(); ++i){
    double model = GBW::sigma_model(x[i], B, omega, phi_BG_GBW);
    double diff = model - sigma_exp[i];

    std::cout 
        << "i=" << i
        << " model=" << model
        << " exp=" << sigma_exp[i]
        << " err=" << error[i]
        << " contrib=" << (diff*diff)/(error[i]*error[i])
        << std::endl;

    chi2 += (diff * diff) / (error[i] * error[i]);
}
        return chi2;
    }
};

void minimizar_chi2(){

    std::vector<int> dataset;
    std::vector<double> W_exp, sigma_exp, error;
    read_sigma_exp("csv/expdata/sigma/phi_fit_data.csv",dataset, W_exp, sigma_exp, error);

    for(size_t i=0; i<W_exp.size(); ++i){
    std::cout 
        << "W=" << W_exp[i]
        << " sigma_exp=" << sigma_exp[i]
        << " err=" << error[i]
        << std::endl;
}



    Chi2 chi2_func{W_exp, sigma_exp, error};

    ROOT::Math::Functor fcn(chi2_func, 2);
    ROOT::Minuit2::Minuit2Minimizer minimizer;

    minimizer.SetFunction(fcn); 

    minimizer.SetVariable(0, "omega", 0.5, 0.01); // chute inicial e passo
    minimizer.SetVariable(1, "B", 1.0, 0.1);    

    minimizer.SetVariableLimits(0, -5.0, 5.0);
    minimizer.SetVariableLimits(1, -0.99, 10.0);

    bool success;

    success = minimizer.Minimize();


if (!success) {
    std::cerr << "Minimização falhou!" << std::endl;
    return;
}

const double* xs = minimizer.X();
const double* errs = minimizer.Errors();

    std::cout << "omega = " << xs[0] << " ± " << errs[0] << std::endl;
    std::cout << "B     = " << xs[1] << " ± " << errs[1] << std::endl;

    std::cout << "Status: " << minimizer.Status() << std::endl;
    double chi2 = minimizer.MinValue();
    double ndof = W_exp.size() - 2;

    std::cout << "chi2/ndof = " << chi2 / ndof << std::endl;
}


double calc_chi2(double B, double omega){
    std::vector<int> dataset;
    std::vector<double> W_exp, sigma_exp, error;

    read_sigma_exp("csv/expdata/sigma/phi_fit_data.csv",
                   dataset, W_exp, sigma_exp, error);

    Chi2 chi2_func{W_exp, sigma_exp, error};

    double par[2] = {omega, B};
    double chi2 = chi2_func(par);

    double ndof = W_exp.size() - 2;

    std::cout << "Chi2: " << chi2
              << ", ndof: " << ndof
              << ", chi2/ndof: " << chi2/ndof
              << std::endl;

    return chi2;
}