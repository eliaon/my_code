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
#include "DGLAP.hpp"
#include "bCGC.h"
#include "LHAPDF.hpp"
#include "IIM.hpp"

using namespace MZ_ipsat;

int main(){
    int threads = 20;
    omp_set_num_threads(threads);// Ajuste conforme o número de núcleos disponíveis
    //draw_wavefunctions(Jpsi_BG);
    //print_B_values(Jpsi_GLC, 0.0);
    //N_plot();
    //plot_sigmaqq();
    //run_rapidez_plot();
    //run_sigma_csv_GBW();
    //plot_sigma_phi("csv/phi_sigma_GBW(10)_Q2=0.csv");
    //plot_sigma(std::string("phi"));    //debug_correc();
    //plot_overlap();
    //printf("Qs = %g\n", QS_bCGC(1e-4,0));
    //dsigma_dump();
    //debug_correc();
    //calc_chi2();
    //minimizar_chi2();
    //run_sigma_A_csv_GBW();
    //print_TA_table();
    //check_TA();
    //plot_sigma_Jpsi_gammaPb_GBW("csv/Jpsi_sigma-gammaPb_GBW(10)_Q2=0.csv");
    //plot_TA_b();
    //rapidez_PbPb_csv();
    //plot_rapidez_PbPb_phi("csv/phi_rapidez_PbPb_5.36e+03GeV.csv", 5360.0);
    //sigma_p_DGLAP_csv(); // GLC e BG idênticos para DGLAP, pois o modelo é só para o proton
    //dipole_DGLAP::N_p_csv(); // gera os arquivos csv para N(r^2) do modelo DGLAP e plota a curva
    //overlap_csv_fc();
    //sigma_ipsat_integrado_csv(); // gera o csv da seção de choque integrada para ipsat e plota a curva
    //IPSAT::dsigma_dump(); // gera os csvs de dsigma/dt para os W escolhidos e plota as curvas
    //compare_N_models(1e-4); // gera os csvs de N(r) para os modelos escolhidos e plota as curvas
    //plot_sigma_Jpsi(bCGC::sigma_p_csv(), "bCGC"); // gera o csv da seção de choque integrada para bCGC e plota a curva para J/psi
    //plot_xf_multiQ2();

    //bCGC::sigma_p_csv();
    IIM::sigma_p_csv();
    //DipoleModel model("CT14lo");
    //model.N_csv();

    //plot_sigma_Jpsi("csv/Jpsi_sigma_gammap_IIM.csv", "IIM");

    return 0;
}