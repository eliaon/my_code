#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <ctime>
#include <Python.h>
#include <filesystem>

#include "matplotlib-cpp/matplotlibcpp.h"
#include "utils.h"
#include "bCGC.h"
#include "GBW.h"
#include "DGLAP.hpp"
#include "ipsat.h"

namespace plt = matplotlibcpp;



// ---------- FUNÇÃO PARA PUXAR A STRING TIMESTAMP (HORÁRIO E DIA/MES/ANO AGORA)









// ------------- PLOTA CURVA PARA dσ/dt PARA UMA DADA ENERGIA EM GeV

void dsigma_dt_curve(const std::string& W)
{
    std::string filename = "csv/Jpsi_dsigma_dt_W=" + W + "GeV.csv";

    std::vector<double> t, glc, bg;

    read_csv(filename, t, glc, bg);

    if(W == "100")
    {
        for(auto& v : glc) v *= 5.0;
        for(auto& v : bg)  v *= 5.0;

        plt::plot(t, glc, {{"label","GLC"}});
        plt::plot(t, bg, {{"label","BG"}});
    }
    else
    {
        plt::plot(t, glc,
            {{"label","GLC (W=" + W + " GeV)"}});
        plt::plot(t, bg,
            {{"label","BG (W=" + W + " GeV)"},{"linestyle","--"}});
    }
}


// ------------- PLOTA DADOS EXPERIMENTAIS PARA dσ/dt DE FOTOPRODUÇÃO DE φ

void plot_dsigma_exp_data_phi()
{
    std::ifstream file("csv/expdata/HEPData-dsigdt_gammap_phi.csv");

    std::vector<double> t;
    std::vector<double> dsdt;
    std::vector<double> err;

    std::string line;

    while(std::getline(file, line))
    {
        // ignorar comentários HEPData
        if(line.empty() || line[0] == '#')
            continue;

        // ignorar header textual
        if(!std::isdigit(line[0]))
            continue;

        std::stringstream ss(line);
        std::string value;

        std::vector<double> row;

        while(std::getline(ss, value, ','))
            row.push_back(std::stod(value));

        if(row.size() < 5)
            continue;

        double t_val = row[0];
        double sigma = row[3] * 1000.0;   // µb → nb
        double error = std::abs(row[4]) * 1000.0;

        t.push_back(t_val);
        dsdt.push_back(sigma);
        err.push_back(error);
    }

    plt::errorbar(t, dsdt, err,
    {{"fmt","o"},
     {"color","black"},
     {"label","HERA (W≈70 GeV)"}});
}

// ------------- PLOTAR AS CURVAS dσ/dt PARA TODOS OS W ESCOLHIDOS

void plot_dsigma_dt(std::string meson)
{
    std::vector<std::string> w_values = {"70"};

    plt::figure_size(700,500);

    for(const auto& W : w_values)
        dsigma_dt_curve(W);

    // plota dados experimentais apenas para o méson phi
    if(meson == "phi")
        plot_dsigma_exp_data_phi();

    plt::xlim(0.0, 2.5);
    plt::ylim(1e-3, 1e4);

    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.xscale('linear')\n"
    "plt.yscale('log')\n"
    );

    plt::xlabel("$|t| (GeV^2)$");
    plt::ylabel("$d\\sigma/dt$ (nb)");

    
    plt::grid(true);
    plt::legend();

    plt::save("plots/dsigma_dt/"+ meson +"dsigma_dt_" + timestamp() + ".png");
    plt::show();
}



// -------------- PLOTA AS CURVAS DE N

void plot_N_models(const std::vector<std::pair<std::string,std::string>>& files, double x)
{
    plt::figure_size(700,500);

    for (const auto& [fname, label] : files)
    {
        std::vector<double> r2, N;

        std::ifstream fin(fname);
        if (!fin.is_open()) {
            std::cerr << "Erro ao abrir: " << fname << std::endl;
            continue;
        }

        std::string line;
        std::getline(fin, line); // pula header

        while (std::getline(fin, line))
        {
            std::stringstream ss(line);
            std::string val1, val2;

            std::getline(ss, val1, ',');
            std::getline(ss, val2, ',');

            r2.push_back(std::stod(val1));
            N.push_back(std::stod(val2));
        }

        plt::plot(r2, N, {{"label", label}});
    }

    plt::xlabel("$r^2$");
    plt::ylabel("$N_p$");

    plt::title("Amplitude de dipolo $N_p(r^2)$ para diferentes modelos (x=" + doubleParaString(x) + ")");

    plt::xlim(0.0001,100.0);
    plt::ylim(0.0001,1.2);

    PyRun_SimpleString( "import matplotlib.pyplot as plt\n"
                        "plt.yscale('log')\n"
                        "plt.xscale('log')\n"

    );

    plt::grid(true);
    plt::legend();

    std::string filename =
        "plots/N/N_models_" + timestamp() + "_x=" + doubleParaString(x) + ".png";

    plt::save(filename);
    plt::show();
}

void compare_N_models(double x)
{
    std::vector<std::pair<std::string,std::string>> files = {
        {IPSAT::N_csv(x), "IPsat"},
        {bCGC::N_csv(x), "bCGC"},
        {GBW::N_csv(x), "GBW"},
        {DGLAP::N_csv(x), "DGLAP"}
    };

    plot_N_models(files, x);
}

// -------------- PLOTA CURVAS DE RAPIDEZ PARA UM SQRT S ESCOLHIDO



void plot_rapidity()
{
    std::string meson = get_meson();
    std::string filename =
        "csv/" + meson + "_rapidez_5.26e+03GeV.csv";

    std::vector<double> Y, rap_GLC, rap_BG;

    read_csv(filename, Y, rap_GLC, rap_BG);

    plt::figure_size(800,600);

    plt::plot(Y, rap_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(Y, rap_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}});

    plt::xlim(-8,8);

    plt::xlabel("Y");
    plt::ylabel("dσ/dY [nb]");

    double s = 5.26e3;

    std::string meson_label = meson;
    if(meson == "Jpsi")
        meson_label = "psi";

    std::stringstream title;
    title << "Distribuição de rapidez do $" << meson_label
          << "$ em pp a "
          << std::fixed << std::setprecision(2)
          << s/1e3 << " TeV";

    plt::title(title.str());

    plt::grid(true);
    plt::legend();

    std::string out =
        "plots/Rapidez/" + meson +
        "_rapidez_5p26TeV_" +
        timestamp() + ".png";

    plt::save(out);
    plt::show();
}

// ------------------ PĹOT OVERLAPS

void plot_Jpsi_overlap()
{
    std::string filename = "csv/Jpsi_overlap_r.csv";

    std::vector<double> r, overlap_GLC, overlap_BG;

    read_csv(filename, r, overlap_GLC, overlap_BG);

    plt::figure_size(1000,600);

    plt::plot(r, overlap_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.8"}});

    plt::plot(r, overlap_BG,
        {{"label","BG"},
         {"color","black"},
         {"linestyle","-"},
         {"linewidth","1.2"}});

    plt::xlim(0.001,1.0);
    plt::ylim(0.0,0.025);

    plt::xlabel("$r$ [fm]");
    plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");

    // escala log em r
    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.xscale('log')\n"
    );

    plt::title("Sobreposição de funções de onda ($\\gamma p \\to J/\\psi p$)");

    plt::grid(true);
    plt::legend();

    std::string out =
        "plots/overlap/Jpsi_overlap_" +
        timestamp() + ".png";

    plt::save(out);
    plt::show();
}

void plot_phi_overlap()
{
    std::string filename = "csv/phi_overlap_r.csv";

    std::vector<double> r, overlap_GLC, overlap_BG;

    read_csv(filename, r, overlap_GLC, overlap_BG);

    plt::figure_size(1000,600);

    plt::plot(r, overlap_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.8"}});

    plt::plot(r, overlap_BG,
        {{"label","BG"},
         {"color","black"},
         {"linestyle","-"},
         {"linewidth","1.2"}});

    plt::xlim(0.01,3.0);
    plt::ylim(0.0,0.01);

    plt::xlabel("$r$ [fm]");
    plt::ylabel("Overlap $r \\phi_V \\phi_{\\gamma}$");

    // escala log em r
    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.xscale('log')\n"
    );

    plt::title("Sobreposição de funções de onda ($\\gamma p \\to \\phi p$)");

    plt::grid(true);
    plt::legend();

    std::string out =
        "plots/overlap/phi_overlap_" +
        timestamp() + ".png";

    plt::save(out);
    plt::show();
}

void plot_overlap()
{
    std::string meson = get_meson();

    if(meson == "Jpsi")
        plot_Jpsi_overlap();

    else if(meson == "phi")
        plot_phi_overlap();

    else
        throw std::runtime_error("Meson desconhecido: " + meson);
}



// ------------------ PLOT SEÇÕES DE CHOQUE INTEGRAIS



void plot_sigma_Jpsi(std::string csv_file, std::string dipolemodel)
{
    int Q2 = 0;

    std::vector<int> dataset;
    std::vector<double> W_exp, sigma_exp, err_exp;

    read_sigma_exp(
    "csv/sigma_gammap_jpsi.csv",
    dataset,          
    W_exp,
    sigma_exp,
    err_exp);

    plt::figure_size(800,600);

    // estilos por experimento
    std::map<int,std::string> exp_map = {
        {0,"H1"},
        {1,"H1"},
        {2,"ALICE"},
        {3,"LHCb"}
    };

    std::map<std::string,std::string> colors = {
        {"H1","blue"},
        {"ALICE","black"},
        {"LHCb","purple"}
    };

    std::map<std::string,std::string> markers = {
        {"H1","o"},
        {"ALICE","s"},
        {"LHCb","^"}
    };

    // agrupar por dataset
    std::map<int,std::vector<int>> groups;

    for(size_t i=0;i<dataset.size();++i)
        groups[dataset[i]].push_back(i);

    for(auto& g : groups)
    {
        int d = g.first;
        std::string exp = exp_map[d];

        std::vector<double> W,s,e;

        for(int idx : g.second)
        {
            W.push_back(W_exp[idx]);
            s.push_back(sigma_exp[idx]);
            e.push_back(err_exp[idx]);
        }

        plt::errorbar(
            W,
            s,
            e,
            {{"label",exp},
             {"marker",markers[exp]},
             {"color",colors[exp]},
             {"linestyle","none"}}
        );
    }

    // --- curvas teóricas ---
    try
    {
        std::vector<double> W,sigma_GLC,sigma_BG;

        read_csv(
            csv_file,
            W,
            sigma_GLC,
            sigma_BG);

        plt::plot(W,sigma_GLC,
            {{"label","GLC"},
             {"color","red"},
             {"linewidth","1.8"}});

        plt::plot(W,sigma_BG,
            {{"label","BG"},
             {"color","red"},
             {"linestyle","--"},
             {"linewidth","1.2"}});
    }
    catch(...)
    {
        std::cerr<<"Falha ao carregar curva teórica\n";
    }

    plt::xlim(20,10000);
    plt::ylim(1e1,2e3);

    plt::xlabel("$W$ [GeV]");
    plt::ylabel("$\\sigma$ [nb]");

    plt::title("$J/\\psi$ produção exclusiva ($\\gamma p \\to J/\\psi p$)");

    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.xscale('log')\n"
        "plt.yscale('log')\n"
        "plt.grid(True, which='both', linestyle='--', alpha=0.6)\n"
        "plt.tight_layout()\n"
    );

    plt::legend();

    std::string plotname = extrair_nome_base(csv_file);
    std::string out =
        "plots/sigma/" +
        plotname + "_" +
        dipolemodel + "_" +
        timestamp() + ".png";

    std::string label = "Arquivo: " + plotname;

PyRun_SimpleString(
    ("import matplotlib.pyplot as plt\n"
     "plt.figtext(0.01, 0.01, '" + label + "', fontsize=8, alpha=0.7)\n").c_str()
);
    plt::save(out);
    plt::show();
}

void plot_sigma_phi(std::string csv, std::string dipolemodel)
{
    int Q2 = 0;
    std::string plotname = extrair_nome_base(csv);

    // --- dados teóricos ---
    std::vector<double> W, sigma_GLC, sigma_BG;
    read_csv(csv, W, sigma_GLC, sigma_BG);

    // --- dados experimentais ---
    std::vector<int> dataset_fixedpoint;
std::vector<double> W_fixedpoint, sigma_fixedpoint, error_fixedpoint;

read_sigma_exp(
    "csv/expdata/sigma/phi_fixedpoint_data(nb).csv",
    dataset_fixedpoint,     
    W_fixedpoint,
    sigma_fixedpoint,
    error_fixedpoint
);


    std::vector<int> dataset_ZEUS;
    std::vector<double> W_ZEUS, sigma_ZEUS, error_ZEUS;

    read_sigma_exp(
        "csv/expdata/sigma/phi_sigma_expdata_ZEUS(1994).csv",
        dataset_ZEUS,
        W_ZEUS,
        sigma_ZEUS,
        error_ZEUS);

    plt::figure_size(800,600);

    plt::plot(W, sigma_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(W, sigma_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}});

    plt::errorbar(W_fixedpoint, sigma_fixedpoint, error_fixedpoint,
        {{"fmt","o"},
         {"color","blue"},
         {"label","Fixed Point"}});

    plt::errorbar(W_ZEUS, sigma_ZEUS, error_ZEUS,
        {{"fmt","s"},
         {"color","green"},
         {"label","ZEUS (1994)"}});

    plt::xlim(8,100);
    plt::ylim(100.0,10000.0);

    plt::xlabel("$W$ [GeV]");
    plt::ylabel("$\\sigma$ [nb]");

    std::stringstream title;
    title << "$\\phi$ produção exclusiva ($\\gamma p \\to \\phi p$)"
          << " — $Q^2=" << Q2 << "$";

    plt::title(title.str());

    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.xscale('log')\n"
        "plt.yscale('log')\n"
        "plt.grid(True, which='both', linestyle='--', alpha=0.6)\n"
        "plt.tight_layout()\n"
    );

    plt::legend();

    std::string out =
        "plots/sigma/" + plotname + "_" +
        std::to_string(Q2) + "_" +
        dipolemodel + "_" +
        timestamp() + ".png";


    std::string label = "Arquivo: " + plotname;

PyRun_SimpleString(
    ("import matplotlib.pyplot as plt\n"
     "plt.figtext(0.01, 0.01, '" + label + "', fontsize=8, alpha=0.7)\n").c_str()
);
    plt::save(out);
    plt::show();
}





void plot_sigma(std::string meson, std::string csv, std::string dipolemodel)
{
    if(meson == "Jpsi")
        plot_sigma_Jpsi(csv, dipolemodel);
    else if(meson == "phi")
        plot_sigma_phi(csv, dipolemodel);
    else
        throw std::runtime_error("Meson desconhecido: " + meson);
}

void plot_sigma_Jpsi_gammaPb_GBW(std::string csv_file)
{
    std::vector<double> W, sigma_GLC, sigma_BG;
    read_csv(csv_file, W, sigma_GLC, sigma_BG);

    plt::clf();

    // =========================
    // Dados experimentais
    // =========================
    std::vector<double> W_CMS, sigma_CMS, error_CMS;
    read_xyscan_csv("csv/expdata/sigma/sigma_psi_gammaPb_CMS.csv",
                    W_CMS, sigma_CMS, error_CMS);

    plt::errorbar(W_CMS, sigma_CMS, error_CMS,
        {{"fmt","o"}, {"color","blue"}, {"label","CMS"}});// CMS Collaboration. Probing Small Bjorken-x Nuclear Gluonic Structure via
                                                        //Coherent J/ψ Photoproduction in Ultraperipheral Pb-Pb Collisions at √sNN = 5.02 TeV.
                                                        //Phys. Rev. Lett., American Physical Society, 2023.

    std::vector<double> W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023;
    read_xyscan_csv("csv/expdata/sigma/sigma_psi_gammaPb_ALICE_2023.csv",  
                    W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023);

    plt::errorbar(W_ALICE_2023, sigma_ALICE_2023, error_ALICE_2023,      //ALICE Collaboration. Energy dependence of coherent photonuclear production of
                                                                        //J/ψ mesons in ultra-peripheral Pb-Pb collisions at √sNN = 5.02 TeV. 2023. 
        {{"fmt","s"}, {"color","green"}, {"label","ALICE(2023)"}});

    std::vector<double> W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024;
    read_xyscan_csv("csv/expdata/sigma/sigma_psi_gammaPb_ALICE_2024.csv",
                    W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024);  // SHATAT, A. Charmonium photoproduction in Pb-Pb collisions with nuclear overlap
                                                                        //measured with ALICE at the LHC. Tese (Theses) — Université Paris-Saclay, 2024. https://theses.hal.science/tel-04797642
    
    plt::errorbar(W_ALICE_2024, sigma_ALICE_2024, error_ALICE_2024,
        {{"fmt","^"}, {"color","orange"}, {"label","ALICE(2024)"}});

    // =========================
    // Banda sombreada
    // =========================

    std::vector<double> sigma_min(W.size()), sigma_max(W.size());

for (size_t i = 0; i < W.size(); ++i) {
    sigma_min[i] = std::min(sigma_GLC[i], sigma_BG[i]);
    sigma_max[i] = std::max(sigma_GLC[i], sigma_BG[i]);
}

// converter para string Python
std::string W_py   = vec_to_pylist(W);
std::string min_py = vec_to_pylist(sigma_min);
std::string max_py = vec_to_pylist(sigma_max);

// plot banda
std::string cmd =
    "import matplotlib.pyplot as plt\n"
    "W = " + W_py + "\n"
    "y1 = " + min_py + "\n"
    "y2 = " + max_py + "\n"
    "plt.fill_between(W, y1, y2, alpha=0.3, color='orange', label='GBW band')\n";

PyRun_SimpleString(cmd.c_str());


    // =========================
    // Estética
    // =========================

    plt::xlabel("W (GeV)");
    plt::ylabel("σ(γ Pb → J/ψ Pb) [nb]");
    plt::title("Produção exclusiva de J/ψ em γ-Pb (GBW)");

    plt::xlim(25,1000);
    plt::ylim(1e3,1e5);

    plt::legend();

    PyRun_SimpleString(
        "ax = plt.gca()\n"
        "ax.set_xscale('log')\n"
        "ax.set_yscale('log')\n"
        "ax.grid(True)\n"
    );

    std::filesystem::create_directories("plots/sigma");
    plt::save("plots/sigma/Jpsi_gammaPb_GBW"+ timestamp() + ".png");

    plt::show();
}

void plot_sigma_phi_gammaPb_GBW(std::string csv_file)
{
    std::vector<double> W, sigma_GLC, sigma_BG;

    read_csv(csv_file, W, sigma_GLC, sigma_BG);

    std::cout << "DEBUG J/psi:\n";
for (int i = 0; i < 5; ++i) {
    std::cout << W[i] << "  "
              << sigma_GLC[i] << "  "
              << sigma_BG[i] << std::endl;
}
double minv = 1e100, maxv = 0;

for (double v : sigma_GLC) {
    if (std::isfinite(v) && v > 0) {
        minv = std::min(minv, v);
        maxv = std::max(maxv, v);
    }
}

std::cout << "min = " << minv << "  max = " << maxv << std::endl;

    plt::clf();
    plt::figure_size(800,600);

    plt::named_plot("GBW (GLC)", W, sigma_GLC, "-");
    plt::named_plot("GBW (Boosted Gaussian)", W, sigma_BG, "--");

    plt::xlabel("W (GeV)");
    plt::ylabel("σ(γ Pb → φ Pb) [nb]");
    plt::title("Produção exclusiva de φ em γ-Pb (GBW)");

    plt::legend();

    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gcf()\n"   // garante que usa a figura atual
    "plt.xscale('log')\n"
    "plt.yscale('log')\n"
    "plt.grid(True)\n"
);

    std::filesystem::create_directories("plots/sigma");
    plt::save("plots/sigma/phi_gammaPb_GBW_"+ timestamp() + ".png");

    plt::show();
}

void create_TA_b_csv(std::string meson)
{
    std::vector<double> W, TA_GLC, TA_BG;

    read_csv("csv/" + meson + "_TA_b.csv", W, TA_GLC, TA_BG);

    std::ofstream out("csv/" + meson + "_TA_b_processed.csv");
    out << "W,TA_GLC,TA_BG\n";
    for(size_t i=0;i<W.size();++i)
        out << W[i] << "," << TA_GLC[i] << "," << TA_BG[i] << "\n";
    out.close();
}

void plot_TA_b(){
    double R = 6.62*5.07; // fm para GeV^-1
    double a = 0.546*5.07; // fm para GeV^-1
    double rho0 = 0.1603/std::pow(5.07, 3); // fm^-3 para GeV^3 

    TA_Table table = precompute_TA(200, 50.0, R, a, rho0); // pré-calcula T_A(b) para Pb-208

    std::vector<double> b = table.b_vals;
    std::vector<double> TA = table.TA_vals;

    plt::figure_size(800,600);
    plt::plot(b, TA, {{"label","T_A(b)"}});
    plt::xlabel("b (GeV^-1)");
    plt::ylabel("T_A(b) (GeV^2)");
    plt::title("Função de perfil nuclear T_A(b) para Pb-208");
    plt::grid(true);
    plt::legend();
    plt::save("plots/TA_b(Pb-208)_" + timestamp() + ".png");
    plt::show();
    
}


void plot_rapidez_PbPb_Jpsi(std::string filename, double sqrt_s)
{
    std::vector<double> Y, rap_GLC, rap_BG;

    read_csv(filename, Y, rap_GLC, rap_BG);

    plt::figure_size(800,600);

    plt::plot(Y, rap_GLC,
        {{"label","GLC"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(Y, rap_BG,
        {{"label","BG"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}});

    plt::xlim(-8,8);

    plt::xlabel("Y");
    plt::ylabel("dσ/dY [nb]");

    double sqrt_s_TeV = sqrt_s / 1e3;

    std::stringstream title;
    title << "Distribuição de rapidez do $J/\\psi$ em Pb-Pb a " << sqrt_s_TeV << " TeV";

    plt::title(title.str());

    plt::grid(true);
    plt::legend();


    std::string out =
        "plots/Rapidez/PbPb-Jpsi_rapidez_+" + doubleParaString(sqrt_s_TeV, 3) + "TeV_" +
        timestamp() + ".png";

    plt::save(out);
    plt::show();
}

void plot_rapidez_PbPb_phi(std::string csv, double sqrt_s)
{
    std::vector<double> Y, rap_GLC, rap_BG;

    std::vector<double> Y_GBW, rap_GLC_GBW, rap_BG_GBW;

    read_csv("csv/phi_rapidez_PbPb_5.36e+03GeV.csv", Y_GBW, rap_GLC_GBW, rap_BG_GBW);
    read_csv(csv, Y, rap_GLC, rap_BG);

    plt::figure_size(800,600);

    // --- dados experimentais ---
    std::vector<double> y, dsigma_dy, error;

    if(sqrt_s == 5.36e3) // só tem dados para 5.36 TeV
    {
    read_rapidity_hepdata(
        "csv/expdata/Rapidez/rapidity_PbPb-phi_5360TeV.csv", // Observation of coherent ϕ(1020) meson photoproduction
        y,                                                   //    in ultraperipheral PbPb collisions at √sNN = 5.36 TeV
        dsigma_dy,                                           //http://dx.doi.org/10.1103/2ssw-wwyy
        error
    );


    plt::errorbar(y, dsigma_dy, error,
        {{"fmt","o"},
         {"color","blue"},
         {"label","CMS (2025)"}});

        plt::plot(Y_GBW, rap_GLC_GBW,
        {{"label","GLC GBW"},
         {"color","orange"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(Y_GBW, rap_BG_GBW,
        {{"label","BG GBW"},
         {"color","orange"},
         {"linestyle","--"},
         {"linewidth","1.2"}});
    }
    plt::plot(Y, rap_GLC,
        {{"label","GLC GBW_fc"},
         {"color","red"},
         {"linestyle","-"},
         {"linewidth","1.8"}});

    plt::plot(Y, rap_BG,
        {{"label","BG GBW_fc"},
         {"color","red"},
         {"linestyle","--"},
         {"linewidth","1.2"}});

    

    plt::ylim(10,130);
    plt::xlim(0.0,1.5);

    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gcf()\n"   // garante que usa a figura atual
    "plt.yscale('log')\n"
    "plt.grid(True)\n"
);

    plt::xlabel("Y");
    plt::ylabel("dσ/dY [mb]");

    double sqrt_s_TeV = sqrt_s / 1e3;
    std::stringstream title;
    title << "Distribuição de rapidez do $\\phi$ em Pb-Pb a " << sqrt_s_TeV << " TeV";

    plt::title(title.str());

    plt::grid(true);
    plt::legend();

    std::string out =
        "plots/Rapidez/PbPb-phi_rapidez_+" + doubleParaString(sqrt_s_TeV, 3) + "TeV_" +
        timestamp() + ".png";

    plt::save(out);
    plt::show();
}

void plot_N_dglap(std::string csv_file1, std::string csv_file2)
{
    std::vector<double> r2_x1, N_dglap_x1;
    std::vector<double> r2_x2, N_dglap_x2;

    read_two_columns(csv_file1, r2_x1, N_dglap_x1);
    read_two_columns(csv_file2, r2_x2, N_dglap_x2);


    plt::figure_size(800,600);
    plt::plot(r2_x1, N_dglap_x1, {{"label","N(r^2) DGLAP x=1e-4"}});
    plt::plot(r2_x2, N_dglap_x2, {{"label","N(r^2) DGLAP x=1e-2"}});

    PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gcf()\n"   // garante que usa a figura atual
    "plt.xscale('log')\n"
    "plt.yscale('log')\n"
    );

    plt::xlim(1e-3,100.0);
    plt::ylim(1e-4,1.2);

    plt::xlabel("$r^2$");
    plt::ylabel("$N_p$");

    plt::title("Função de dipolo N(r^2) do modelo DGLAP");

    plt::grid(true);
    plt::legend();

    std::string out =
        "plots/N/N_dglap_" +
        timestamp() + ".png";

    plt::save(out);
    plt::show();
}

void plot_overlap_fc(std::string csv_file, std::string csv_file_fc, std::string meson)
{
    std::vector<double> r1, r2;
    std::vector<double> overlap_BG, overlap_GLC;
    std::vector<double> overlap_BG_fc, overlap_GLC_fc;

    read_csv(csv_file, r1, overlap_BG, overlap_GLC);
    read_csv(csv_file_fc, r2, overlap_BG_fc, overlap_GLC_fc);

    if(r1.empty() || r2.empty()) {
    throw std::runtime_error("Vetores vazios — falha no read_csv");
}
if(meson == "Jpsi") {
    // -------- SEM FC --------

plt::plot(r1, overlap_BG, {{"label","BG"}});
plt::plot(r1, overlap_GLC, {{"label","GLC"}});

PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gca().set_xscale('log')\n"
);

plt::title("Sem fator de correção");
plt::xlabel("$r$ [fm]");
plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");


plt::xlim(0.0001, 1.0);
plt::ylim(0.0, 0.025);

plt::legend();
plt::grid(true);



plt::save("plots/overlap/sem_fc_" + meson + "_" + timestamp() + ".png");

// -------- COM FC --------
plt::figure();

plt::plot(r2, overlap_BG_fc, {{"label","BG (fc)"}});
plt::plot(r2, overlap_GLC_fc, {{"label","GLC (fc)"}});

PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gca().set_xscale('log')\n"
);

plt::title("Com fator de correção");
plt::xlabel("$r$ [fm]");
plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");
plt::legend();
plt::grid(true);

plt::xlim(0.0001, 1.0);
plt::ylim(0.0, 0.025);


plt::save("plots/overlap/com_fc_" + meson + "_" + timestamp() + ".png");

plt::show();
}
else if(meson == "phi") {
    // -------- SEM FC --------
plt::figure();

plt::plot(r1, overlap_BG, {{"label","BG"}});
plt::plot(r1, overlap_GLC, {{"label","GLC"}});
PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gca().set_xscale('log')\n"
);

plt::title("Sem fator de correção");
plt::xlabel("$r$ [fm]");
plt::xlim(0.01, 3.0);
plt::ylim(0.0, 0.01);
plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");


plt::legend();
plt::grid(true);

plt::save("plots/overlap/sem_fc_" + meson + "_" + timestamp() + ".png");

// -------- COM FC --------
plt::figure();

plt::plot(r2, overlap_BG_fc, {{"label","BG (fc)"}});
plt::plot(r2, overlap_GLC_fc, {{"label","GLC (fc)"}});

PyRun_SimpleString(
    "import matplotlib.pyplot as plt\n"
    "plt.gca().set_xscale('log')\n"
);
plt::title("Com fator de correção");
plt::xlabel("$r$ [fm]");
plt::ylabel("Overlap $r \\Psi_V \\Psi_{\\gamma}$");
plt::xlim(0.01, 3.0);
plt::ylim(0.0, 0.01);


plt::legend();
plt::grid(true);

plt::save("plots/overlap/com_fc_" + meson + "_" + timestamp() + ".png");

plt::show();
}
}


















































