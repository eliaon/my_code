#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <ctime>
#include <Python.h>
#include "matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;



// ---------- FUNÇÃO PARA PUXAR A STRING TIMESTAMP (HORÁRIO E DIA/MES/ANO AGORA)

std::string timestamp(void)
{
    std::time_t now = std::time(nullptr);
    std::tm* local = std::localtime(&now);

    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", local);
    return std::string(buffer);
}
// ------------ LÊ CSV COM 2 COLUNAS, PARA N 
void read_two_columns(
    const std::string& filename,
    std::vector<double>& x,
    std::vector<double>& y)
{
    std::ifstream file(filename);
    std::cout << "Reading file: " << filename << std::endl;
    if(!file)
        throw std::runtime_error("Cannot open file: " + filename);

    std::string line;

    std::getline(file,line); // header

    while(std::getline(file,line))
    {
        if(line.empty()) continue;

        std::stringstream ss(line);
        std::string a,b;

        std::getline(ss,a,',');
        std::getline(ss,b,',');

        x.push_back(std::stod(a));
        y.push_back(std::stod(b));
    }
}

// ------------ LÊ CSV'S COM 1 INDEPENDENTE E DUAS DEPENDENTES, IDEAL PARA COMPARAÇÕES
// ------------- BOOSTED GAUSSIAN - GAUSSIAN LIGHT CONE

void read_csv(
    const std::string& filename,
    std::vector<double>& x,
    std::vector<double>& glc,
    std::vector<double>& bg)
{
    std::ifstream file(filename);

    if(!file)
        throw std::runtime_error("Cannot open file: " + filename);

    std::string line;

    std::getline(file, line); // header
    
    while (std::getline(file, line))
    {
        if(line.empty()) continue;

        std::stringstream ss(line);
        std::string a,b,c;

        std::getline(ss, a, ',');
        std::getline(ss, b, ',');
        std::getline(ss, c, ',');

        x.push_back(std::stod(a));
        glc.push_back(std::stod(b));
        bg.push_back(std::stod(c));
    }
}

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

// -------------- CARREGA SET DE VALORES DE CURVAS DE N

bool load_set(
    const std::string& xstr,
    std::vector<double>& r2_ipsat,
    std::vector<double>& N_ipsat,
    std::vector<double>& r2_gbw,
    std::vector<double>& N_gbw,
    std::vector<double>& r2_iim,
    std::vector<double>& N_iim)
{
    std::cout << "Loading x=" << xstr << std::endl;
    read_two_columns("csv/N_ipsat_x=" + xstr + ".csv", r2_ipsat, N_ipsat);
    read_two_columns("csv/N_GBW_x=" + xstr + ".csv",  r2_gbw,  N_gbw);
    read_two_columns("csv/N_IIM_x=" + xstr + ".csv",  r2_iim,  N_iim);
    if(r2_ipsat.empty() || r2_gbw.empty() || r2_iim.empty())
{
    std::cerr << "Empty dataset for x=" << xstr << std::endl;
    return false;
}
return true;
}

// -------------- PLOTA AS CURVAS DE N

void plot_N_models()
{
    
    std::vector<std::string> xvals = {"1e-4","1e-2"};

    for(const auto& xstr : xvals)
    {
        plt::figure_size(700,500);

        std::vector<double> r2_ipsat,N_ipsat;
        std::vector<double> r2_gbw,N_gbw;
        std::vector<double> r2_iim,N_iim;

        load_set(xstr,
                 r2_ipsat,N_ipsat,
                 r2_gbw,N_gbw,
                 r2_iim,N_iim);

        plt::plot(r2_ipsat,N_ipsat, {{"label","IPsat"}});
        plt::plot(r2_gbw,N_gbw, {{"label","GBW"},{"linestyle","--"}});
        plt::plot(r2_iim,N_iim, {{"label","IIM"},{"linestyle",":"}});

        plt::title("x=" + xstr);

        plt::xlim(0.0,100.0);
        plt::ylim(0.0,1.2);

        plt::xlabel("$r^2$");
        plt::ylabel("$N_p$");

        plt::grid(true);
        plt::legend();

        std::string filename =
            "plots/N/N_gbw_iim_ipsat_x=" + xstr + "_" + timestamp() + ".png";

        plt::save(filename);
        plt::show();
        
    }
}

// -------------- PLOTA CURVAS DE RAPIDEZ PARA UM SQRT S ESCOLHIDO

std::string get_meson()
{
    std::cout << "Escreva o meson (Jpsi, phi): ";
    std::string meson;
    std::cin >> meson;

    if (meson != "Jpsi" && meson != "phi")
    {
        std::cerr << "Meson inválido. Use 'Jpsi' ou 'phi'." << std::endl;
        exit(1);
    }

    return meson;
}

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

void read_sigma_exp(
    const std::string& filename,
    std::vector<int>& dataset,
    std::vector<double>& W,
    std::vector<double>& sigma,
    std::vector<double>& error)
{
    std::ifstream file(filename);
    if(!file)
        throw std::runtime_error("Cannot open file: " + filename);

    std::string line;
    std::getline(file,line); // header

    while(std::getline(file,line))
    {
        if(line.empty()) continue;

        std::stringstream ss(line);
        std::string a,b,c,d;

        std::getline(ss,a,',');
        std::getline(ss,b,',');
        std::getline(ss,c,',');
        std::getline(ss,d,',');

        dataset.push_back(std::stoi(a));
        W.push_back(std::stod(b));
        sigma.push_back(std::stod(c));
        error.push_back(std::stod(d));
    }
}

void plot_sigma_Jpsi()
{
    int Q2 = 0;

    std::vector<int> dataset;
    std::vector<double> W_exp, sigma_exp, err_exp;

    read_sigma_exp(
        "csv/sigma_gammap_jpsi.csv",
        dataset, W_exp, sigma_exp, err_exp);

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
            "csv/Jpsi_sigma_Q2=0.csv",
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

    std::string out =
        "plots/sigma/sigma_Jpsi_Q2=0_" +
        timestamp() + ".png";

    plt::save(out);
    plt::show();
}

void plot_sigma_phi()
{
    int Q2 = 0;

    std::string filename =
        "csv/phi_sigma_Q2=" + std::to_string(Q2) + ".csv";

    std::vector<double> W, sigma_GLC, sigma_BG;

    read_csv(filename, W, sigma_GLC, sigma_BG);

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

    // ponto experimental H1
    std::vector<double> W_exp = {70};
    std::vector<double> sigma_exp = {960};

    plt::scatter(W_exp, sigma_exp, 30.0,
        {{"color","blue"},{"label","H1"}});

    plt::xlim(30,1000);
    plt::ylim(10,10000);

    plt::xlabel("$W$ [GeV]");
    plt::ylabel("$\\sigma$ [nb]");

    std::stringstream title;
    title << "$\\phi$ produção exclusiva ($\\gamma p \\to \\phi p$)"
          << " — $Q^2=" << Q2 << "$";

    plt::title(title.str());

    // escalas log
    PyRun_SimpleString(
        "import matplotlib.pyplot as plt\n"
        "plt.xscale('log')\n"
        "plt.yscale('log')\n"
        "plt.grid(True, which='both', linestyle='--', alpha=0.6)\n"
        "plt.tight_layout()\n"
    );

    plt::legend();

    std::string out =
        "plots/sigma/sigma_phi_Q2=" +
        std::to_string(Q2) + "_" +
        timestamp() + ".png";

    plt::save(out);
    plt::show();
}

void plot_sigma(std::string meson)
{
    if(meson == "Jpsi")
        plot_sigma_Jpsi();
    else if(meson == "phi")
        plot_sigma_phi();
    else
        throw std::runtime_error("Meson desconhecido: " + meson);
}

