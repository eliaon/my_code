#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <vector>
#include <iomanip>
#include <sstream>
#include <boost/math/special_functions/bessel.hpp>

// ---------------- Constantes globais ----------------
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double alfem = 1.0 / 137.0;   // constante de estrutura fina
const double Nc    = 3.0;           // número de cores
const double lambda = 0.29;       // parâmetro do modelo GBW
const double gamma_s = 0.46; // parâmetro do modelo bCGC
const double CFAC = 5.07; // conversão fm para GeV^-1

// Conversão de unidades
const double GeV2_to_nb = 3.89379e5; // 1 GeV^-2 = 3.89379×10^5 nb
// Parâmetros de seção de choque
//const double sigma0_GeV2 = 0.23; // Constante sigma0 em GeV^-2. (23.0 mb = 0.23 GeV^-2)
// ---------------- Classe Meson ----------------
class Meson {
public:
    std::string meson;
    std::string nome;
    double MV;   // massa do méson (GeV)
    double mf;   // massa do quark (GeV)
    double ef;   // carga efetiva
    double NT;   // normalização transversa
    double NL;   // normalização longitudinal
    double R2T;  // parâmetro transverso GLC
    double R2L;  // parâmetro longitudinal GLC
    double R2;   // parâmetro BG
    bool isGLC;  // se é GLC ou BG

    // Construtor Gaus-LC
    Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
          double NT_, double R2T_, double NL_, double R2L_)
        : meson(m), nome(n), MV(MV_), mf(mf_), ef(ef_), NT(NT_), NL(NL_),
          R2T(R2T_), R2L(R2L_), R2(0.0), isGLC(true) {
            //if m = Jpsi {
            //    B = 4.5 + 2* 0.164 *log(W2/95);
           //     else if m = phi {
            //        B = 0.60 * (14.0 / pow((Q2 + MV*MV), 0.26) + 1.0);
           // }
          }

    // Construtor Boosted Gaussian
    Meson(std::string m, std::string n, double MV_, double mf_, double ef_,
          double NT_, double NL_, double R2_)
        : meson(m), nome(n), MV(MV_), mf(mf_), ef(ef_), NT(NT_), NL(NL_),
          R2T(0.0), R2L(0.0), R2(R2_), isGLC(false) {}
};

// ---------------- função para calcular B(Q2) ----------------
inline double calc_B(double x, double Q2, const Meson& M) {
    double W2 = M.MV*M.MV/x;
    if (M.meson == "Jpsi"){
        double B1 = 4.5 + 2* 0.164 *log(W2/95);
        return B1;
    }
        else if (M.meson == "phi"){
            double B2 = 0.60 * (14.0 / pow((Q2 + M.MV*M.MV), 0.26) + 1.0);
            return B2;}
            else {
                std::cerr << "Mésen desconhecido para cálculo de B: " << M.meson << std::endl;
                return 0.0;
            }
}
// ---------------- phi_T e derivada ----------------
inline double phi_T(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        double zz = z * (1.0 - z);
        return M.NT * (zz * zz) * exp(-(r * r) / (2.0 * M.R2T));
    } else {
        double zz = z * (1.0 - z);
        double part1 = M.NT * zz;
        double arg1  = -(M.mf*M.mf * M.R2) / (8.0 * zz);
        double arg2  = -(2.0 * zz * r * r) / M.R2;
        double arg3  = (M.mf*M.mf * M.R2) / 2.0;
        return part1 * exp(arg1 + arg2 + arg3);
    }
}

inline double dphi_dr(double r, double z, const Meson& M)
{
    if (M.isGLC) {
        return -(r / M.R2T) * phi_T(r, z, M);
    } else {
        double zz = z * (1.0 - z);
        return -((4.0 * r * zz) / M.R2) * phi_T(r, z, M);
    }
}

// ----------------psi_V psi_T ----------------
inline double psi_Vpsi_T(double z, double r,double Q2,  const Meson& M)
{
    double Mf2 = M.mf * M.mf;
    double EPS2 = z * (1.0 - z) * Q2 + Mf2;
    double EPS  = sqrt(EPS2);

    double K0 = boost::math::cyl_bessel_k(0, EPS * r);
    double K1 = boost::math::cyl_bessel_k(1, EPS * r);

    double PHIT  = phi_T(r, z, M);
    double DPHIT = dphi_dr(r, z, M);

    double ZZZ = z*z + (1.0 - z)*(1.0 - z);
    double ANORM = M.ef * sqrt(4.0 * M_PI * alfem) * Nc / (M_PI * z * (1.0 - z));

    return ANORM * ((Mf2 * K0 * PHIT) - (ZZZ * EPS * K1 * DPHIT));
}

// ----------------psi_V*psi_L ----------------
inline double psi_Vpsi_L(double z, double r,double Q2,  const Meson& M)
{
    double Mf2 = M.mf * M.mf;
    double EPS2 = z * (1.0 - z) * Q2 + Mf2;
    double EPS  = sqrt(EPS2);

    double K0 = boost::math::cyl_bessel_k(0, EPS * r);

    double PHIL  = phi_T(r, z, M);

    double ANORM = M.ef * sqrt(4.0 * M_PI * alfem) * Nc / (M_PI * z * (1.0 - z));

    return ANORM * (2.0 * sqrt(Q2) * z * (1.0 - z) * K0 * PHIL);
}

// ---------------- integração numérica (Simpson) ----------------
double integrate_simpson(const std::function<double(double)>& f,
                         double a, double b, int n = 1000)
{
    if (n < 2) n = 2;
    if (n % 2 != 0) ++n;
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        const double x = a + i * h;
        s += f(x) * ((i % 2 == 0) ? 2.0 : 4.0);
    }
    return s * h / 3.0;
}
// Implementação das rotinas de integração SGS0 e SGS8 (Gauss-Legendre 8-point)
double sgs8(double A, double B, const std::function<double(double)>& F) {
    // regras e pontos/weights conforme Fortran
    double C = 0.5 * (A + B);
    double H = 0.5 * (B - A);
    double S = 0.0;
    double X, tmp;

    X = 0.96028985 * H;
    S = 0.10122853 * (F(C + X) + F(C - X));
    X = 0.79666647 * H;
    S += 0.22238103 * (F(C + X) + F(C - X));
    X = 0.52553240 * H;
    S += 0.31370664 * (F(C + X) + F(C - X));
    X = 0.18343464 * H;
    S += 0.36268378 * (F(C + X) + F(C - X));
    return S * H;
}

double sgs0(double A, double B, double EPS, const std::function<double(double)>& F) {
    double S = 0.0;
    double U = A;
    // Loop: enquanto U < B, processar um subintervalo [U,V]
    while (U < B) {
        double V = B;
        // primeiro SF = SGS8(F,U,V)
        double SF = sgs8(U, V, F);
        bool refine = true;
        while (refine) {
            double C = 0.5 * (U + V);
            double SL = sgs8(U, C, F);
            double SG = sgs8(C, V, F);
            double SP = SL + SG;
            double abb = abs(SF);
            if (abb == 0.0) abb = 1.0;
            double SA = abs(SP - SF) / (abb * EPS);
            if (SA >= 1.0) {
                // precisa refinar o intervalo [U,V] -> reduzir V para C e repetir
                V = C;
                SF = SL;
                // continuar refinando
            } else {
                // aceitável
                refine = false;
                // somar SP ao resultado e pular para o próximo subintervalo
                S += SP;
                U = V;
            }
            // caso V tenha chegado a U por aproximação numérica, evitar loop infinito
            if (V - U <= 1e-15) {
                // adiciona SF e avança para evitar travamento
                S += SF;
                U = V;
                refine = false;
            }
        } // end refine
    } // end while U<B
    return S;
}

// ---------------- overlap integrado em z ----------------
double overlap_r(double r, double Q2, const Meson& M, int Nz = 200) {
    auto fz = [r, Q2, &M](double z) {
        return psi_Vpsi_T(z, r, Q2, M) + psi_Vpsi_L(z, r, Q2, M);
    };
    return integrate_simpson( fz, 1e-6, 1.0 - 1e-6, Nz );
}
//---------------- QS2 --------------------

double QS2_bCGC(double x, double b,
            double x0 = 1.84e-6, double Q0sq = 1.0)
{
    double b2 = b * b;
    double B_bCGC = 7.5;
    double bexp = std::exp(-b2 / (2.0 *B_bCGC)); // B=5.0 GeV^-2
    return Q0sq * pow(x0 / x, lambda)*bexp;
}

double QS2_GBW(double x, double x_0, double Q0sq=1.0)
{
    return Q0sq * pow(x_0 / x, lambda);
}   
// ---------------- Modelo GBW ----------------
double N_GBW(double r, double x,
             double x0 = 3e-4)
{
    double sigma0 = 20.3* 0.003894; // gev^-2
    double Qs2 = QS2_GBW(x,  x0);
    double arg = (r * r) * Qs2 / 4.0;
    return sigma0*(1.0 - exp(-arg));
}
// ---------------- N_bCGC ----------------
// ---------------- N_bCGC ----------------
double N_bCGC(double r, double x, double b,
              double x0 = 1.84e-6, double lambda_bCGC = 0.119)
{
    double Qs2 = QS2_bCGC(x, b, x0, 1.0);
    double N_0 = 0.558;
    double c = N_0*gamma_s/(1-N_0);
    double A = -c*c/log(1.0-N_0);
    double B = 0.5*std::pow(1.0-N_0, -1.0/c);
    double rQs = r * sqrt(Qs2);
    double Y = log(1.0/x);
    double K = 9.9;

    if (rQs <= 2.0) {
        double term1 = N_0*pow(rQs*rQs/4.0, 2.0*(gamma_s + log(2.0/rQs)/(K*lambda_bCGC*Y)));
        return term1;
    } else {
        
        return 1.0 - exp(-A * log(B * rQs * rQs)*log(B * rQs * rQs) );
    }
}
//----------------- plot n_bCGC -----------------------
void N_plot(std::string N_method)
{   
    if (N_method != "bCGC" && N_method != "GBW") {
        std::cerr << "Método N desconhecido: " << N_method << std::endl;
        return;
    }

    if (N_method == "bCGC") {
        double x = 1e-4;
        auto Nrb = [&](double r, double b) {
            return N_bCGC(r, x, b);
        };
        auto Nr_integrated =  [&](double r) {
            // integrar em b
            auto fb = [&](double b) {
                return Nrb(r, b) * 2.0 * M_PI * b; // incluir fator 2πb do jacobiano
            };
            // integrar em b de 0 a 2 GeV^-1 usando 200 subdivisões
            double Ib = integrate_simpson(fb, 0.0, 2.0, 200);
            return Ib;
        };
        std::string filename = "csv/N_bCGC_r.csv";
        std::ofstream fout(filename);
        fout << "r,N_bCGC\n";
        const int Npoints = 1000;
        double rmin = 1e-4, rmax = 100.0;

        for (int i = 0; i < Npoints; ++i) {
            double frac = static_cast<double>(i) / (Npoints - 1);
            double r = rmin * pow(rmax / rmin, frac);

            double N_bCGC = Nr_integrated(r);
            fout << r*r << "," << N_bCGC << "\n";
        }
        fout.close();
        std::cout << "Arquivo '" << filename << "' gerado." << std::endl;
        return;
    }

    if (N_method == "GBW") {
        double x = 1e-4;
        // For GBW the local amplitude does not depend on impact parameter b,
        // so call N_GBW with (r,x) and ignore b in the lambda.
        auto Nrb = [&](double r, double b) {
            (void)b;
            return N_GBW(r, x);
        };
        auto Nr_integrated =  [&](double r) {
            // integrar em b
            auto fb = [&](double b) {
                return Nrb(r, b) * 2.0 * M_PI * b; // incluir fator 2πb do jacobiano
            };
            // integrar em b de 0 a 2 GeV^-1 usando 200 subdivisões
            double Ib = integrate_simpson(fb, 0.0, 2.0, 200);
            return Ib;
        };
        std::string filename = "csv/N_gbw_r.csv";
        std::ofstream fout(filename);
        fout << "r,N_gbw\n";
        const int Npoints = 1000;
        double rmin = 1e-4, rmax = 100.0;

        for (int i = 0; i < Npoints; ++i) {
            double frac = static_cast<double>(i) / (Npoints - 1);
            double r = rmin * pow(rmax / rmin, frac);

            double N_gbw = Nr_integrated(r);
            fout << r*r << "," << N_gbw << "\n";
        }
        fout.close();
        std::cout << "Arquivo '" << filename << "' gerado." << std::endl;
        return;
    }
}
// ---------------- Rg --------------------------

const double Rg = (pow(2, 2*lambda+3.0)*tgamma(lambda + 2.5))/(sqrt(M_PI)*tgamma(lambda + 4.0));

// ---------------- densidade de gluons efetiva ----------------
double xg(double x, double mu2)
{
    // Parametrização simples para a densidade de gluons
    // xg(x, mu2) = A_g * x^{-lambda_g} * (1 - x)^5.6
    const double A_g = 1.0;      // Normalização
    const double lambda_g = 0.2; // Exponente de x

    return A_g * pow(x, -lambda_g) * pow(1 - x, 5.6);
}

// ---------------- proton shape function ----------------
double T(double b)
{
    const double B_G = 4.0; // GeV^-2
    return (1.0 / (2.0 * M_PI * B_G)) * exp(-b * b / (2.0 * B_G));
}


// ---------------- Amplitude (Integral) ----------------
double calculate_amplitude(double x, double Q2, std::string N_method, const Meson& M,
                           int Nr = 600, int Nz = 200,
                           double rmin = 1e-4, double rmax = 10.0)
{
    const double bmax = 10.0; // GeV^-1

    if (N_method == "bCGC") {
        auto fr = [&](double r) {
            auto fb = [&](double b) {
                double Ov = overlap_r(r, Q2, M, Nz);
                double Nrval = N_bCGC(r, x, b);
                return Ov * Nrval * b; // <-- REMOVIDO 2π aqui
            };
            double Ib = integrate_simpson(fb, 0.0, bmax, 200);
            return Ib * r; // <-- agora inclui r de d²r = 2π r dr
        };

        double Ir = integrate_simpson(fr, rmin, rmax, Nr);
        return 2.0 * M_PI * Ir; // <-- 2π total de d²r
    }
    else if (N_method == "GBW") {
        auto fr = [&](double r) {
            double Ov = overlap_r(r, Q2, M, Nz);
            double Nrval = N_GBW(r, x);
            return Ov * Nrval * r; // r de d²r = 2π r dr
        };
        double Ir = integrate_simpson(fr, rmin, rmax, Nr);
        return 2.0 * M_PI * Ir;
    }
    else {
        std::cerr << "Método N desconhecido: " << N_method << std::endl;
        return 0.0;
    }
}
double calculate_deltinha(double x, double Q2, const Meson& M, const std::string& N_method,
                          int Nr = 600, int Nz = 200,
                          double rmin = 1e-4, double rmax = 50.0)
{
    const double epsilon = 1e-3; // passo relativo

    double x_plus  = x * (1.0 + epsilon);
    double x_minus = x * (1.0 - epsilon);

    // Evitar x fora do domínio físico
    if (x_minus <= 0) x_minus = x * 0.5;

    // Calcular amplitudes
    double A0   = calculate_amplitude(x, Q2, N_method, M, Nr, Nz, rmin, rmax);
    double A_p  = calculate_amplitude(x_plus, Q2, N_method, M,  Nr, Nz, rmin, rmax);
    double A_m  = calculate_amplitude(x_minus, Q2, N_method, M,  Nr, Nz, rmin, rmax);

    // Usar módulo para evitar log de número negativo
    A0 = std::abs(A0);
    A_p = std::abs(A_p);
    A_m = std::abs(A_m);

    // Evitar log(0)
    if (A0 <= 0 || A_p <= 0 || A_m <= 0) {
        std::cerr << "Atenção: amplitude não positiva em x=" << x << std::endl;
        return lambda; // fallback
    }

    // Derivada numérica central: d ln A / d ln(1/x) = - d ln A / d ln x
    // Mas como ln(1/x) = -ln x, então:
    // delta = d ln A / d ln(1/x) = - d ln A / d ln x
    // Usando diferenças finitas:
    // delta ≈ [ln A(x_minus) - ln A(x_plus)] / [ln(x_plus) - ln(x_minus)]

    double logA_m = std::log(A_m);
    double logA_p = std::log(A_p);
    double logx_p = std::log(x_plus);
    double logx_m = std::log(x_minus);

    double delta = (logA_m - logA_p) / (logx_p - logx_m);

    // Garantir valor físico
    if (delta < 0) delta = lambda;
    if (std::isnan(delta) || std::isinf(delta)) delta = lambda;

    return delta;
}
// ---------------- sigma(x)----------------
// Esta função retorna a seção de choque em GeV^-2.
double sigma_x(double x, double Q2 , const Meson& M, std::string N_method,
               int Nr = 600, int Nz = 200, 
               double rmin = 1e-4, double rmax = 10.0)
{
    // A seção de choque é proporcional ao quadrado da amplitude.
    double amplitude = calculate_amplitude(x, Q2, N_method , M,  Nr, Nz, rmin, rmax);
    std::cout << "Amplitude calculada: " << amplitude << std::endl;
    double B_val = calc_B(x, Q2,  M); // Unidades de GeV^-2
    std::cout << "B calculado: " << B_val << std::endl;
    double deltinha = calculate_deltinha(x, Q2, M, N_method, Nr, Nz, rmin, rmax);
    //std::cout << "Deltinha calculado: " << deltinha << std::endl;
    double termo_RG = std::pow(2.0,2.0*deltinha +3.0) * tgamma(deltinha + 2.5) 
                        / (sqrt(M_PI) * tgamma(deltinha + 4.0));
    //std::cout << "Termo Rg calculado: " << termo_RG << std::endl;
    double beta = tan(M_PI * deltinha / 2.0);
    //std::cout << "Beta calculado: " << beta << std::endl;
    double fator_beta = 1 + beta * beta;
    //std::cout << "Fator beta calculado: " << fator_beta << std::endl;
    std::cout << "quadrado da amplitude: " << amplitude * amplitude << std::endl;

    return (amplitude * amplitude) / (16.0 * M_PI * B_val)
           * termo_RG * termo_RG
           //* pow(1-x, 5.26)
           * fator_beta;
           //* (1 + tan(M_PI * deltinha / 2) * tan(M_PI * deltinha / 2));
};

void calculate_sigma(double Q2 , const Meson& M_GLC, const Meson& M_BG, std::string N_method) {
    std::ostringstream oss;
    oss <<"csv/"<< M_GLC.meson << "_sigma_Q2=" << std::fixed << std::setprecision(1) << Q2 << ".csv";
    std::string filename = oss.str();


        // Opcional: remover zeros finais para valores como 3.2
        size_t pos = filename.find_last_not_of('0');
        if (pos != std::string::npos && filename[pos] == '.') {
            filename.erase(pos, 1); // Remove o '.' e o '0'
        } else if (pos != std::string::npos && filename.back() == '0' && filename[pos] != '.') {
            // Caso como 3.20, remove o zero
            filename.erase(pos + 1);
        }

        std::ofstream fout(filename);
        fout << "W,sigma_GLC,sigma_BG\n";

        const int Nx = 40;
        double xmin = 1e-4, xmax = 1e-1;
        for (int i = 0; i < Nx; ++i) {
            double frac = static_cast<double>(i) / (Nx - 1);
            double x = xmin * pow(xmax / xmin, frac);

            double W = sqrt((M_GLC.MV*M_BG.MV)/x);

            // resultado em GeV^-2 → converter para nb
            double sigma_glc = sigma_x(x, Q2, M_GLC, N_method) * GeV2_to_nb;
            double sigma_bg  = sigma_x(x, Q2, M_BG, N_method)  * GeV2_to_nb;

            fout << W << "," << sigma_glc << "," << sigma_bg << "\n";
            std::cout << "Q2=" << Q2 << "  x="<<x<<"  W="<<W<<" GeV  sigma_GLC="<<sigma_glc
                      <<" nb  sigma_BG="<<sigma_bg<<" nb\n";
        }

        fout.close();
        std::cout << "Arquivo '" << filename << "' gerado." << std::endl;
}

void plot_overlap(const Meson& M_GLC, const Meson& M_BG, std::string N_method, int Nz = 200)
{
    std::string filename = M_GLC.meson + "_overlap_r.csv";
    std::ofstream fout(filename);
    fout << "r,overlap_GLC,overlap_BG\n";
    const int Npoints = 1000;
    double rmin = 1e-4, rmax = 10.0;
    double Q2 = 0.0;

    for (int i = 0; i < Npoints; ++i) {
        double frac = static_cast<double>(i) / (Npoints - 1);
        double r = rmin * pow(rmax / rmin, frac);

        double overlap_glc = overlap_r(r, Q2, M_GLC, Nz);
        double overlap_bg  = overlap_r(r, Q2, M_BG, Nz);

        fout << r/CFAC << "," << overlap_glc << "," << overlap_bg << "\n";
    }
    fout.close();
    std::cout << "Arquivo '" << filename << "' gerado." << std::endl;
}
// ---------------- main ----------------
int main()
{
    double qJ = 2.0/3.0;
    double qS = 1.0/3.0;
    Meson Jpsi_GLC("Jpsi", "Jpsi_GLC", 3.097, 1.4, qJ, 1.23, 6.5, 0.0, 0.0);
    Meson Jpsi_BG ("Jpsi", "Jpsi_BG",  3.097, 1.4, qJ, 0.578, 0.575, 2.3);
    Meson phi_GLC("phi", "phi_GLC", 1.019, 0.14, qS, 4.75, 16.0, 0.0, 0.0);
    Meson phi_BG("phi", "phi_BG", 1.019, 0.14, qS, 0.919, 0.825, 11.2);

    double Q2 = 0;
    //N_plot("GBW");
    //plot_overlap(Jpsi_GLC, Jpsi_BG, "GBW");
    //calculate_sigma(Q2, phi_GLC, phi_BG, "bCGC");
    double x = 1e-3;
    std::cout << "Calculando sigma para x=" << x << " ou W=" << sqrt((Jpsi_GLC.MV*Jpsi_BG.MV)/x) << " e Q2=" << Q2 << " GeV^2" << std::endl;
    double sigma = sigma_x(x, Q2, Jpsi_BG, "GBW");
    std::cout << "Sigma calculada: " << sigma << " GeV^-2" << std::endl;
    std::cout << "Sigma calculada: " << sigma * GeV2_to_nb << " nb" << std::endl;

    return 0;
}