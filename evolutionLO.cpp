#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>

// ====================== CONSTANTES ======================
const double Nc = 3.0;
const double Nf = 4.0;
const double Lambda2 = 0.04; // (GeV^2) ~ (0.2 GeV)^2

// ====================== ALPHA_S (1-loop) ======================
double alpha_s(double mu2)
{
    double beta0 = (11.0 * Nc - 2.0 * Nf) / 3.0;
    return 4.0 * M_PI / (beta0 * std::log(mu2 / Lambda2));
}

// ====================== SPLITTING FUNCTION P_gg ======================
double Pgg(double z)
{
    double eps = 1e-6;
    if (z > 1.0 - eps) z = 1.0 - eps; // regularização

    return 6.0 * ( (1.0 / (1.0 - z)) + (1.0 / z) - 2.0 + z * (1.0 - z) );
}

// ====================== DISTRIBUIÇÃO INICIAL ======================
double xg_initial(double x)
{
    // parametrização típica: xg = A x^{-lambda} (1-x)^n
    double A = 1.0;
    double lambda = 0.3;
    double n = 4.0;

    return A * std::pow(x, -lambda) * std::pow(1.0 - x, n);
}

// ====================== CONVOLUÇÃO ======================
double convolution(double x, double mu2,
                   const std::function<double(double)>& xg)
{
    int Nz = 200;
    double zmin = x;
    double zmax = 1.0;

    double dz = (zmax - zmin) / Nz;
    double sum = 0.0;

    for (int i = 0; i < Nz; ++i)
    {
        double z = zmin + (i + 0.5) * dz;

        double xp = x / z;
        if (xp > 1.0) continue;

        sum += (1.0 / z) * Pgg(z) * xg(xp);
    }

    return sum * dz;
}

// ====================== EVOLUÇÃO EM MU² ======================
void evolve(std::vector<double>& xgrid,
            std::vector<double>& xg_vals,
            double mu2_initial,
            double mu2_final,
            int Nsteps)
{
    double t0 = std::log(mu2_initial);
    double t1 = std::log(mu2_final);
    double dt = (t1 - t0) / Nsteps;

    for (int step = 0; step < Nsteps; ++step)
    {
        double mu2 = std::exp(t0 + step * dt);

        std::vector<double> xg_new = xg_vals;

        for (size_t i = 0; i < xgrid.size(); ++i)
        {
            double x = xgrid[i];

            auto interp = [&](double xp) {
                // interpolação linear simples
                for (size_t j = 0; j < xgrid.size() - 1; ++j)
                {
                    if (xp >= xgrid[j] && xp <= xgrid[j+1])
                    {
                        double t = (xp - xgrid[j]) / (xgrid[j+1] - xgrid[j]);
                        return xg_vals[j] * (1 - t) + xg_vals[j+1] * t;
                    }
                }
                return 0.0;
            };

            double conv = convolution(x, mu2, interp);

            double deriv = (alpha_s(mu2) / (2.0 * M_PI)) * conv;

            xg_new[i] += deriv * dt;
        }

        xg_vals = xg_new;
    }
}
