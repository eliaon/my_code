// integration.hpp
#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include <functional>
#include <cmath>  // para std::abs

// ---------------- Integração numérica (Simpson) ----------------
double integrate_simpson(const std::function<double(double)>& f,
                         double a, double b, int n = 1000);

// ---------------- Gauss-Legendre 8-pontos (SGS8) ----------------
double sgs8(double A, double B, const std::function<double(double)>& F);

// ---------------- Integração adaptativa (SGS0) ----------------
double sgs0(double A, double B, double EPS, const std::function<double(double)>& F);

#endif // INTEGRATION_HPP