// integration.cpp
#include "integration.hpp"
#include <cmath>  // std::abs

double integrate_simpson(const std::function<double(double)>& f,
                         double a, double b, int n)
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

double sgs8(double A, double B, const std::function<double(double)>& F) {
    double C = 0.5 * (A + B);
    double H = 0.5 * (B - A);
    double S = 0.0;

    double X = 0.96028985 * H;
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
    const double MIN_STEP = 1e-15;

    while (U < B) {
        double V = B;
        double SF = sgs8(U, V, F);
        bool refine = true;

        while (refine) {
            double C = 0.5 * (U + V);
            double SL = sgs8(U, C, F);
            double SG = sgs8(C, V, F);
            double SP = SL + SG;
            double abb = std::abs(SF);
            if (abb == 0.0) abb = 1.0;
            double SA = std::abs(SP - SF) / (abb * EPS);

            if (SA >= 1.0) {
                V = C;
                SF = SL;
            } else {
                S += SP;
                U = V;
                refine = false;
            }

            if (V - U <= MIN_STEP) {
                S += SF;
                U = V;
                refine = false;
            }
        }
    }
    return S;
}