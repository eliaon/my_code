/*
 * IPsat dipole amplitude fitted to HERA data
 *
 * Reference: H. Mäntysaari, P. Zurita,
 * Phys.Rev. D98 (2018) 036002 (arXiv:1804.05311)
 *
 * This code solves DGLAP using the fitted parametrization
 * for the initial condition, supporting both non-linear IPsat
 * and linearized IPnonsat (fitted separately)
 *
 * Questions and comments: Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "dipoleamplitude.hpp"
#include <string>
#include <cmath>
#include <iostream>

using namespace MZ_ipsat;

DipoleAmplitude::DipoleAmplitude(IPSAT_PARAMETRIZATION mode)
{
    if (mode == MZ_IPSAT)
    {
        C = 2.2894; lambda_g = 0.08289; A_g = 2.1953; mc = 1.3528;
        mu0=std::sqrt(1.1);
        saturation=true;
    }
    else if (mode == MZ_IPNONSAT)
    {
        C = 4.2974; lambda_g = -0.006657; A_g = 3.0391; mc = 1.3504;
        mu0=std::sqrt(1.1);
        saturation=false;
    }
    
    B_p = 4.0;
    Nc = 3;
    mb = 4.75;
    mt = 175;
    
    // Init alphas(M_Z=91.1876 GeV) = 0.1183
    alphas = new AlphaStrong(0, 1.0, 91.1876, 0.1183, mc, mb, mt);
    // DGLAP_Solver will take care of deleting alphas when it is deleted
    cppdglap = new EvolutionLO_gluon(alphas);
}

DipoleAmplitude::DipoleAmplitude(double C_, double mu0_, double lambda_g_, double A_g_, double mc_, double mb_, double mt_)
{
    C=C_; mu0=mu0_; lambda_g = lambda_g_; A_g = A_g_; mc = mc_; mb = mb_; mt = mt_;
    saturation = true;
    B_p=4.0;    // All fits are done with fixed B_p=4
    Nc=3;
    mb=4.75;
    mt=175;
    
    // Init alphas(M_Z=91.1876 GeV) = 0.1183
    alphas = new AlphaStrong(0, 1.0, 91.1876, 0.1183, mc, mb, mt);
    // DGLAP_Solver will take care of deleting alphas when it is deleted
    cppdglap = new EvolutionLO_gluon(alphas);
}

/*
 * Initialize lookup table
 */
void DipoleAmplitude::EnableLookupTable()
{
    cppdglap->generateLookupTable(mu0, 0, A_g, lambda_g, 0, 0);
    cppdglap->useLookupTable(true);
}

void DipoleAmplitude::DisableLookupTable()
{
    cppdglap->useLookupTable(false);
}

DipoleAmplitude::~DipoleAmplitude()
{
    delete cppdglap;
}

double DipoleAmplitude::Alphas_xg(double x, double musqr)
{
    double as_xg=0;
    int coupling = 0;
    double As=0;
    double lambdas=0; //singlet
    
    return cppdglap->alphasxG(x, musqr, mu0, 0, A_g, lambda_g, As, lambdas);
}

double DipoleAmplitude::Alphas(double Q)
{
    return alphas->value(Q);
}

double DipoleAmplitude::xg(double x, double musqr)
{
    return cppdglap->xG(x, musqr, mu0, 0, A_g, lambda_g, 0 , 0);
}

double DipoleAmplitude::N(double r, double xbj, double b)
{
    double musqr = mu0*mu0 + C / (r*r);
    double exponent = M_PI*M_PI / (2.0 * Nc) * r*r * Alphas_xg(xbj, musqr) * Tp(b);
    
    if (!saturation)     // IPnonsat
        return exponent;
    else
        return 1.0 - std::exp(-exponent);

}



double DipoleAmplitude::Tp(double b)
{
    return 1.0/(2.0*M_PI*B_p) * std::exp(-b*b / (2.0*B_p));
}


