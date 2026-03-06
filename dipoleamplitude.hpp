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

#include "dglap_cpp/AlphaStrong.h"
#include "dglap_cpp/EvolutionLO_nocoupling.h"

namespace MZ_ipsat
{
    
enum IPSAT_PARAMETRIZATION
{
    MZ_IPSAT,       // Standard IPsat
    MZ_IPNONSAT     // Linearized IPsat
};

class DipoleAmplitude
{
    public:
        /*
         * Constructors: user can specify his own initial condition, or use
         * fixed fit results (MZ_IPSAT and MZ_IPNONSAT)
         */
        DipoleAmplitude(IPSAT_PARAMETRIZATION mode=MZ_IPSAT);
        DipoleAmplitude(double C_, double mu0_, double lambda_g_, double A_g_, double mc_, double mb_=4.75, double mt_=175 );
        ~DipoleAmplitude();
    
        /*
         * Enable lookup table
         * Faster evaluation of the dipole, but not really exact
         */
        void EnableLookupTable();
    
        /*
         * Disable lookup table
         */
        void DisableLookupTable();
    
        /*
         * Compute alphas(mu^2) * xg(x, mu^2), musqr in GeV^2
         */
        double Alphas_xg(double x, double musqr);
    
        // Strong coupling, Q in GeV
        double Alphas(double Q);
    
        // Only xg, calculated internally as Alphas_xg / Alphas, so not as effective!
        double xg(double x, double musqr);
    
        /*
         * Dipole amplitude
         * [r] = GeV^-1, dipole size
         * [b]=Gev^-1, impact parameter
         */
        double N(double r, double xbj, double b);

    
        /*
         * Proton profile, normalized to unity
         * \int d^2 b T_p = 1
         *
         * [b] = GeV^-1
         */
        double Tp(double b);
    
        // Other methods
        double GetMu0() { return mu0; }
        double GetLambdaG() { return lambda_g; }
        double GetAg() { return A_g; }
        double GetMc() { return mc; }
        double GetMb() { return mb; }
        double GetMt() { return mt; }
        bool GetSaturation() { return saturation; }
        void SetSaturation(bool s) { saturation = s; }
    
    private:
        EvolutionLO_gluon *cppdglap;
        AlphaStrong *alphas;
        bool saturation;   // True for ipsat, false for ipnonsat
        double InitAlphas();
        double C;
        double mu0;
        double lambda_g;
        double A_g;
        double mc;
        double mb;
        double mt;
        double B_p; // Proton size, GeV^-2
        double alphas_mur;      // alpha_s at mu0
        int Nc;
};
    
};
