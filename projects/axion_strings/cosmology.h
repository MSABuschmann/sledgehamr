#ifndef PROJECTS_AXION_STRINGS_COSMOLOGY_H_
#define PROJECTS_AXION_STRINGS_COSMOLOGY_H_

#include <sledgehamr.h>

namespace axion_strings {

ADD_SCALARS(Psi1, Psi2)
ADD_CONJUGATE_MOMENTA(Pi1, Pi2)

/** @brief TODO
 */
AMREX_FORCE_INLINE
double a_prime2(amrex::Array4<amrex::Real const> const& state, const int i,
        const int j, const int k, const int lev, const double time,
        const double dt, const double dx,
        const std::vector<double>& params) {
    double Psi1    = state(i, j, k, Scalar::Psi1);
    double Psi2    = state(i, j, k, Scalar::Psi2);
    double Pi1     = state(i, j, k, Scalar::Pi1);
    double Pi2     = state(i, j, k, Scalar::Pi2);
    double r2      = Psi1*Psi1 + Psi2*Psi2;
    double prime_a = (Psi1*Pi2 - Psi2*Pi1)/r2;
    return prime_a*prime_a;
}

AMREX_FORCE_INLINE
double a_prime_screened(amrex::Array4<amrex::Real const> const& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx,
        const std::vector<double>& params) {
    double Psi1    = state(i, j, k, Scalar::Psi1);
    double Psi2    = state(i, j, k, Scalar::Psi2);
    double Pi1     = state(i, j, k, Scalar::Pi1);
    double Pi2     = state(i, j, k, Scalar::Pi2);
    double prime_a = (Psi1*Pi2 - Psi2*Pi1);
    return prime_a;
}

AMREX_FORCE_INLINE
double r_prime2(amrex::Array4<amrex::Real const> const& state, const int i,
        const int j, const int k, const int lev, const double time,
        const double dt, const double dx,
        const std::vector<double>& params) {
    double Psi1    = state(i, j, k, Scalar::Psi1);
    double Psi2    = state(i, j, k, Scalar::Psi2);
    double Pi1     = state(i, j, k, Scalar::Pi1);
    double Pi2     = state(i, j, k, Scalar::Pi2);
    double r2      = Psi1*Psi1 + Psi2*Psi2;
    double prime_r = (Psi1*Pi1 + Psi2*Pi2);
    return prime_r*prime_r/r2;
}

/** @brief Class to simulate axion strings.
 */
class Cosmology {
  public:
    void Init(sledgehamr::Sledgehamr* owner);

    bool CreateLevelIf(const int lev, const double time);

    double Mr(const double eta) {
        return std::sqrt(2. * lambda) * eta;
    }

    double H(const double eta) {
        return 1./eta;
    }

    double StringWidth(const int lev, const double eta) {
        return 1./(Mr(eta) * sim->GetDx(lev));
    }

    double RefinementTime(const int lev) {
        return sim->GetDimN(lev) / (sqrt(2.*lambda) 
                * string_width_threshold * sim->GetL());
    }

    double Log(double eta) {
        if (eta <= 0)
            return -DBL_MAX;
        return std::log( Mr(eta) / H(eta) );
    }

    double LogTruncated(const double eta) {
        double log = Log(eta);
        if (log < spectra_log_min) 
            return 0;
        else
            return log;
    }

  private:
    void ParseVariables();
    void PrintRefinementTimes();
    void SetProjections();
    void SetSpectra();

    double string_width_threshold;
    double spectra_log_min = 5;
    const double lambda = 1;

    sledgehamr::Sledgehamr* sim;
};

}; // namespace axion_strings

#endif // PROJECTS_AXION_STRINGS_COSMOLOGY_H_
