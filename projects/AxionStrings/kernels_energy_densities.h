#ifndef PROJECTS_AXION_STRINGS_KERNELS_ENERGY_DENSITIES_H_
#define PROJECTS_AXION_STRINGS_KERNELS_ENERGY_DENSITIES_H_

#include "setup.h"

namespace AxionStrings {

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

}; // namespace AxionStrings

#endif // PROJECTS_AXION_STRINGS_KERNELS_ENERGY_DENSITIES_H_
