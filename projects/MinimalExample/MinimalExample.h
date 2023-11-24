#pragma once

#include <sledgehamr.h>

namespace MinimalExample {

SLEDGEHAMR_ADD_SCALARS(Psi1, Psi2)
SLEDGEHAMR_ADD_CONJUGATE_MOMENTA(Pi1, Pi2)

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Rhs(const amrex::Array4<double>& rhs,
         const amrex::Array4<const double>& state,
         const int i, const int j, const int k, const int lev,
         const double time, const double dt, const double dx,
         const double* params) {
    // Fetch field values.
    double Psi1 = state(i, j, k, Scalar::Psi1);
    double Psi2 = state(i, j, k, Scalar::Psi2);
    double Pi1  = state(i, j, k, Scalar::Pi1);
    double Pi2  = state(i, j, k, Scalar::Pi2);

    double eta = time;

    // Compute Laplacians.
    constexpr int order = 2;
    double laplacian_Psi1 = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Scalar::Psi1, dx*dx);
    double laplacian_Psi2 = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Scalar::Psi2, dx*dx);

    // Compute potential.
    double potential = eta*eta*( Psi1*Psi1 + Psi2*Psi2 - 1. ) + 0.56233;

    // Full EOM.
    rhs(i, j, k, Scalar::Psi1) =  Pi1;
    rhs(i, j, k, Scalar::Psi2) =  Pi2;
    rhs(i, j, k, Scalar::Pi1)  = -Pi1*2./eta + laplacian_Psi1 - Psi1*potential;
    rhs(i, j, k, Scalar::Pi2)  = -Pi2*2./eta + laplacian_Psi2 - Psi2*potential;
}

SLEDGEHAMR_FINISH_SETUP

class MinimalExample : public sledgehamr::Sledgehamr {
  public:
    SLEDGEHAMR_INITIALIZE_PROJECT(MinimalExample)
};

}; // namespace MinimalExample
