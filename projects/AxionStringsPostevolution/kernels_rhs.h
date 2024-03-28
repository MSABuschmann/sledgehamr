#ifndef PROJECTS_AXION_STRINGS_POSTEVOLUTION_KERNELS_RHS_H_
#define PROJECTS_AXION_STRINGS_POSTEVOLUTION_KERNELS_RHS_H_

#include <sledgehamr_utils.h>
#include "../AxionStrings/cosmology.h"

namespace AxionStringsPostevolution {
using namespace AxionStrings;

/** @brief Function that calculates the RHS of the EOM at a single cell.
 * @param   rhs     Container to be filled with RHS.
 * @param   state   Data from which to calculate RHS (current state).
 * @param   i       i-th cell index.
 * @param   j       j-th cell index.
 * @param   k       k-th cell index.
 * @param   lev     Current level.
 * @param   time    Current time.
 * @param   dt      Time step size.
 * @param   dx      Grid spacing.
 */
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

    // Compute EOM.
    const double eta_0 = params[0];
    const double eta_pre = params[1];
    const double eta_sq = eta_pre*eta_pre;
    const double lambda = eta_0*eta_0;
    const double drag = eta_0;

    double potential_pre = lambda/eta_sq*( Psi1*Psi1 + Psi2*Psi2 - 1. );
    double rhs_Psi1_pre = -Pi1*3./eta_pre + laplacian_Psi1/eta_sq/drag
                                 - Psi1*potential_pre;
    double rhs_Psi2_pre = -Pi2*3./eta_pre + laplacian_Psi2/eta_sq/drag
                                 - Psi2*potential_pre;

    double potential_post = eta*eta*( Psi1*Psi1 + Psi2*Psi2 - 1. );
    double rhs_Psi1_post = -Pi1*2./eta + laplacian_Psi1 - Psi1*potential_post;
    double rhs_Psi2_post = -Pi2*2./eta + laplacian_Psi2 - Psi2*potential_post;

    double frac = params[2];
    rhs(i, j, k, Scalar::Psi1) =  Pi1;
    rhs(i, j, k, Scalar::Psi2) =  Pi2;
    rhs(i, j, k, Scalar::Pi1)  = (1.-frac)*rhs_Psi1_pre + frac*rhs_Psi1_post;
    rhs(i, j, k, Scalar::Pi2)  = (1.-frac)*rhs_Psi2_pre + frac*rhs_Psi2_post;
}

}; // namespace AxionStringsPostevolution

#endif // PROJECTS_AXION_STRINGS_POSTEVOLUTION_KERNELS_RHS_H_
