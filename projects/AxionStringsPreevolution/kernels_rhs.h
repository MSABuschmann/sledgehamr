#ifndef PROJECTS_AXION_STRINGS_PREEVOLUTION_KERNELS_RHS_H_
#define PROJECTS_AXION_STRINGS_PREEVOLUTION_KERNELS_RHS_H_

#include "../AxionStrings/cosmology.h"
#include <sledgehamr_utils.h>

namespace AxionStringsPreevolution {
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
Rhs(const amrex::Array4<double> &rhs, const amrex::Array4<const double> &state,
    const int i, const int j, const int k, const int lev, const double time,
    const double dt, const double dx, const double *params) {
    // Fetch field values.
    double Psi1 = state(i, j, k, Scalar::Psi1);
    double Psi2 = state(i, j, k, Scalar::Psi2);
    double Pi1 = state(i, j, k, Scalar::Pi1);
    double Pi2 = state(i, j, k, Scalar::Pi2);

    double eta = time;

    // Compute Laplacians.
    constexpr int order = 2;
    double laplacian_Psi1 = sledgehamr::utils::Laplacian<order>(
        state, i, j, k, Scalar::Psi1, dx * dx);
    double laplacian_Psi2 = sledgehamr::utils::Laplacian<order>(
        state, i, j, k, Scalar::Psi2, dx * dx);

    // Compute EOM.
    const double eta_0 = params[0];
    const double eta_sq = eta * eta;
    const double lambda = eta_0 * eta_0;
    const double drag = std::sqrt(eta_0);
    // const double drag = eta_0;

    double potential = lambda / eta_sq * (Psi1 * Psi1 + Psi2 * Psi2 - 1.);

    rhs(i, j, k, Scalar::Psi1) = Pi1;
    rhs(i, j, k, Scalar::Psi2) = Pi2;
    rhs(i, j, k, Scalar::Pi1) =
        -Pi1 * 3. / eta + laplacian_Psi1 / eta_sq / drag - Psi1 * potential;
    rhs(i, j, k, Scalar::Pi2) =
        -Pi2 * 3. / eta + laplacian_Psi2 / eta_sq / drag - Psi2 * potential;
}

}; // namespace AxionStringsPreevolution

#endif // PROJECTS_AXION_STRINGS_PREEVOLUTION_KERNELS_RHS_H_
