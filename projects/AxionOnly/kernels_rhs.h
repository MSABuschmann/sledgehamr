#ifndef PROJECTS_AXION_ONLY_KERNELS_RHS_H_
#define PROJECTS_AXION_ONLY_KERNELS_RHS_H_

#include "setup.h"
#include <sledgehamr_utils.h>

namespace AxionOnly {

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
 * @param   params  Optional parameters.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
Rhs(const amrex::Array4<double> &rhs, const amrex::Array4<const double> &state,
    const int i, const int j, const int k, const int lev, const double time,
    const double dt, const double dx, const double *params) {
    // Fetch field values.
    double theta = state(i, j, k, Scalar::theta);
    double dtheta = state(i, j, k, Scalar::dtheta);
    double eta = time;

    // Compute Laplacians.
    constexpr int order = 2;
    double laplacian_theta = sledgehamr::utils::Laplacian<order>(
        state, i, j, k, Scalar::theta, dx * dx);

    // Compute EOM.
    double eta_grow = params[0];
    double potential = eta_grow * eta * eta * std::sin(theta);

    rhs(i, j, k, Scalar::theta) = dtheta;
    rhs(i, j, k, Scalar::dtheta) =
        -dtheta * 2. / eta + laplacian_theta - potential;
}

}; // namespace AxionOnly

#endif // PROJECTS_AXION_ONLY_KERNELS_RHS_H_
