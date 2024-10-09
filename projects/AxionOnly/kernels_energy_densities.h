#ifndef PROJECTS_AXION_ONLY_KERNELS_ENERGY_DENSITIES_H_
#define PROJECTS_AXION_ONLY_KERNELS_ENERGY_DENSITIES_H_

#include "setup.h"

namespace AxionOnly {

/** @brief Kernel function computing axion energy density a'^2.
 * @param   state   Data from which to calculate RHS (current state).
 * @param   i       i-th cell index.
 * @param   j       j-th cell index.
 * @param   k       k-th cell index.
 * @param   lev     Current level.
 * @param   time    Current time.
 * @param   dt      Time step size.
 * @param   dx      Grid spacing.
 * @param   params  Optional parameters.
 * @return a'^2.
 */
AMREX_FORCE_INLINE
double theta_prime2(amrex::Array4<amrex::Real const> const &state, const int i,
                    const int j, const int k, const int lev, const double time,
                    const double dt, const double dx,
                    const std::vector<double> &params) {
    double dtheta = state(i, j, k, Scalar::dtheta);
    return dtheta * dtheta;
}

}; // namespace AxionOnly

#endif // PROJECTS_AXION_ONLY_KERNELS_ENERGY_DENSITIES_H_
