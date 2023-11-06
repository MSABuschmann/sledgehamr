#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_RHS_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_RHS_H_

#include <sledgehamr_utils.h>
#include "setup.h"

namespace FirstOrderPhaseTransition {

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
    double quadratic = params[0];
    double cubic     = params[1];
    double quartic   = params[2];
    double Phi       = state(i, j, k, Scalar::Phi);
    double potential = quadratic*Phi + cubic*Phi*Phi + quartic*Phi*Phi*Phi;

    constexpr int order = 2;
    double laplacian_Phi = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Scalar::Phi, dx*dx);

    rhs(i, j, k, Scalar::Phi)  = state(i, j, k, Scalar::dPhi);
    rhs(i, j, k, Scalar::dPhi) = laplacian_Phi + potential;
}

template<> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void GravitationalWavesRhs<true>(const amrex::Array4<double>& rhs,
        const amrex::Array4<const double>& state, const int i, const int j,
        const int k, const int lev, const double time, const double dt,
        const double dx, const double* params) {
    // Compute Laplacians.
    double dx2 = dx * dx;
    constexpr int order = 2;
    double laplacian_u_xx = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_xx, dx2);
    double laplacian_u_yy = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_yy, dx2);
    double laplacian_u_zz = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_zz, dx2);
    double laplacian_u_xy = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_xy, dx2);
    double laplacian_u_xz = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_xz, dx2);
    double laplacian_u_yz = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_yz, dx2);

    // Compute gradients.
    constexpr int order2 = 2;
    double grad_x_Phi = sledgehamr::utils::Gradient<order2>(
            state, i, j, k, Scalar::Phi, dx, 'x');
    double grad_y_Phi = sledgehamr::utils::Gradient<order2>(
            state, i, j, k, Scalar::Phi, dx, 'y');
    double grad_z_Phi = sledgehamr::utils::Gradient<order2>(
            state, i, j, k, Scalar::Phi, dx, 'z');

    // Compute EOM.
    rhs(i, j, k, Gw::u_xx) = state(i, j, k, Gw::du_xx);
    rhs(i, j, k, Gw::u_yy) = state(i, j, k, Gw::du_yy);
    rhs(i, j, k, Gw::u_zz) = state(i, j, k, Gw::du_zz);
    rhs(i, j, k, Gw::u_xy) = state(i, j, k, Gw::du_xy);
    rhs(i, j, k, Gw::u_xz) = state(i, j, k, Gw::du_xz);
    rhs(i, j, k, Gw::u_yz) = state(i, j, k, Gw::du_yz);
    rhs(i, j, k, Gw::du_xx) = laplacian_u_xx + grad_x_Phi*grad_x_Phi;
    rhs(i, j, k, Gw::du_yy) = laplacian_u_yy + grad_y_Phi*grad_y_Phi;
    rhs(i, j, k, Gw::du_zz) = laplacian_u_zz + grad_z_Phi*grad_z_Phi;
    rhs(i, j, k, Gw::du_xy) = laplacian_u_xy + grad_x_Phi*grad_y_Phi;
    rhs(i, j, k, Gw::du_xz) = laplacian_u_xz + grad_x_Phi*grad_z_Phi;
    rhs(i, j, k, Gw::du_yz) = laplacian_u_yz + grad_y_Phi*grad_z_Phi;
}

}; // namespace FirstOrderPhaseTransition

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_RHS_H_
