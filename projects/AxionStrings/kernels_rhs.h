#ifndef PROJECTS_AXION_STRINGS_KERNELS_RHS_H_
#define PROJECTS_AXION_STRINGS_KERNELS_RHS_H_

#include <sledgehamr_utils.h>
#include "setup.h"

namespace AxionStrings {

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
    double potential = eta*eta*( Psi1*Psi1 + Psi2*Psi2 - 1. ) + 0.56233;

    rhs(i, j, k, Scalar::Psi1) =  Pi1;
    rhs(i, j, k, Scalar::Psi2) =  Pi2;
    rhs(i, j, k, Scalar::Pi1)  = -Pi1*2./eta + laplacian_Psi1 - Psi1*potential;
    rhs(i, j, k, Scalar::Pi2)  = -Pi2*2./eta + laplacian_Psi2 - Psi2*potential;
}

template<> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void GravitationalWavesRhs<true>(const amrex::Array4<double>& rhs,
        const amrex::Array4<const double>& state, const int i, const int j,
        const int k, const int lev, const double time, const double dt,
        const double dx, const double* params) {
    // Fetch field values.
    double du_xx = state(i, j, k, Gw::du_xx);
    double du_yy = state(i, j, k, Gw::du_yy);
    double du_zz = state(i, j, k, Gw::du_zz);
    double du_xy = state(i, j, k, Gw::du_xy);
    double du_xz = state(i, j, k, Gw::du_xz);
    double du_yz = state(i, j, k, Gw::du_yz);

    double eta = time;

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
    double grad_x_Psi1 = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Psi1, dx, 'x');
    double grad_y_Psi1 = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Psi1, dx, 'y');
    double grad_z_Psi1 = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Psi1, dx, 'z');

    double grad_x_Psi2 = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Psi2, dx, 'x');
    double grad_y_Psi2 = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Psi2, dx, 'y');
    double grad_z_Psi2 = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Psi2, dx, 'z');

    // Compute EOM.
    rhs(i, j, k, Gw::u_xx) = du_xx;
    rhs(i, j, k, Gw::u_yy) = du_yy;
    rhs(i, j, k, Gw::u_zz) = du_zz;
    rhs(i, j, k, Gw::u_xy) = du_xy;
    rhs(i, j, k, Gw::u_xz) = du_xz;
    rhs(i, j, k, Gw::u_yz) = du_yz;

    rhs(i, j, k, Gw::u_xx)  = - du_xx*2./eta + laplacian_u_xx
                              + grad_x_Psi1*grad_x_Psi1
                              + grad_x_Psi2*grad_x_Psi2;
    rhs(i, j, k, Gw::u_yy)  = - du_yy*2./eta + laplacian_u_yy
                              + grad_y_Psi1*grad_y_Psi1
                              + grad_y_Psi2*grad_y_Psi2;
    rhs(i, j, k, Gw::u_zz)  = - du_zz*2./eta + laplacian_u_zz
                              + grad_z_Psi1*grad_z_Psi1
                              + grad_z_Psi2*grad_z_Psi2;
    rhs(i, j, k, Gw::u_xy)  = - du_xy*2./eta + laplacian_u_xy
                              + grad_x_Psi1*grad_y_Psi1
                              + grad_x_Psi2*grad_y_Psi2;
    rhs(i, j, k, Gw::u_xz)  = - du_xz*2./eta + laplacian_u_xz
                              + grad_x_Psi1*grad_z_Psi1
                              + grad_x_Psi2*grad_z_Psi2;
    rhs(i, j, k, Gw::u_yz)  = - du_yz*2./eta + laplacian_u_yz
                              + grad_y_Psi1*grad_z_Psi1
                              + grad_y_Psi2*grad_z_Psi2;
}

}; // namespace AxionStrings

#endif  // PROJECTS_AXION_STRINGS_KERNELS_RHS_H_
