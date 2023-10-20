#ifndef PROJECTS_AXION_STRINGS_H_
#define PROJECTS_AXION_STRINGS_H_

#include <sledgehamr.h>
#include <sledgehamr_utils.h>
#include "cosmology.h"

namespace axion_strings {

//ADD_SCALARS(Psi1, Psi2)
//ADD_CONJUGATE_MOMENTA(Pi1, Pi2)

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
 
/** @brief Checks for zero-crossings between two points in the complex scalar
 *         field.
 * @param   Psi1_1  \Psi_1 of first point.
 * @param   Psi2_1  \Psi_2 of first point.
 * @param   Psi1_2  \Psi_1 of second point.
 * @param   Psi2_2  \Psi_2 of second point.
 * @return Sign of slope of zero-crossing. 0 if no crossing.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ZeroXing(double Psi1_1, double Psi2_1, double Psi1_2, double Psi2_2) {
    if (Psi2_1 * Psi2_2 >= 0) return 0;
    if (Psi2_1 * Psi1_2 - Psi1_1 * Psi2_2 > 0) return 1;
    return -1;
}

/** @brief Computes the winding factor along a given axis. Will be non-zero if
 *         plaquette is pierced by a string.
 * @return  Winding factor.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis1(const amrex::Array4<const double>& state,
                 const int i, const int j, const int k) {
    return ZeroXing(state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2),
                    state(i+1,j  ,k  ,Scalar::Psi1),
                    state(i+1,j  ,k  ,Scalar::Psi2))
         + ZeroXing(state(i+1,j  ,k  ,Scalar::Psi1),
                    state(i+1,j  ,k  ,Scalar::Psi2),
                    state(i+1,j+1,k  ,Scalar::Psi1),
                    state(i+1,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state(i+1,j+1,k  ,Scalar::Psi1),
                    state(i+1,j+1,k  ,Scalar::Psi2),
                    state(i  ,j+1,k  ,Scalar::Psi1),
                    state(i  ,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state(i  ,j+1,k  ,Scalar::Psi1),
                    state(i  ,j+1,k  ,Scalar::Psi2),
                    state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis2(const amrex::Array4<const double>& state,
                 const int i, const int j, const int k) {
    return ZeroXing(state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2),
                    state(i+1,j  ,k  ,Scalar::Psi1),
                    state(i+1,j  ,k  ,Scalar::Psi2))
         + ZeroXing(state(i+1,j  ,k  ,Scalar::Psi1),
                    state(i+1,j  ,k  ,Scalar::Psi2),
                    state(i+1,j  ,k+1,Scalar::Psi1),
                    state(i+1,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state(i+1,j  ,k+1,Scalar::Psi1),
                    state(i+1,j  ,k+1,Scalar::Psi2),
                    state(i  ,j  ,k+1,Scalar::Psi1),
                    state(i  ,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state(i  ,j  ,k+1,Scalar::Psi1),
                    state(i  ,j  ,k+1,Scalar::Psi2),
                    state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WindingAxis3(const amrex::Array4<const double>& state,
                 const int i, const int j, const int k) {
    return ZeroXing(state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2),
                    state(i  ,j+1,k  ,Scalar::Psi1),
                    state(i  ,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state(i  ,j+1,k  ,Scalar::Psi1),
                    state(i  ,j+1,k  ,Scalar::Psi2),
                    state(i  ,j+1,k+1,Scalar::Psi1),
                    state(i  ,j+1,k+1,Scalar::Psi2))
         + ZeroXing(state(i  ,j+1,k+1,Scalar::Psi1),
                    state(i  ,j+1,k+1,Scalar::Psi2),
                    state(i  ,j  ,k+1,Scalar::Psi1),
                    state(i  ,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state(i  ,j  ,k+1,Scalar::Psi1),
                    state(i  ,j  ,k+1,Scalar::Psi2),
                    state(i  ,j  ,k  ,Scalar::Psi1),
                    state(i  ,j  ,k  ,Scalar::Psi2));
}

/** @brief Function that tags individual cells for refinement.
 * @return  Boolean value as to whether cell should be refined or not.
 */
template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool TagCellForRefinement<true>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double* params) {
    // Check all three plaquettes (in positive index direction) for string
    // piercings.
    if (WindingAxis1(state, i, j, k) != 0) return true;
    if (WindingAxis2(state, i, j, k) != 0) return true;
    if (WindingAxis3(state, i, j, k) != 0) return true;

    return false;
}

/** @brief Modifies the truncation error criteria for Pi1 and Pi2 from its
 *         default \tau > \tau_{crit} to \tau * \Delta t_{\ell} > \tau_{crit}.
 * @param   truncation_error    \tau
 * @return f(\tau) for criteria f(\tau) > \tau_{crit}.
 */
template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Scalar::Pi1>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Scalar::Pi2>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_xx>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_yy>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_zz>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_xy>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_xz>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_yz>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

FINISH_SLEDGEHAMR_SETUP

/** @brief Class to simulate axion strings.
 */
class axion_strings : public sledgehamr::Sledgehamr {
  public:
    START_PROJECT(axion_strings)

    void Init() override;
    bool CreateLevelIf(const int lev, const double time) override;

  private:
    Cosmology cosmo;
};

}; // namespace axion_strings

#endif // PROJECTS_AXION_STRINGS_H_
