#ifndef PROJECTS_AXION_STRINGS_KERNELS_TAGGING_H_
#define PROJECTS_AXION_STRINGS_KERNELS_TAGGING_H_

#include "setup.h"

namespace AxionStrings {

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
};

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

}; // namespace AxionStrings

#endif // PROJECTS_AXION_STRINGS_KERNELS_TAGGING_H_
