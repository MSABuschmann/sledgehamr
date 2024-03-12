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
 * @param   state   Current state.
 * @param   i       i-th cell index.
 * @param   j       j-th cell index.
 * @param   k       k-th cell index.
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
#define AXION_STRING_TRUNCATION_MODIFIER(x) \
template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE \
double TruncationModifier<x>(const amrex::Array4<const double>& state, \
        const int i, const int j, const int k, const int lev, \
        const double time, const double dt, const double dx, \
        const double truncation_error, const double* params) { \
    return truncation_error * dt; \
}

AXION_STRING_TRUNCATION_MODIFIER(Scalar::Pi1)
AXION_STRING_TRUNCATION_MODIFIER(Scalar::Pi2)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_xx)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_yy)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_zz)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_xy)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_xz)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_yz)

}; // namespace AxionStrings

#endif // PROJECTS_AXION_STRINGS_KERNELS_TAGGING_H_
