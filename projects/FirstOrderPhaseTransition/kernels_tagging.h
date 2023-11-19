#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_TAGGING_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_TAGGING_H_

#include "setup.h"

namespace FirstOrderPhaseTransition {

/** @brief Modifies the truncation error criteria for Pi1 and Pi2 from its
 *         default \tau > \tau_{crit} to \tau * \Delta t_{\ell} > \tau_{crit}.
 * @param   truncation_error    \tau
 * @return f(\tau) for criteria f(\tau) > \tau_{crit}.
 */
#define FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER(x) \
template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE \
double TruncationModifier<x>(const amrex::Array4<const double>& state, \
        const int i, const int j, const int k, const int lev, \
        const double time, const double dt, const double dx, \
        const double truncation_error, const double* params) { \
    return truncation_error / params[x]; \
}

#define FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER_DT(x) \
template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE \
double TruncationModifier<x>(const amrex::Array4<const double>& state, \
        const int i, const int j, const int k, const int lev, \
        const double time, const double dt, const double dx, \
        const double truncation_error, const double* params) { \
    return truncation_error * dt / params[x]; \
}

FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER(Scalar::Phi)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER(Gw::u_xx)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER(Gw::u_yy)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER(Gw::u_zz)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER(Gw::u_xy)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER(Gw::u_xz)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER(Gw::u_yz)

FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER_DT(Scalar::dPhi)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER_DT(Gw::du_xx)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER_DT(Gw::du_yy)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER_DT(Gw::du_zz)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER_DT(Gw::du_xy)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER_DT(Gw::du_xz)
FIRST_ORDER_PHASE_TRANSITION_TRUNCATION_MODIFIER_DT(Gw::du_yz)

}; // namespace FirstOrderPhaseTransition

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_TAGGING_H_
