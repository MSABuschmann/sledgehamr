#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_TAGGING_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_TAGGING_H_

#include "setup.h"

namespace FirstOrderPhaseTransition {

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

AXION_STRING_TRUNCATION_MODIFIER(Scalar::dPhi)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_xx)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_yy)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_zz)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_xy)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_xz)
AXION_STRING_TRUNCATION_MODIFIER(Gw::du_yz)

}; // namespace FirstOrderPhaseTransition

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_TAGGING_H_
