#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_MISC_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_MISC_H_

#include "setup.h"
#include "bubbles.h"

namespace FirstOrderPhaseTransition {

/** @brief TODO
 */
AMREX_FORCE_INLINE
double dPhi2(amrex::Array4<double const> const& state, const int i,
        const int j, const int k, const int lev, const double time,
        const double dt, const double dx,
        const std::vector<double>& params) {
    double dPhi = state(i, j, k, Scalar::dPhi);
    return dPhi * dPhi;
}

AMREX_FORCE_INLINE
double Distance(const double a, const double b, const double L) {
    const double dist = fabs(a-b);
    return dist < L/2. ? dist : L - dist;
}

AMREX_FORCE_INLINE
void AddBubble(const int i, const int j, const int k, const double dx,
               const double L, amrex::Array4<amrex::Real> const& fab,
               const Bubble& bubble) {
    const double Dx = Distance(i*dx,bubble.x, L);
    const double Dy = Distance(j*dx,bubble.y, L);
    const double Dz = Distance(k*dx,bubble.z, L);

    double D = std::sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
    double pos = bubble.GetPos(D);
    if( pos == -1 ) {
        return;
    }
    int ind = pos;
    double frac = pos - (double)ind;

    for (int n = 0; n < fab.nComp(); ++n) {
        int n0 = n;
        if (n == Scalar::Phi) n0 = 0;
        else if (n == Scalar::dPhi) n0 = 1;
        else continue;

        double val = bubble.GetVal(n0, ind, frac);
        fab(i, j, k, n) += val;
    }
}

}; // namespace FirstOrderPhaseTransition

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_KERNELS_MISC_H_
