#pragma once

#include <sledgehamr.h>

namespace BetaDecay {

SLEDGEHAMR_ADD_SCALARS(Phi)
SLEDGEHAMR_ADD_CONJUGATE_MOMENTA(dPhi)

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
Rhs(const amrex::Array4<double> &rhs, const amrex::Array4<const double> &state,
    const int i, const int j, const int k, const int lev, const double time,
    const double dt, const double dx, const double *params) {}

SLEDGEHAMR_FINISH_SETUP

class BetaDecay : public sledgehamr::Sledgehamr {
  public:
    SLEDGEHAMR_INITIALIZE_PROJECT(BetaDecay)
};

}; // namespace BetaDecay
