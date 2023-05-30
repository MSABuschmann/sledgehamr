#ifndef AMREX_INTEGRATOR_H_
#define AMREX_INTEGRATOR_H_

#include "integrator.h"

namespace sledgehamr {

class AMReXIntegrator : public Integrator {
  protected:
     virtual void Integrate(LevelData& mf_old, LevelData& mf_new, amrex::Geometry& geom, const int lev, const double dt);
};

}; // namespace sledgehamr

#endif // AMREX_INTEGRATOR_H_
