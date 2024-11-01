#ifndef SLEDGEHAMR_INTEGRATOR_LEAPFROG_H_
#define SLEDGEHAMR_INTEGRATOR_LEAPFROG_H_

#include "integrator.h"

namespace sledgehamr {

/** @brief Implementation of the leap-frog integration scheme (kick-drift-kick).
 */
class IntegratorLeapfrog : public Integrator {
    using Integrator::Integrator;

  protected:
    virtual void Integrate(LevelData &mf_old, LevelData &mf_new, const int lev,
                           const double dt, const double dx) override;

  private:
    void DebugPrint(amrex::MultiFab &mf, const char *msg);
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_INTEGRATOR_LEAPFROG_H_
