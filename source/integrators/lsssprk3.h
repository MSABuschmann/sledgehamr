#ifndef INTEGRATOR_LSSSPRK3_H_
#define INTEGRATOR_LSSSPRK3_H_

#include "integrator.h"

namespace sledgehamr {

class IntegratorLSSSPRK3 : public Integrator {
    using Integrator::Integrator;

  protected:
    virtual void Integrate(LevelData& mf_old, LevelData& mf_new, const int lev,
                           const double dt, const double dx) override;
};

}; // namespace sledgehamr

#endif // INTEGRATOR_LSSSPRK3_H_