#ifndef SLEDGEHAMR_INTEGRATOR_LSSSPRK3_H_
#define SLEDGEHAMR_INTEGRATOR_LSSSPRK3_H_

#include "integrator.h"

namespace sledgehamr {

/** @brief Implementation of the low-storage strong stability preserving third
 *         order Runge-Kutta integration scheme (LSSSPRK3).
 */
class IntegratorLsssprk3 : public Integrator {
    using Integrator::Integrator;

  protected:
    virtual void Integrate(LevelData& mf_old, LevelData& mf_new, const int lev,
                           const double dt, const double dx) override;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_INTEGRATOR_LSSSPRK3_H_
