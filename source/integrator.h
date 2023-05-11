#ifndef SLEDGEHAMR_INTEGRATOR_H_
#define SLEDGEHAMR_INTEGRATOR_H_

#include "sledgehamr.h"

namespace sledgehamr {

class Sledgehamr;

/** @brief Abstract base class that handles the time integration for a
 *         single level.
 *         TODO: Currently not actually abstract, needs to be expanded for other
 *               integrators.
 */
class Integrator {
  public:
    Integrator(Sledgehamr* owner);

    /** @brief Advances a single level.
     * @param   lev Level to be advanced.
     */
    virtual void Advance(const int lev);

  protected:
    /** @brief TODO
     */
    virtual void Integrate(LevelData& mf_old, LevelData& mf_new, const int lev,
                           const double dt, const double dx) = 0;

    /** @brief Pointer to owner on whose data this class operates.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_INTEGRATOR_H_
