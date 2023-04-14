#ifndef SLEDGEHAMR_INTEGRATOR_H_
#define SLEDGEHAMR_INTEGRATOR_H_

#include "sledgehamr.h"

namespace sledgehamr {

class Sledgehamr;

/** @brief Pure virtual base class that handles the time integration for a
 *         single level.
 *         TODO: Currently not pure virtual, needs to be expanded for other
 *               integrators.
 */
class Integrator {
  public:
    Integrator(Sledgehamr* owner);

    /** @brief Advances a single level.
     * @param   lev Level to be advanced.
     */
    virtual void Advance(int lev);

  private:

    /** @brief Pointer to owner on whose data this class operates.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_INTEGRATOR_H_
