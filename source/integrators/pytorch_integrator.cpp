#include <torch/script.h>

#include "pytorch_integrator.h"

namespace sledgehamr {

/** @brief Advances one level by one time step using the low-storage strong
 *         stability preserving third order Runge-Kutte integration scheme
 *         (LSSSPRK3).
 * @param   mf_old  Current state.
 * @param   mf_new  New state after advancement.
 * @param   lev     Current level.
 * @param   dt      Time step size.
 * @param   dx      Grid spacing.
 */
void IntegratorPytorch::Integrate(LevelData &mf_old, LevelData &mf_new,
                                  const int lev, const double dt,
                                  const double dx) {}

}; // namespace sledgehamr
