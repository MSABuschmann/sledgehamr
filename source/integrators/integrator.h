#ifndef SLEDGEHAMR_INTEGRATOR_H_
#define SLEDGEHAMR_INTEGRATOR_H_

#include "level_data.h"
#include "sledgehamr.h"

namespace sledgehamr {

/** @brief Enum containing all valid integration types. These are the values
 *         that are being used in the inputs file under 'integrator.type'.
 */
enum IntegratorType {
    AmrexRkButcherTableau = 0,
    AmrexForwardEuler = 1,
    AmrexTrapezoid = 2,
    AmrexSsprk3 = 3,
    AmrexRk4 = 4,
    Lsssprk3 = 10,
    Leapfrog = 11,
    RknButcherTableau = 20,
    Rkn4 = 21,
    Rkn5 = 22
};

class Sledgehamr;

/** @brief Abstract base class that handles the time integration for a
 *         single level.
 */
class Integrator {
  public:
    Integrator(Sledgehamr *owner) : sim(owner){};
    virtual void Advance(const int lev);
    static std::string Name(IntegratorType type);
    static void DebugMessage(amrex::MultiFab &mf, std::string msg);

  protected:
    /** @brief Purely virtual function that advances one level by one time step
     *         using an integration scheme of your choice.
     * @param   mf_old  Current state.
     * @param   mf_new  New state after advancement.
     * @param   lev     Current level.
     * @param   dt      Time step size.
     * @param   dx      Grid spacing.
     */
    virtual void Integrate(LevelData &mf_old, LevelData &mf_new, const int lev,
                           const double dt, const double dx) = 0;

    /** @brief Pointer to the simulation.
     */
    Sledgehamr *sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_INTEGRATOR_H_
