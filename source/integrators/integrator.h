#ifndef SLEDGEHAMR_INTEGRATOR_H_
#define SLEDGEHAMR_INTEGRATOR_H_

#include "sledgehamr.h"
#include "level_data.h"

namespace sledgehamr {

/** @brief
 */
enum IntegratorType {
    AmrexRkButcherTableau = 0,
    AmrexForwardEuler = 1,
    AmrexTrapezoid = 2,
    AmrexSsprk3 = 3,
    AmrexRk4 = 4,
    Lsssprk3 = 10,
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
    Integrator(Sledgehamr* owner);

    /** @brief Advances a single level.
     * @param   lev Level to be advanced.
     */
    virtual void Advance(const int lev);

    /** @brief
     */
    static std::string Name(IntegratorType type);
    static void DebugMessage(amrex::MultiFab& mf, std::string msg);

  protected:
    /** @brief
     */
    virtual void Integrate(LevelData& mf_old, LevelData& mf_new, const int lev,
                           const double dt, const double dx) = 0;

    /** @brief Pointer to owner on whose data this class operates.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_INTEGRATOR_H_
