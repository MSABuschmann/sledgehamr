#ifndef SLEDGEHAMR_TIMESTEPPER_H_
#define SLEDGEHAMR_TIMESTEPPER_H_

#include "sledgehamr.h"
#include "local_regrid/local_regrid.h"
#include "regrid_scheduler.h"
#include "integrators/integrator.h"

namespace sledgehamr {

class SledgeHAMR;
class Integrator;
class LocalRegrid;
class RegridScheduler;

/** @brief Class that takes care of the sub-cycling in time algorithm and
 *         schedules regrids when necessary.
 */
class TimeStepper {
  public:
    TimeStepper(Sledgehamr* owner);
    void Advance(int lev);

    /** @brief Vector of regridding intervals at each level.
     */
    std::vector<double> regrid_dt;

    /** @brief Pointer to integration module.
     */
    std::unique_ptr<Integrator> integrator;

    /** @brief Pointer to the local regrid module.
     */
    std::unique_ptr<LocalRegrid> local_regrid;

    /** @brief Pointer to the regrid scheduler.
     */
    std::unique_ptr<RegridScheduler> scheduler;

    /** @brief Vector of times at which a given level has been regridded last.
     */
    std::vector<double> last_regrid_time;

private:
    void SynchronizeLevels(int lev);
    void SynchronizeTimes();

    void PreAdvanceMessage(int lev);
    void PostAdvanceMessage(int lev, double duration);
    std::string LevelMessage(int lev, int istep);

    void ScheduleRegrid(int lev);
    void DoRegridIfScheduled(int lev);
    void NoShadowRegrid(int lev);
    void DoRegrid(int lev, double time);

    void ParseParams();
    void SetIntegrator();

    /** @brief Whether we want to force write output at the beginning of the
     *         sim.
     */
    bool output_of_initial_state = true;

    /** @brief Whether we are running a semi-static sim.
     */
    bool semistatic_sim = false;

    /** @brief Pointer to the simulation.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_TIMESTEPPER_H_
