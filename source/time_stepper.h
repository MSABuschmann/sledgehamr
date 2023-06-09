#ifndef SLEDGEHAMR_TIMESTEPPER_H_
#define SLEDGEHAMR_TIMESTEPPER_H_

#include "sledgehamr.h"
#include "local_regrid.h"
#include "integrators/integrator.h"

namespace sledgehamr {

class SledgeHAMR;
class Integrator;
class LocalRegrid;

/** @brief Class that takes care of the sub-cycling in time algorithm as well as
 *         the scheduling regrids.
 */
class TimeStepper {
  public:
    TimeStepper(Sledgehamr* owner);
    ~TimeStepper ();

    /** @brief Recursive function and the core of the sub-cycling in time
     *         algorithm. Advances a given level by dt[level].
     * @param   lev Level to be advanced.
     */
    void Advance(int lev);

    /** @brief Vector of regridding intervals at each level.
     */
    std::vector<double> regrid_dt;

    /** @brief Pointer to integration module.
     */
    Integrator* integrator;
 
private:
    /** @brief Synchronizes two levels by averaging down. Computes truncation
     *         errors if regrid is scheduled.
     * @param   lev Synchronizes lev with either lev-1 or lev+1 depending on the
     *              regrid status.
     */
    void SynchronizeLevels(int lev);

    /** @brief Synchronizes the time of all levels. This is needed to avoid
     *         de-synchronization due to floating point precission errors during
     *         t -> t + dt operations.
     */
    void SynchronizeTimes();

    /** @brief Prints message before and after each level has been advanced.
     * @param   lev Level that has/will be advanced.
     */
    void PreAdvanceMessage(int lev);
    void PostAdvanceMessage(int lev, double duration);
    std::string LevelMessage(int lev, int istep);

    /** @brief Schedules a regrid ahead of such that we can compute truncations
     *         errors first.
     * @param   lev Level at which to tag cells.
     */
    void ScheduleRegrid(int lev);

    /** @brief Checks whether a regrid has been scheduled and triggers a regrid
     *         if so.
     * @param   lev Current level.
     */
    void DoRegridIfScheduled(int lev);

    /** @brief Checks if a regrid at a given step and level has been scheduled
     *         and returns its index. Will return -1 if not found.
     * @param   vec     Regrid schedule of a level.
     * @param   target  Time step to check.
     * @return  Index of target within vec. -1 if not found.
     */
    int GetIndexOfScheduledRegrid(std::vector<int>& vec, int target);

    /** @brief Performs a regrid if needed if no truncation tags are to be
     *         performed.
     * @param   lev Level at which to tag cells.
     */
    void NoShadowRegrid(int lev);

    /** @brief Performs the actual regrid, either local or global as
     *         appropriate.
     * @param   lev     Level at which to tag cells.
     * @param   time    Current time.
     */
    void DoRegrid(int lev, double time);

    /** @brief Vector of times at which a given level has been regridded last.
     */
    std::vector<double> last_regrid_time;

    /** @brief Vector for each level containing the step numbers at we want to
     *         perform a regrid. Needed such that we can compute truncation
     *         errors in time.
     */
    std::vector<std::vector<int> > scheduled_regrids;

    /** @brief Vector containing for each scheduled regrid which level is the
     *         coarsest.
     */
    std::vector<int> regrid_level;

    /** @brief Pointer to other modules.
     */
    LocalRegrid* local_regrid;
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_TIMESTEPPER_H_
