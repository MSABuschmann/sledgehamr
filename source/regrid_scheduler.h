#ifndef SLEDGEHAMR_REGRID_SCHEDULER_H_
#define SLEDGEHAMR_REGRID_SCHEDULER_H_

#include "sledgehamr.h"

namespace sledgehamr {

class SledgeHAMR;

/** @brief Simple class that holds a single scheduled regrid. This is used by
 *         the RegridScheduler.
 */
class ScheduledRegrid {
  public:
    /** @brief Does nothing but fetch meta data.
     * @param   lev         Lowest level to be regridded.
     * @param   t_regrid    Time at which we want to regrid.
     */
    ScheduledRegrid(int lev, double t_regrid)
        : lowest_level{lev}, t{t_regrid} {};

    /** @brief Checks whether this regrid has been scheduled at a given time.
     * @param  time    Time to match.
     * @return Whether times match.
     */
    bool operator==(const double time) const {
        double teps = t*1e-12;
        return (t > time - teps && t < time + teps);
    }

    /** @brief Lowest level to be regridded.
     */
    int lowest_level;

    /** @brief Regrid time.
     */
    double t;
};

/** @brief A class that helps us keep track of when and on what levels we want
 *         to regrid. That way we can adjust our workflow accordingly by e.g.
 *         computing truncation error estimates ahead of time or creating
 *         shadow level.
 */
class RegridScheduler {
  public:
    /** @brief Does nothing.
     */
    RegridScheduler() {};

    void Schedule(int lev, double t);
    bool DoRegrid(int lev, double t) const;
    bool NeedTruncationError(int lev, double t) const;
    void DidRegrid(double t);

  private:
    int FindSchedule(double t) const;

    /** @brief Vector containing all scheduled regrids.
     */
    std::vector<ScheduledRegrid> schedule;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_REGRID_SCHEDULER_H_
