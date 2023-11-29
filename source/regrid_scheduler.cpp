#include "regrid_scheduler.h"

namespace sledgehamr {

/** @brief Schedules a regrid.
 * @param   lev Lowest level to be regridded.
 * @param   t   Regrid time.
 */
void RegridScheduler::Schedule(int lev, double t) {
    int id = FindSchedule(t);

    if (id < 0) {
        schedule.emplace_back(lev, t);
    } else {
        schedule[id].lowest_level = std::min(lev, schedule[id].lowest_level);
    }
}

/** @brief Whether we have scheduled a regrid for this level and time.
 * @param   lev Level.
 * @param   t   Time.
 * @return Whether it has been scheduled or not.
 */
bool RegridScheduler::DoRegrid(int lev, double t) const {
    int id = FindSchedule(t);

    if (id < 0)
        return false;

    return (lev == schedule[id].lowest_level);
}

/** @brief Whether we have scheduled a regrid for this level and time and need
 *         to compute truncation errors.
 * @param   lev Level.
 * @param   t   Time.
 * @return Whether it has been scheduled or not.
 */
bool RegridScheduler::NeedTruncationError(int lev, double t) const {
    int id = FindSchedule(t);
        
    if (id < 0)
        return false;

    return (lev >= schedule[id].lowest_level);
}

/** @brief Removes all scheduled regrids at the given time.
 * @param   t   Time.
 */
void RegridScheduler::DidRegrid(double t) {
    int id = FindSchedule(t);

    if (id >= 0) {
        schedule.erase(schedule.begin() + id);
    }
}

/** @brief Finds the internal ID of a regrid scheduled at the given time.
 * @param   t   Time.
 * @return Internal ID of scheduled regrid. Returns -1 if none is found.
 */
int RegridScheduler::FindSchedule(double t) const {
    int id = -1;

    if (!schedule.empty())
        id = std::distance(schedule.begin(), std::find(schedule.begin(),
                           schedule.end(), t));

    if (id == schedule.size())
        id = -1;

    return id;
}

}; // namespace sledgehamr
