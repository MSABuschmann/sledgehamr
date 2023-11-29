#ifndef SLEDGEHAMR_PERFORMANCE_MONITOR_H_
#define SLEDGEHAMR_PERFORMANCE_MONITOR_H_

#include "sledgehamr.h"
#include "timer.h"

namespace sledgehamr {

class PerformanceMonitor {
  public:
    PerformanceMonitor(Sledgehamr* owner);

    void Start(int id, int offset = 0);
    double Stop(int id, int offset = 0);
    void Log(hid_t file_id);

    /** @brief Returns whether the performance monitor is active.
     */
    bool IsActive() const {
        return active;
    };

    /** @brief IDs of timers.
     */
    int idx_total = -1;
    int idx_rhs = -1;
    int idx_fill_patch = -1;
    int idx_fill_intermediate_patch = -1;
    int idx_average_down = -1;
    int idx_truncation_error = -1;
    int idx_tagging = -1;
    int idx_local_regrid = -1;
    int idx_global_regrid = -1;
    int idx_read_input = -1;
    int idx_output = -1;

    /** @brief Vector of all timers.
     */
    std::vector<Timer> timer;

  private:
    std::vector<int> TimerArgsort(std::vector<Timer> timers);

    /** @brief Interval at which we want to print out times.
     */
    double interval = -1;

    /** @brief Is the performance monitor active.
     */
    bool active = false;

    /** @brief Pointer to the simulation.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_PERFORMANCE_MONITOR_H_
