#ifndef SLEDGEHAMR_PERFORMANCE_MONITOR_H_
#define SLEDGEHAMR_PERFORMANCE_MONITOR_H_

#include "sledgehamr.h"
#include "timer.h"

namespace sledgehamr {

class PerformanceMonitor {
  public:
    PerformanceMonitor(Sledgehamr* owner);

    void Log(hid_t file_id);

    bool IsActive() const {
        return active;
    };

    // Per level.
    std::vector<Timer> timer_rhs;
    std::vector<Timer> timer_fill_patch;
    std::vector<Timer> timer_fill_intermediate_patch;
    std::vector<Timer> timer_average_down;
    std::vector<Timer> timer_tagging;
    std::vector<Timer> timer_local_regrid;
    std::vector<Timer> timer_global_regrid;

    // For each output type.
    std::vector<Timer> timer_output;
    Timer timer_total("total");

  private:
    double interval = 0;
    bool active = false;
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_PERFORMANCE_MONITOR_H_
