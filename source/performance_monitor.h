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

    void Start(int id, int offset = 0);
    double Stop(int id, int offset = 0);

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

    std::vector<Timer> timer;

  private:
    std::vector<int> TimerArgsort(std::vector<Timer> timers);

    double interval = 0;
    bool active = false;
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_PERFORMANCE_MONITOR_H_
