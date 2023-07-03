#ifndef SLEDGEHAMR_PERFORMANCE_MONITOR_H_
#define SLEDGEHAMR_PERFORMANCE_MONITOR_H_

#include "sledgehamr.h"

namespace sledgehamr {

class PerformanceMonitor {
  public:
    PerformanceMonitor(Sledgehamr* owner);

    void Log(hid_t file_id);

    bool IsActive() const {
        return active;
    };

  private:
    double interval = 0;
    bool active = false;
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_PERFORMANCE_MONITOR_H_
