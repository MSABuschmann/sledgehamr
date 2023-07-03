#include "performance_monitor.h"

namespace sledgehamr {

PerformanceMonitor::PerformanceMonitor(Sledgehamr* owner)
  : sim{owner} {
    amrex::ParmParse pp("output");
    pp.query("interval_performance_monitor", interval);
    active = (interval > 0);

    if (active) {
        timer_total.Start();
    }
}

void PerformanceMonitor::Log(hid_t file_id) {

}

}; // namespace sledgehamr
