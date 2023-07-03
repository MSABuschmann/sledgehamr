#include "performance_monitor.h"

namespace sledgehamr {

PerformanceMonitor::PerformanceMonitor(Sledgehamr* owner)
  : sim{owner} {
    amrex::ParmParse pp("output");
    pp.query("interval_performance_monitor", interval);
    active = (interval > 0);
}

}; // namespace sledgehamr
