#include "performance_monitor.h"
#include "sledgehamr_utils.h"

namespace sledgehamr {

PerformanceMonitor::PerformanceMonitor(Sledgehamr* owner)
  : sim{owner} {
    amrex::ParmParse pp("output");
    pp.query("interval_performance_monitor", interval);
    active = (interval > 0);

    if (!active)
        return;

    idx_total = timer.size();
    timer.push_back(Timer("Total time"));
    timer[idx_total].Start();

    // + 1 for shadow level support.
    idx_rhs = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.push_back(Timer("::RHS " + post));
    }
}

void PerformanceMonitor::Start(int id, int offset) {
    if (active)
        timer[id + offset].Start();
}

double PerformanceMonitor::Stop(int id, int offset) {
    if (active) {
        timer[id + offset].Stop();
        return timer[id + offset].GetLastDurationSeconds();
    } else {
        return -DBL_MAX;
    }
}

void PerformanceMonitor::Log(hid_t file_id) {
    amrex::Print() << " ------------ PERFORMANCE ------------\n";
    for(Timer& t : timer) {
        double d = t.GetTotalTimeSeconds();
        if (d != 0) {
            amrex::Print() << t.GetName() << "\t\t" << d << " s\n";
        }
    }
    amrex::Print() << "--------------------------------------" << std::endl;
}

}; // namespace sledgehamr
