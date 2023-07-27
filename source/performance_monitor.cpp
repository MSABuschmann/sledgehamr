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
        timer.push_back(Timer("::Rhs " + post));
    }

    idx_fill_patch = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.push_back(Timer("LevelSynchronizer::FillPatch " + post));
    }

    idx_fill_intermediate_patch = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.push_back(
                Timer("LevelSynchronizer::FillIntermediatePatch " + post));
    }

    idx_average_down = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.push_back(Timer("LevelSynchronizer::AverageDownTo " + post));
    }

    idx_truncation_error = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.push_back(
                Timer("LevelSynchronizer::ComputeTruncationErrors " + post));
    }

    idx_tagging = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.push_back(Timer("Sledgehamr::ErrorEst " + post));
    }

    idx_local_regrid = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.push_back(Timer("LocalRegrid::AttemptRegrid " + post
                            + " (and higher)"));
    }

    idx_global_regrid = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.push_back(Timer("AmrCore::regrid " + post + " (and higher)"));
    }

    idx_read_input = timer.size();
    if (sim->restart_sim)
        timer.push_back(Timer("IOModule::RestartSim"));
    else
        timer.push_back(Timer("Sledgehamr::InitFromScratch"));

    idx_output = timer.size();
    for(OutputModule& out : sim->io_module->output) {
        timer.push_back(Timer("OutputModule::Write " + out.GetName()));
    }
}

void PerformanceMonitor::Start(int id, int offset) {
//    if (id == idx_fill_patch)
//        amrex::AllPrint() << "perf: " << id << " " << offset << " " << active << " " << timer.size() << std::endl;

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

std::vector<int> PerformanceMonitor::TimerArgsort(std::vector<Timer> timers) {
    std::vector<int> idx(timers.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
        [&timers](int i1, int i2) {
            return   timers[i1].GetTotalTimeSeconds()
                   > timers[i2].GetTotalTimeSeconds();
        });

    return idx;
}

void PerformanceMonitor::Log(hid_t file_id) {
    std::vector<int> idx = TimerArgsort(timer);

    amrex::Print() << " ------------ PERFORMANCE ------------\n";
    for(int i : idx) {
        double d = timer[i].GetTotalTimeSeconds();
        if (d != 0) {
            amrex::Print() << std::left << std::setw(60) << timer[i].GetName()
                           << d << "s\n";
        }
    }
    amrex::Print() << "--------------------------------------" << std::endl;
}

}; // namespace sledgehamr
