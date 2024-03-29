#include "performance_monitor.h"
#include "sledgehamr_utils.h"

namespace sledgehamr {

/** @brief Initalizes timers for all tracked functions.
 * @param   owner   Pointer to simulation.
 */
PerformanceMonitor::PerformanceMonitor(Sledgehamr* owner)
  : sim{owner} {
    amrex::ParmParse pp("output.performance_monitor");
    pp.query("interval", interval);
    active = (interval > 0);

    if (!active)
        return;

    idx_total = timer.size();
    timer.emplace_back("Total time");
    timer[idx_total].Start();

    // + 1 for shadow level support.
    idx_rhs = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.emplace_back("::Rhs " + post);
    }

    idx_fill_patch = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.emplace_back("LevelSynchronizer::FillPatch " + post);
    }

    idx_fill_intermediate_patch = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.emplace_back("LevelSynchronizer::FillIntermediatePatch " + post);
    }

    idx_average_down = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.emplace_back("LevelSynchronizer::AverageDownTo " + post);
    }

    idx_truncation_error = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.emplace_back(
                "LevelSynchronizer::ComputeTruncationErrors " + post);
    }

    idx_tagging = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.emplace_back("Sledgehamr::ErrorEst " + post);
    }

    idx_local_regrid = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.emplace_back("LocalRegrid::AttemptRegrid " + post
                           + " (and higher)");
    }

    idx_global_regrid = timer.size() + 1;
    for(int lev = -1; lev <= sim->max_level; ++lev) {
        std::string post = utils::LevelName(lev);
        timer.emplace_back("AmrCore::regrid " + post + " (and higher)");
    }

    idx_read_input = timer.size();
    if (sim->restart_sim)
        timer.emplace_back("IOModule::RestartSim");
    else
        timer.emplace_back("Sledgehamr::InitFromScratch");

    idx_output = timer.size();
    for(OutputModule& out : sim->io_module->output) {
        timer.emplace_back("OutputModule::Write " + out.GetName());
    }
}

/** @brief Starts a timer.
 * @param   id      ID of timer to be started.
 * @param   offset  Offset to be added to the timer ID.
 */
void PerformanceMonitor::Start(int id, int offset) {
    if (active) 
        timer[id + offset].Start();
}

/** @brief Stops a timer.
 * @param   id      ID of timer to be stopped.
 * @param   offset  Offset to be added to the timer ID.
 * @return Time passed since start of timer in seconds. If performance monitor
 *         has not been active then -DBL_MAX will be returned.
 */
double PerformanceMonitor::Stop(int id, int offset) {
    if (active) {
        timer[id + offset].Stop();
        return timer[id + offset].GetLastDurationSeconds();
    } else {
        return -DBL_MAX;
    }
}

/** @brief Sorts all timers by total time passed.
 * @param   timers  Vector of timers.
 * @return Vector of indices that would sort the timers.
 */
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

/** @brief Prints total time passed of all timers.
 * @param  file_id HDF5 file to log the times. Currently not implemented.
 */
void PerformanceMonitor::Log(hid_t file_id) {
    std::vector<int> idx = TimerArgsort(timer);

    amrex::Print() << " ------------------------ PERFORMANCE"
                   << " ------------------------------------\n";
    for(int i : idx) {
        double d = timer[i].GetTotalTimeSeconds();
        if (d != 0) {
            amrex::Print() << std::left << std::setw(60) << timer[i].GetName()
                           << d << "s\n";
        }
    }
    amrex::Print() << " ------------------------------------"
                   << "-------------------------------------" << std::endl;
}

}; // namespace sledgehamr
