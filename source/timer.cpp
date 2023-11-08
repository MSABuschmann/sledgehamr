#include "timer.h"

namespace sledgehamr {

void Timer::Start() {
    if (is_running)
        return;

    //amrex::ParallelDescriptor::Barrier();
    start_time = std::chrono::steady_clock::now();
    is_running = true;
}

void Timer::CheckClock() {
    stop_time = std::chrono::steady_clock::now();
    last_duration_micro =
            std::chrono::duration_cast<std::chrono::microseconds>(
                stop_time - start_time);
}

void Timer::Stop() {
    if (!is_running)
        return;

    //amrex::ParallelDescriptor::Barrier();
    CheckClock();
    total_micro +=  static_cast<double>(last_duration_micro.count());
    is_running = false;
}

double Timer::GetTotalTimeSeconds() {
    double extra_micro = 0;
    if (is_running) {
        CheckClock();
        extra_micro = static_cast<double>(last_duration_micro.count());
    }

    return (total_micro + extra_micro)/1e6;
};

double Timer::GetLastDurationSeconds() {
    if (is_running)
        CheckClock();

    return static_cast<double>(last_duration_micro.count())/1e6;
};

}; // namespace sledgehamr
