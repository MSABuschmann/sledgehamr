#include "timer.h"

namespace sledgehamr {

void Timer::Start() {
    if (running)
        return;

    amrex::ParallelDescriptor::Barrier();
    start = std::chrono::steady_clock::now();
    running = true;
}

void Timer::Check() {
    amrex::ParallelDescriptor::Barrier();
    stop = std::chrono::steady_clock::now();
    last_duration_micro = 
            std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
}

void Timer::Stop() {
    if (!running)
        return;

    Check();
    total_micro +=  static_cast<double>(last_duration_micro.count());
    running = false;
}

double Timer::GetTotalTimeSeconds() {
    double extra_micro = 0;
    if (running) {
        Check();
        extra_micro = static_cast<double>(last_duration_micro.count());
    }

    return (total_micro + extra_micro)/1e6;
};

double Timer::GetLastDurationSeconds() {
    if (running)
        Check();

    return static_cast<double>(last_duration_micro.count())/1e6;
};

}; // namespace sledgehamr
