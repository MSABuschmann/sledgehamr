#include "timer.h"

namespace sledgehamr {

/** @brief Starts the timer by setting the start timer point.
 */
void Timer::Start() {
    if (is_running)
        return;

    //amrex::ParallelDescriptor::Barrier();
    start_time = std::chrono::steady_clock::now();
    is_running = true;
}

/** @brief Stops the timer and checks time passed since it started running.
 */
void Timer::Stop() {
    if (!is_running)
        return;

    //amrex::ParallelDescriptor::Barrier();
    CheckClock();
    total_micro +=  static_cast<double>(last_duration_micro.count());
    is_running = false;
}

/** @brief Returns total time in seconds the timer was running across multiple
 *         start/stop cycles.
 *  return Time in seconds.
 */
double Timer::GetTotalTimeSeconds() {
    double extra_micro = 0;
    if (is_running) {
        CheckClock();
        extra_micro = static_cast<double>(last_duration_micro.count());
    }

    return (total_micro + extra_micro)/1e6;
};

/** @brief Same as GetTotalTimeSeconds() but only time since the timer was 
 *         last started.
 */
double Timer::GetLastDurationSeconds() {
    if (is_running)
        CheckClock();

    return static_cast<double>(last_duration_micro.count())/1e6;
};

/** @brief Computes time passed between now and the last starting time point
           internally.
 */
void Timer::CheckClock() {
    stop_time = std::chrono::steady_clock::now();
    last_duration_micro =
            std::chrono::duration_cast<std::chrono::microseconds>(
                stop_time - start_time);
}

}; // namespace sledgehamr
