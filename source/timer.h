#ifndef SLEDGEHAMR_TIMER_H_
#define SLEDGEHAMR_TIMER_H_

#include <chrono>

#include <AMReX_AmrCore.H>

namespace sledgehamr {

/** @brief This class keeps track of passed time just like stopwatch.
 */
class Timer {
  public:
    /** @brief Give the timer a display name during construction.
     *  @param  my_name Display name of timner.
     */
    Timer(std::string my_name) : name{my_name} {};

    void Start();
    void Stop();

    double GetTotalTimeSeconds();
    double GetLastDurationSeconds();

    /** @brief Returns display name of timer.
     */
    std::string GetName() const {
        return name;
    };

    /** @brief Returns whether the timer is currently running.
     */
    bool IsRunning() const {
        return is_running;
    };

  private:
    void CheckClock();

    /** @brief Time point of when pointer has been started last.
     */
    std::chrono::steady_clock::time_point stop_time;

    /** @brief Time point of when pointer has been stopped last.
     */
    std::chrono::steady_clock::time_point stop_time;

    /** @brief Time in milliseconds from the last starting time point to when
     *         Check() has been called last.
     */
    std::chrono::microseconds last_duration_micro;

    /** @brief Total time in milliseconds that the timer has been running
     *         overall (across multiple start/stop cycles).
     */
    double total_micro = 0;

    /** @brief Display name of timer.
     */
    const std::string name = "Unknown Timer";

    /** @brief Flag whether timer is currently running.
     */
    bool is_running = false;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_TIMER_H_
