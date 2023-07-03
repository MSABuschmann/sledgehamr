#ifndef SLEDGEHAMR_TIMER_H_
#define SLEDGEHAMR_TIMER_H_

#include <chrono>

#include <AMReX_AmrCore.H>

namespace sledgehamr {

class Timer {
  public:
    Timer(std::string my_name) : name{my_name} {};

    void Start();
    void Stop();

    double GetTotalTimeSeconds();
    double GetLastDurationSeconds();

    std::string GetName() const {
        return name;
    };

  private:
    void Check();

    std::chrono::steady_clock::time_point start, stop;
    std::chrono::microseconds last_duration_micro;
    double total_micro = 0;
    const std::string name = "Unknown Timer";
    bool running = false;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_TIMER_H_
