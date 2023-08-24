#ifndef SLEDGEHAMR_REGRID_SCHEDULER_H_
#define SLEDGEHAMR_REGRID_SCHEDULER_H_

#include "sledgehamr.h"

namespace sledgehamr {

class SledgeHAMR;

class ScheduledRegrid {
  public:
    ScheduledRegrid(int lev, double t_regrid)
        : lowest_level{lev}, t{t_regrid} {};

    bool operator==(const double time) const {
        double teps = t*1e-12;
        return (t > time - teps && t < time + teps);
    }

    int lowest_level;
    double t;
};

class RegridScheduler {
  public:
    RegridScheduler() {};

    void Schedule(int lev, double t);
    bool DoRegrid(int lev, double t) const;
    bool NeedTruncationError(int lev, double t) const;
    void DidRegrid(double t);

  private:
    int FindSchedule(double t) const;

    std::vector<ScheduledRegrid> schedule;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_REGRID_SCHEDULER_H_
