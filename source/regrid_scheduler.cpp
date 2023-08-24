#include "regrid_scheduler.h"

namespace sledgehamr {

void RegridScheduler::Schedule(int lev, double t) {
    int id = FindSchedule(t);

    //amrex::Print() << " -------------------------------------------------" 
    //                << "schedule: " << lev << " @ " << t << " " << id << std::endl;

    if (id < 0) {
        schedule.emplace_back(lev, t);
    } else {
        schedule[id].lowest_level = std::min(lev, schedule[id].lowest_level);
    }
}

bool RegridScheduler::DoRegrid(int lev, double t) const {
    int id = FindSchedule(t);
    if (id < 0)
        return false;

    return (lev == schedule[id].lowest_level);
}

bool RegridScheduler::NeedTruncationError(int lev, double t) const {
    int id = FindSchedule(t);
        
    if (id < 0)
        return false;

    return (lev >= schedule[id].lowest_level);
}

void RegridScheduler::DidRegrid(double t) {
    int id = FindSchedule(t);

    if (id >= 0) {
        schedule.erase(schedule.begin() + id);
    }
}

int RegridScheduler::FindSchedule(double t) const {
    int id = -1;

    if (!schedule.empty())
        id = std::distance(schedule.begin(), std::find(schedule.begin(),
                           schedule.end(), t));

    if (id == schedule.size())
        id = -1;

    //amrex::Print() << " -------------------------------------------------" 
    //                << "find: " << t << " " << id << std::endl;

    return id;
}

}; // namespace sledgehamr
