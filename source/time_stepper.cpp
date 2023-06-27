#include "time_stepper.h"
#include "sledgehamr_utils.h"

#include "integrators/amrex_integrators.h"
#include "integrators/lsssprk3.h"
#include "integrators/rkn.h"

namespace sledgehamr {

TimeStepper::TimeStepper(Sledgehamr* owner) {
    sim = owner;
    local_regrid = new LocalRegrid(sim);

    // Initialize the correct integrator.
    amrex::ParmParse pp_inte("integrator");
    int inte_type;
    pp_inte.get("type", inte_type);
    IntegratorType integrator_type = static_cast<IntegratorType>(inte_type);

    amrex::Print() << "Integrator type: "
                   << Integrator::Name(integrator_type)
                   << std::endl;

    switch (inte_type) {
        case AmrexRkButcherTableau:
            //[[fallthrough]];
        case AmrexForwardEuler:
            //[[fallthough]];
        case AmrexTrapezoid:
            //[[fallthrough]];
        case AmrexSsprk3:
            //[[fallthrough]];
        case AmrexRk4:
            integrator = new IntegratorAMReX(sim);
            break;
        case Lsssprk3:
            integrator = new IntegratorLsssprk3(sim);
            break;
        case RknButcherTableau:
            //[[fallthrough]]; 
        case Rkn4:
            //[[fallthrough]];
        case Rkn5:
            integrator = new IntegratorRkn(sim, integrator_type);
            break;
        default:
            amrex::Abort("#error: Unknown integration type: "
                        + std::to_string(inte_type));
            break;
    }

    // Set regridding intervals.
    amrex::ParmParse pp_amr("amr");

    double reg_dt = 1e99;
    pp_amr.query("regrid_dt", reg_dt);

    for (int lev=0; lev<=sim->max_level; ++lev) {
        regrid_dt.push_back( reg_dt / pow(2, lev) );
        last_regrid_time.push_back( sim->t_start );
    }

    scheduled_regrids.resize( sim->max_level+1 );
}

TimeStepper::~TimeStepper() {
    delete integrator;
    delete local_regrid;
}

void TimeStepper::Advance(int lev) {
    // Perform or schedule regrids.
    if (sim->shadow_hierarchy) {
        // Schedule regrids ahead of time if needed such that we can compute
        // truncation errors in time. In this case regrids will be performed at
        // the end of two time steps.
        ScheduleRegrid(lev);
    } else {
        // Invoke regridding routine at the beginning of a time step if we do
        // not use a shadow hierarchy.
        NoShadowRegrid(lev);
    }

    PreAdvanceMessage(lev);

    // Advance this level.
    utils::sctp timer = utils::StartTimer();
    integrator->Advance(lev);

    PostAdvanceMessage(lev, utils::DurationSeconds(timer));

    // Advance any finer levels twice.
    if (lev != sim->finest_level) {
        Advance(lev+1);
        Advance(lev+1);
    }

    // Synchronize this level with finer/coarser levels.
    SynchronizeLevels(lev);

    // regrid if needed
    DoRegridIfScheduled(lev);

    // Synchronize times to avoid any floating point precision errors from
    // advancing times on each level separately.
    if (lev == 0)
        SynchronizeTimes();
}

void TimeStepper::SynchronizeLevels(int lev) {
    // Check if regrid has been scheduled so we can decide whether we need to
    // compute truncation error erstimate on top of averaging down. Value of
    // 'index' will be -1 of no regrid has been scheduled.
    int istep = sim->grid_new[lev].istep;
    int index = GetIndexOfScheduledRegrid(scheduled_regrids[lev], istep);

    if (lev < sim->finest_level) {
        if (sim->shadow_hierarchy && index != -1) {
            // A regrid is scheduled so we do not synchronize levels quite yet
            // in order to compute truncation errors first.
        } else {
            // Average lev+1 onto lev.
            sim->level_synchronizer->AverageDownTo(lev);
        }
    }

    if (lev >= 1 - sim->shadow_hierarchy && index != -1) {
        // Compute truncation errors for level lev and average down between lev
        // and lev-1.
        sim->level_synchronizer->ComputeTruncationErrors(lev);
    }
}

void TimeStepper::SynchronizeTimes() {
    for (int lev=1; lev <= sim->finest_level; ++lev)
        sim->grid_new[lev].t = sim->grid_new[0].t;
}

void TimeStepper::PreAdvanceMessage(int lev) {
    std::string level_message = LevelMessage(lev, sim->grid_new[lev].istep);

    long ncells = sim->CountCells(lev);
    double coverage_fraction = (double)ncells / pow(sim->dimN[lev],3)*100;
    int nba = sim->grid_new[lev].boxArray().size();

    amrex::Print() << std::left << std::setw(50) << level_message
                   << "Advancing " << ncells << " cells in " << nba
                   << " boxes ... " << "(" << coverage_fraction
                   << "\% coverage)" << std::endl;
}

void TimeStepper::PostAdvanceMessage(int lev, double duration) {
    std::string level_message = LevelMessage(lev, sim->grid_new[lev].istep-1);

    amrex::Print() << std::left << std::setw(50) << level_message
                   << "Advanced to t=" << sim->grid_new[lev].t << " by "
                   << "dt=" << sim->dt[lev] << " in " << duration << "s."
                   << std::endl;
}

std::string TimeStepper::LevelMessage(int lev, int istep) {
    std::string level_name = sledgehamr::utils::LevelName(lev);
    std::string out = "  ";
    for (int i=1;i<=lev;++i)
        out += "| ";
    out += "Level " + std::to_string(lev) + " (" + level_name + ") step #"
          + std::to_string(istep);
    return out;
}

void TimeStepper::ScheduleRegrid(int lev) {
    double time  = sim->grid_new[lev].t;
    int    istep = sim->grid_new[lev].istep;

    // Check if a future regrid has already been scheduled by a coarser level.
    int index = GetIndexOfScheduledRegrid(scheduled_regrids[lev], istep + 1);
    if (index != -1) return;

    // Regrid changes level "lev+1" so we don't regrid on max_level.
    if (lev >= sim->max_level) return;

    // Do not regrid at the end of even time steps as we cannot compute
    // truncation errors otherwise.
    if (istep%2 == 0) return;

    // Check user requirement if we want to invoke a new level. Pass it the
    // level to be created and the time by which the next regrid could be
    // performed if we were to skip this regrid.
    if (!sim->CreateLevelIf(lev+1, time + 3.*sim->dt[lev])) return;

    // Check if enough time since last regrid has passed. We add 3*dt[lev] since
    // we do not want to violate this criteria next time around in case we skip
    // this regrid. Only relevant if local regrid module has not requested an
    // early global regrid.
    if (time + 3.*sim->dt[lev] <= last_regrid_time[lev] + regrid_dt[lev] &&
        !local_regrid->do_global_regrid[lev]) return;

    // Passed all criteria, now schedule regrid. Make sure we schedule the
    // computation of truncation errors at this and all finer levels to be
    // performed at the end of this levels time step and therefore several
    // timesteps away from now for finer levels.
    for (int k = lev; k <= sim->finest_level; ++k) {
        scheduled_regrids[k].push_back(sim->grid_new[k].istep + pow(2,k-lev));
    }

    // Print message.
    std::string level_message = LevelMessage(lev, istep);
    amrex::Print() << std::left << std::setw(50) << level_message
                   << "Regrid scheduled for after time step #"
                   << scheduled_regrids[lev].back()-1 << std::endl;

    if (lev == 0) {
        std::string level_message = LevelMessage(-1, 0);
        amrex::Print() << std::left << std::setw(50) << level_message
                       << "Advancing shadow level." << std::endl;

        sim->CreateShadowLevel();
    } else {
        // Tell the coarser level as well we need truncation error estimates.
        scheduled_regrids[lev-1].push_back( sim->grid_new[lev-1].istep );
    }

    // Mark this as the coarsest level to be regridded.
    regrid_level.push_back(lev);
}

void TimeStepper::DoRegridIfScheduled(int lev) {
    double time  = sim->grid_new[lev].t;
    int    istep = sim->grid_new[lev].istep;

    // No regrid has been scheduled.
    int index = GetIndexOfScheduledRegrid(scheduled_regrids[lev], istep);
    if (index == -1) return;

    // Check if this is really the coarsest level scheduled for a regrid this
    // time.
    if (regrid_level[index] != lev) return;

    // Passed all criteria, remove regrid from schedule.
    for (int k = lev; k <= sim->finest_level; ++k)
        scheduled_regrids[k].erase(scheduled_regrids[k].begin() + index);

    if (lev > 0) {
        scheduled_regrids[lev-1].erase(scheduled_regrids[lev-1].begin()+index);
    }
    regrid_level.erase(regrid_level.begin() + index);

    // Actually do the regrid if we made it this far.
    DoRegrid(lev, time);
}

int TimeStepper::GetIndexOfScheduledRegrid (std::vector<int>& vec, int target) {
    int index = -1;

    if (!vec.empty())
        index = std::distance(vec.begin(), std::find(vec.begin(),
                                                     vec.end(), target));

    if (index == vec.size()) index = -1;

    return index;
}

void TimeStepper::NoShadowRegrid(int lev) {
    double time = sim->grid_new[lev].t;

    // TODO: Implement semi-static case.

    // Skip regrid right at the beginning of the sim. Allowed to be overwritten
    // if no truncation errors are used (TODO).
    if (time == sim->t_start) return;

    // Regrid changes level "lev+1" so we don't regrid on max_level.
    if (lev >= sim->max_level) return;

    // Check if enough time since last regrid has passed. We add dt[lev] since
    // we do not want to violate this criteria next time around in case we skip
    // this regrid. Only relevant if local regrid module has not requested an
    // early global regrid.
    if (time + sim->dt[lev] <= last_regrid_time[lev] + regrid_dt[lev] &&
        !local_regrid->do_global_regrid[lev]) return;

    // Check user requirement if we want to invoke a new level. Pass it the
    // level to be created and the time by which the next regrid could be
    // performed if we were to skip this regrid.
    if (!sim->CreateLevelIf(lev+1, time + sim->dt[lev])) return;

    // Actually do regrid if we made it this far.
    DoRegrid(lev, time);
}

void TimeStepper::DoRegrid(int lev, double time) {
    // Try local regrid first.
    utils::sctp timer = utils::StartTimer();
    bool successfull = local_regrid->AttemptRegrid(lev);

    amrex::Print() << "Local regrid took " << utils::DurationSeconds(timer)
                   << "s." << std::endl;

    // Do global regrid if local regrid failed.
    if (!successfull) {
        amrex::Print() << std::endl << "Perform global regrid at level "
                       << lev+1 << " and higher." << std::endl;

        timer = utils::StartTimer();
        sim->regrid(lev, time);
        local_regrid->DidGlobalRegrid(lev);

        amrex::Print() << "Global regrid took "
                       << utils::DurationSeconds(timer) << "s." << std::endl;
    }

    // Update las regrid times for all levels that have been regridded.
    for (int l=lev; l <= sim->finest_level; ++l) {
        last_regrid_time[l] = time;
    }

    amrex::Print() << std::endl;
}

}; // namespace sledgehamr
