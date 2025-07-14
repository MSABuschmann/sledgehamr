#include "time_stepper.h"
#include "sledgehamr_utils.h"

#include "integrators/amrex_integrators.h"
#include "integrators/leapfrog.h"
#include "integrators/lsssprk3.h"
#include "integrators/rkn.h"

namespace sledgehamr {

/** @brief Initalizes all unique instances for regridding and initializes the
 *         integrator.
 * @param   owner   Pointer to the simulation.
 */
TimeStepper::TimeStepper(Sledgehamr *owner) : sim(owner) {
    local_regrid = std::make_unique<LocalRegrid>(sim);
    scheduler = std::make_unique<RegridScheduler>();
    ParseParams();
    SetIntegrator();
}

/** @brief Initalizes the integrator class. Will abort the entire program if
 *         no corresponding integrator could be found.
 */
void TimeStepper::SetIntegrator() {
    // Initialize the correct integrator.
    amrex::ParmParse pp_inte("integrator");
    int inte_type;
    pp_inte.get("type", inte_type);
    IntegratorType integrator_type = static_cast<IntegratorType>(inte_type);

    amrex::Print() << "Integrator type: " << Integrator::Name(integrator_type)
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
        integrator = std::make_unique<IntegratorAMReX>(sim);
        break;
    case Lsssprk3:
        integrator = std::make_unique<IntegratorLsssprk3>(sim);
        break;
    case Leapfrog:
        integrator = std::make_unique<IntegratorLeapfrog>(sim);
        break;
    case RknButcherTableau:
        //[[fallthrough]];
    case Rkn4:
        //[[fallthrough]];
    case Rkn5:
        integrator = std::make_unique<IntegratorRkn>(sim, integrator_type);
        break;
    default:
        amrex::Abort("#error: Unknown integration type: " +
                     std::to_string(inte_type));
        break;
    }
}

/** @brief Parses all parameters needed for the sub-cycling algorithm.
 */
void TimeStepper::ParseParams() {
    // Set regridding intervals.
    amrex::ParmParse pp_amr("amr");

    double reg_dt = DBL_MAX;
    pp_amr.query("regrid_dt", reg_dt);
    pp_amr.query("semistatic_sim", semistatic_sim);

    for (int lev = 0; lev <= sim->max_level; ++lev) {
        regrid_dt.push_back(reg_dt / pow(2, lev));
        last_regrid_time.push_back(sim->t_start);
    }

    amrex::ParmParse pp_out("output");
    pp_out.query("output_of_initial_state", output_of_initial_state);
}

/** @brief Recursive function and the core of the sub-cycling in time
 *         algorithm. Advances a given level by dt[level].
 * @param   lev Level to be advanced.
 */
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

    // For self-consistency invoke function only if we don't have a shadow
    // level.
    if (lev == 0 && !sim->shadow_level.isDefined())
        sim->BeforeTimestep(sim->grid_new[lev].t);

    if (sim->grid_new[0].t == sim->t_start && output_of_initial_state)
        sim->io_module->Write(true);

    // Advance this level.
    PreAdvanceMessage(lev);
    utils::sctp timer = utils::StartTimer();
    integrator->Advance(lev);
    PostAdvanceMessage(lev, utils::DurationSeconds(timer));

    // Advance any finer levels twice.
    if (lev != sim->finest_level) {
        Advance(lev + 1);
        Advance(lev + 1);
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

/** @brief Synchronizes two levels by averaging down. Computes truncation
 *         errors if regrid is scheduled.
 * @param   lev Synchronizes lev with either lev-1 or lev+1 depending on the
 *              regrid status.
 */
void TimeStepper::SynchronizeLevels(int lev) {
    // Check if regrid has been scheduled so we can decide whether we need to
    // compute truncation error erstimate on top of averaging down. Value of
    // 'index' will be -1 of no regrid has been scheduled.
    bool need_truncation_errors =
        scheduler->NeedTruncationError(lev, sim->grid_new[lev].t);

    if (lev < sim->finest_level) {
        if (sim->shadow_hierarchy && need_truncation_errors) {
            // A regrid is scheduled so we do not synchronize levels quite yet
            // in order to compute truncation errors first.
        } else {
            // Average lev+1 onto lev.
            sim->level_synchronizer->AverageDownTo(lev);
        }
    }

    if (lev >= 1 - sim->shadow_hierarchy && need_truncation_errors) {
        // Compute truncation errors for level lev and average down between lev
        // and lev-1.
        sim->level_synchronizer->ComputeTruncationErrors(lev);
    }
}

/** @brief Synchronizes the time of all levels. This is needed to avoid
 *         de-synchronization due to floating point precission errors during
 *         t -> t + dt operations.
 */
void TimeStepper::SynchronizeTimes() {
    for (int lev = 1; lev <= sim->finest_level; ++lev)
        sim->grid_new[lev].t = sim->grid_new[0].t;
}

/** @brief Prints message just before a level has been advanced.
 * @param   lev Level that will be advanced.
 */
void TimeStepper::PreAdvanceMessage(int lev) {
    std::string level_message = LevelMessage(lev, sim->grid_new[lev].istep);

    long ncells = sim->CountCells(lev);
    double coverage_fraction = (double)ncells / pow(sim->dimN[lev], 3) * 100;
    int nba = sim->grid_new[lev].boxArray().size();

    amrex::Print() << std::left << std::setw(50) << level_message
                   << "Advancing " << ncells << " cells in " << nba
                   << " boxes ... " << "(" << coverage_fraction
                   << "\% coverage)" << std::endl;
}

/** @brief Prints message right after a level has been advanced.
 * @param   lev Level that has be advanced.
 */
void TimeStepper::PostAdvanceMessage(int lev, double duration) {
    std::string level_message = LevelMessage(lev, sim->grid_new[lev].istep - 1);

    amrex::Print() << std::left << std::setw(50) << level_message
                   << "Advanced to t=" << sim->grid_new[lev].t << " by "
                   << "dt=" << sim->dt[lev] << " in " << duration << "s."
                   << " (" << amrex::ParallelDescriptor::second()
                   << "s since start)" << std::endl;
}

/** @brief Returns a display message containing a level and its step number.
 * @param   lev     Level to appear in message.
 * @param   istep   Corresponding step number.
 * @return Display message.
 */
std::string TimeStepper::LevelMessage(int lev, int istep) {
    std::string level_name = sledgehamr::utils::LevelName(lev);
    std::string out = "  ";
    for (int i = 1; i <= lev; ++i)
        out += "| ";
    out += "Level " + std::to_string(lev) + " (" + level_name + ") step #" +
           std::to_string(istep);
    return out;
}

/** @brief Schedules a regrid ahead of such that we can compute truncations
 *         errors first if needed.
 * @param   lev Level at which to tag cells.
 */
void TimeStepper::ScheduleRegrid(int lev) {
    double time = sim->grid_new[lev].t;
    int istep = sim->grid_new[lev].istep;

    // Check if a future regrid has already been scheduled by a coarser level.
    if (scheduler->DoRegrid(lev, time + sim->dt[lev])) {
        return;
    }

    // Sanity check we can actually compute a shadow level. Prevents a regrid on
    // level 1 right after a restart.
    if (sim->grid_old[lev].t == time) {
        return;
    }

    // Regrid changes level "lev+1" so we don't regrid on max_level.
    if (lev >= sim->max_level) {
        return;
    }

    // Do not regrid at the end of even time steps as we cannot compute
    // truncation errors otherwise. Not relevant for coarse level as we
    // can create shadow level whenever.
    if (istep % 2 == 0 && lev > 0) {
        return;
    }

    // Again, since we can create shadow levels whenever we could regrid
    // earlier for coarse level.
    double time_next_opportunity =
        lev > 0 ? time + 3. * sim->dt[lev] : time + 2. * sim->dt[lev];

    // Check user requirement if we want to invoke a new level. Pass it the
    // level to be created and the time by which the next regrid could be
    // performed if we were to skip this regrid.
    if (!sim->DoCreateLevelIf(lev + 1, time_next_opportunity)) {
        return;
    }

    // Check if enough time since last regrid has passed. We add 3*dt[lev] since
    // we do not want to violate this criteria next time around in case we skip
    // this regrid. Only relevant if local regrid module has not requested an
    // early global regrid.
    if (time_next_opportunity <= last_regrid_time[lev] + regrid_dt[lev] &&
        !local_regrid->do_global_regrid[lev]) {
        return;
    }

    // This is to avoid regridding on coarse right after restarting from a
    // checkpoint. We do not have a valid grid_old yet to evolve the needed
    // shadow level.
    if (lev == 0 && sim->grid_old[lev].t == -DBL_MAX) {
        return;
    }

    // Passed all criteria, now schedule regrid. Make sure we schedule the
    // computation of truncation errors at this and all finer levels to be
    // performed at the end of this levels time step and therefore several
    // timesteps away from now for finer levels.
    scheduler->Schedule(lev, time + sim->dt[lev]);

    // Print message.
    std::string level_message = LevelMessage(lev, istep);
    amrex::Print() << std::left << std::setw(50) << level_message
                   << "Regrid scheduled for after time step #" << istep << "."
                   << std::endl;

    if (lev == 0) {
        std::string level_message = LevelMessage(-1, 0);
        amrex::Print() << std::left << std::setw(50) << level_message
                       << "Advancing shadow level." << std::endl;

        sim->CreateShadowLevel();
    }
}

/** @brief Checks whether a regrid has been scheduled and triggers a regrid
 *         if so.
 * @param   lev Current level.
 */
void TimeStepper::DoRegridIfScheduled(int lev) {
    double time = sim->grid_new[lev].t;
    int istep = sim->grid_new[lev].istep;

    if (!scheduler->DoRegrid(lev, time))
        return;

    // Actually do the regrid if we made it this far.
    DoRegrid(lev, time);
    scheduler->DidRegrid(time);
}

/** @brief Performs a regrid if needed if no truncation tags are to be
 *         performed.
 * @param   lev Level at which to tag cells.
 */
void TimeStepper::NoShadowRegrid(int lev) {
    double time = sim->grid_new[lev].t;

    // Regrid changes level "lev+1" so we don't regrid on max_level.
    if (lev >= sim->max_level && !semistatic_sim)
        return;

    // Check if enough time since last regrid has passed. We add dt[lev] since
    // we do not want to violate this criteria next time around in case we skip
    // this regrid. Only relevant if local regrid module has not requested an
    // early global regrid.
    if (time + sim->dt[lev] <= last_regrid_time[lev] + regrid_dt[lev] &&
        !local_regrid->do_global_regrid[lev])
        return;

    // Check user requirement if we want to invoke a new level. Pass it the
    // level to be created and the time by which the next regrid could be
    // performed if we were to skip this regrid.
    if (!sim->DoCreateLevelIf(lev + 1, time + sim->dt[lev]))
        return;

    // Actually do regrid if we made it this far.
    DoRegrid(lev, time);
}

/** @brief Performs the actual regrid, either local or global as
 *         appropriate.
 * @param   lev     Level at which to tag cells.
 * @param   time    Current time.
 */
void TimeStepper::DoRegrid(int lev, double time) {
    if (semistatic_sim) {
        sim->level_synchronizer->IncreaseCoarseLevelResolution();
        return;
    }

    // Before regridding we want to attempt to write output to allow for
    // truncation errors estimates to be written to file. Those would be
    // corrupted or non-exist at any other stage of the simulation.
    if (lev == 0)
        sim->io_module->Write();

    // Try local regrid first.
    utils::sctp timer = utils::StartTimer();
    sim->performance_monitor->Start(sim->performance_monitor->idx_local_regrid,
                                    lev);
    bool successfull = local_regrid->AttemptRegrid(lev);

    sim->performance_monitor->Stop(sim->performance_monitor->idx_local_regrid,
                                   lev);
    amrex::Print() << "Local regrid took " << utils::DurationSeconds(timer)
                   << "s." << std::endl;

    // Do global regrid if local regrid failed.
    if (!successfull) {
        amrex::Print() << std::endl
                       << "Perform global regrid at level " << lev + 1
                       << " and higher." << std::endl;

        timer = utils::StartTimer();
        sim->performance_monitor->Start(
            sim->performance_monitor->idx_global_regrid, lev);

        sim->regrid(lev, time);

        sim->performance_monitor->Stop(
            sim->performance_monitor->idx_global_regrid, lev);

        local_regrid->DidGlobalRegrid(lev);

        amrex::Print() << "Global regrid took " << utils::DurationSeconds(timer)
                       << "s." << std::endl;
    }

    // Update las regrid times for all levels that have been regridded.
    for (int l = lev; l <= sim->finest_level; ++l) {
        last_regrid_time[l] = time;
    }

    amrex::Print() << std::endl;
}

}; // namespace sledgehamr
