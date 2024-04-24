#include <type_traits>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include "fill_level.h"
#include "hdf5_utils.h"
#include "sledgehamr.h"
#include "sledgehamr_utils.h"

namespace sledgehamr {

/** @brief Creates instances of submodules and reads input parameters.
 */
Sledgehamr::Sledgehamr() {
    amrex::Print() << "\nStarting sledgehamr..." << std::endl;

    ParseInput();

    // Initialize modules.
    time_stepper = std::make_unique<TimeStepper>(this);
    io_module = std::make_unique<IOModule>(this);

    grid_new.resize(max_level + 1);
    grid_old.resize(max_level + 1);

    for (int lev = 0; lev <= max_level; ++lev) {
        dimN.push_back(coarse_level_grid_size * pow(2, lev));
        dx.push_back(L / (double)dimN[lev]);
        dt.push_back(dx[lev] * cfl);
    }

    DoPrerunChecks();
}

/** @brief Initalizes data from scratch or from checkpoint file. Sets up the
 *         remaining modules and updates them.
 */
void Sledgehamr::InitSledgehamr() {
    if (with_gravitational_waves)
        gravitational_waves = std::make_unique<GravitationalWaves>(this);

    // Initialize here and not in the SledgeHAMR constructor such that it knows
    // about the number of scalar fields during construction. Necessary so it
    // can initialize boundary conditions.
    level_synchronizer = std::make_unique<LevelSynchronizer>(this);
    performance_monitor = std::make_unique<PerformanceMonitor>(this);

    ParseInputScalars();

    if (nerrors > 0) {
        amrex::ParallelDescriptor::Barrier();
        amrex::Abort("Found " + std::to_string(nerrors) + " error(s)");
    }

    if (no_simulation) {
        return;
    }

    performance_monitor->Start(performance_monitor->idx_read_input);

    if (restart_sim) {
        io_module->RestartSim();
    } else {
        InitFromScratch(t_start);
    }

    if (increase_coarse_level_resolution)
        level_synchronizer->IncreaseCoarseLevelResolution();

    performance_monitor->Stop(performance_monitor->idx_read_input);

    // Initialize project
    Init();

    io_module->UpdateOutputModules();
}

/** @brief Starts the evolution
 */
void Sledgehamr::Evolve() {
    if (no_simulation)
        return;

    amrex::Print() << "Starting evolution!" << std::endl;

    // Main loop over time.
    while (!StopRunning(grid_new[0].t)) {
        // Advance all levels starting at lev=0. This performs an entire
        // shadow/coarse level time step.
        amrex::Print() << std::endl;
        utils::sctp timer = utils::StartTimer();

        time_stepper->Advance(0);

        // Write any output if requested.
        amrex::Print() << "Full step took " << utils::DurationSeconds(timer)
                       << "s.\n"
                       << std::endl;
        io_module->Write();

#ifdef AMREX_MEM_PROFILING
        std::ostringstream ss;
        ss << "[COARSE STEP " << grid_new[0].istep << "]";
        amrex::MemProfiler::report(ss.str());
#endif
    }

    // Force write at the end of simulation.
    io_module->Write(true);

    amrex::Print() << "Finished!" << std::endl;
}

/** @brief Make a new level from scratch using provided BoxArray and
 *         DistributionMapping. Only used during initialization. Overrides
 *         the pure virtual function in amrex::AmrCore.
 * @param   lev     Level to be created.
 * @param   time    Time of new grid.
 * @param   ba      New amrex::BoxArray.
 * @param   dm      New amrex::DistributionMapping.
 */
void Sledgehamr::MakeNewLevelFromScratch(int lev, amrex::Real time,
                                         const amrex::BoxArray &ba,
                                         const amrex::DistributionMapping &dm) {
    const int ncomp = scalar_fields.size();

    // Define lowest level from scratch.
    grid_new[lev].define(ba, dm, ncomp, nghost, time);
    grid_old[lev].define(ba, dm, ncomp, nghost);

    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // Fill current level lev with initial state data.
    FillLevel fill_level(this, lev);
    fill_level.FromInitialStateFile();
}

/** @brief Make a new level using provided BoxArray and DistributionMapping,
 *         and fills it with interpolated coarse level data. Overrides the
 *         pure virtual function in amrex::AmrCore.
 * @param   lev     Level to be created.
 * @param   time    Time of new grid.
 * @param   ba      New amrex::BoxArray.
 * @param   dm      New amrex::DistributionMapping.
 */
void Sledgehamr::MakeNewLevelFromCoarse(int lev, amrex::Real time,
                                        const amrex::BoxArray &ba,
                                        const amrex::DistributionMapping &dm) {
    const int ncomp = grid_new[lev - 1].nComp();
    const int nghost = grid_new[lev - 1].nGrow();

    grid_new[lev].define(ba, dm, ncomp, nghost, time);
    grid_old[lev].define(ba, dm, ncomp, nghost);

    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // Fill new level with coarse level data.
    level_synchronizer->FillCoarsePatch(lev, time, grid_new[lev]);
}

/** @brief Remake a new level using provided BoxArray and
 *         DistributionMapping, and fills it with interpolated coarse level
 *         data. Overrides the pure virtual function in amrex::AmrCore.
 * @param   lev     Level to be remade.
 * @param   time    Time of new grid.
 * @param   ba      New amrex::BoxArray.
 * @param   dm      New amrex::DistributionMapping.
 */
void Sledgehamr::RemakeLevel(int lev, amrex::Real time,
                             const amrex::BoxArray &ba,
                             const amrex::DistributionMapping &dm) {
    const int ncomp = grid_new[lev].nComp();
    const int nghost = grid_new[lev].nGrow();

    // Remake new_grid and fill with data.
    LevelData new_state(ba, dm, ncomp, nghost, grid_new[lev].t);
    new_state.istep = grid_new[lev].istep;
    level_synchronizer->FillPatch(lev, time, new_state);
    std::swap(new_state, grid_new[lev]);
    new_state.clear();

    // Remake old_grid.
    grid_old[lev].clear();
    grid_old[lev].define(ba, dm, ncomp, nghost);
}

/** @brief Delete level data. Overrides the pure virtual function in
 *         amrex::AmrCore.
 * @param   lev Level to be deleted.
 */
void Sledgehamr::ClearLevel(int lev) {
    grid_new[lev].clear();
    grid_old[lev].clear();
}

/** @brief Tag cells for refinement. Overrides the pure virtual function in
 *         amrex::AmrCore.
 * @param   lev         Level on which cells are tagged.
 * @param   time        Time of said level.
 * @param   ngrow       Grid growth factor.
 * @param   ntags_user  Counts number of user-defined tags.
 */
void Sledgehamr::ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time,
                          int ngrow) {
    // Skip regrid right at the beginning of the sim. Allowed to be overridden
    // if no truncation errors are used.
    if (time == t_start && shadow_hierarchy)
        return;

    // Also skip if level is not supposed to be created yet.
    if (!DoCreateLevelIf(lev + 1, time))
        return;

    performance_monitor->Start(performance_monitor->idx_tagging, lev);
    utils::sctp timer = utils::StartTimer();

    if (tagging_on_gpu)
        DoErrorEstGpu(lev, tags, time);
    else
        DoErrorEstCpu(lev, tags, time);

    amrex::Print() << "  Tagging took " << utils::DurationSeconds(timer) << "s."
                   << std::endl;
    performance_monitor->Stop(performance_monitor->idx_tagging, lev);
}

/** @brief Creates a shadow level and evolves it by one time step. Needed
 *         to compute truncation errors on the coarse level.
 */
void Sledgehamr::CreateShadowLevel() {
    const int ncomp = scalar_fields.size();
    const double time = grid_old[0].t;
    amrex::BoxArray ba = grid_old[0].boxArray();
    ba.coarsen(2);
    amrex::DistributionMapping dm = grid_old[0].DistributionMap();

    shadow_level.define(ba, dm, ncomp, nghost);
    shadow_level_tmp.define(ba, dm, ncomp, nghost, time);

    shadow_level_geom = geom[0];
    shadow_level_geom.coarsen(amrex::IntVect(2, 2, 2));

    amrex::average_down(grid_old[0], shadow_level_tmp, geom[0],
                        shadow_level_geom, 0, ncomp, refRatio(0));

    time_stepper->integrator->Advance(-1);
}

/** @brief Checks whether we want to save the coarse level box layout for
 *         chunking of the initial state.
 */
void Sledgehamr::DoPrerunChecks() {
    if (get_box_layout_nodes > 0)
        DetermineBoxLayout();
}

/** @brief Saves the coarse level box layout.
 */
void Sledgehamr::DetermineBoxLayout() {
    amrex::Print() << "Get box layout for " << get_box_layout_nodes
                   << " nodes and exit ..." << std::endl;

    amrex::Box bx(amrex::IntVect(0),
                  amrex::IntVect(coarse_level_grid_size - 1));
    amrex::BoxArray ba(bx);
    ChopGrids(0, ba, get_box_layout_nodes);
    io_module->WriteBoxArray(ba);
    no_simulation = true;
}

/** @brief Performs error estimation (i.e. tagging) on CPUs.
 * @param  lev     Tagging level.
 * @param  tags    Container to save tags.
 * @param  time    Current time.
 */
void Sledgehamr::DoErrorEstCpu(int lev, amrex::TagBoxArray &tags, double time) {
    // Current state.
    const amrex::MultiFab &state = grid_new[lev];

    // State containing truncation errors if they have been calculated.
    const LevelData &state_te = grid_old[lev];

    // Initialize tag counters.
    long ntags_total = 0;
    long ntags_user = 0;
    std::vector<long> ntags_trunc(scalar_fields.size(), 0);

    std::vector<double> params_tag, params_mod;
    SetParamsTagCellForRefinement(params_tag, time, lev);
    if (shadow_hierarchy)
        SetParamsTruncationModifier(params_mod, time, lev);

        // Loop over boxes and cells.
#pragma omp parallel reduction(+ : ntags_total) reduction(+ : ntags_user)      \
    reduction(vec_long_plus : ntags_trunc)
    for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
        const amrex::Box &tilebox = mfi.tilebox();
        const amrex::Array4<double const> &state_fab = state.array(mfi);
        const amrex::Array4<double const> &state_fab_te = state_te.array(mfi);
        const amrex::Array4<char> &tag_arr = tags.array(mfi);

        // Tag with or without truncation errors.
        if (shadow_hierarchy && state_te.contains_truncation_errors) {
            TagWithTruncationCpu(state_fab, state_fab_te, tag_arr, tilebox,
                                 time, lev, &ntags_total, &ntags_user,
                                 &(ntags_trunc[0]), params_tag, params_mod);
        } else if (shadow_hierarchy) {
            std::string msg = "Trying to tag using truncation errors but no ";
            msg += "truncation errors are computed on level ";
            msg += std::to_string(lev) + "!";
            amrex::Abort(msg.c_str());
        } else {
            TagWithoutTruncationCpu(state_fab, tag_arr, tilebox, time, lev,
                                    &ntags_total, params_tag);
        }
    }

    // Collect all tags across MPI ranks.
    amrex::ParallelDescriptor::ReduceLongSum(ntags_total, 0);

    if (shadow_hierarchy) {
        amrex::ParallelDescriptor::ReduceLongSum(ntags_user, 0);
        amrex::ParallelDescriptor::ReduceLongSum(&(ntags_trunc[0]),
                                                 ntags_trunc.size(), 0);
    }

    // Print statistics.
    long ncells = CountCells(lev);
    double ftotal = (double)ntags_total / (double)ncells;
    double fuser = (double)ntags_user / (double)ncells;
    amrex::Print() << "  Tagged cells at level " << lev << ": " << ntags_total
                   << " of " << ncells << " (" << ftotal * 100. << "\%)"
                   << std::endl;

    if (shadow_hierarchy) {
        amrex::Print() << "    User-defined tags: " << ntags_user << std::endl;

        for (int i = 0; i < scalar_fields.size(); ++i) {
            amrex::Print() << "    Truncation error tags on "
                           << scalar_fields[i]->name << ": " << ntags_trunc[i]
                           << std::endl;
        }
    }
}

/** @brief Same as Sledgehamr::DoErrorEstCpu but on GPUs. Will not keep track
 *         of tagging statistics.
 */
void Sledgehamr::DoErrorEstGpu(int lev, amrex::TagBoxArray &tags, double time) {
    // Current state.
    const amrex::MultiFab &state = grid_new[lev];

    // State containing truncation errors if they have been calculated.
    const amrex::MultiFab &state_te = grid_old[lev];

    std::vector<double> params_tag, params_mod;
    SetParamsTagCellForRefinement(params_tag, time, lev);
    if (shadow_hierarchy)
        SetParamsTruncationModifier(params_mod, time, lev);

        // Loop over boxes and cells.
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(state, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box &tilebox = mfi.tilebox();
        const amrex::Array4<double const> &state_fab = state.array(mfi);
        const amrex::Array4<double const> &state_fab_te = state_te.array(mfi);
        const amrex::Array4<char> &tag_arr = tags.array(mfi);

        // Tag with or without truncation errors.
        if (shadow_hierarchy) {
            TagWithTruncationGpu(state_fab, state_fab_te, tag_arr, tilebox,
                                 time, lev, params_tag, params_mod);
        } else {
            TagWithoutTruncationGpu(state_fab, tag_arr, tilebox, time, lev,
                                    params_tag);
        }
    }

    amrex::Print() << "  Tagged cells at level " << lev << "." << std::endl;
}

/** @brief Parses and tests all relevant input parameters.
 */
void Sledgehamr::ParseInput() {
    amrex::ParmParse pp("");

    pp.query("input.do_parameter_check_for_n_mpi_ranks", check_mpi_ranks);
    do_thorough_checks = check_mpi_ranks > 0;
    if (!do_thorough_checks)
        check_mpi_ranks = amrex::ParallelDescriptor::NProcs();
    else
        no_simulation = true;

    utils::ErrorState validity =
        (utils::ErrorState)utils::IsPowerOfTwo(check_mpi_ranks);
    std::string param_name = "#MPI ranks";
    std::string error_msg = param_name + " needs to be a power of 2!";
    std::string warning_msg = "";
    utils::AssessParam(validity, param_name, check_mpi_ranks, error_msg,
                       warning_msg, nerrors, do_thorough_checks);

    param_name = "input.restart";
    pp.query(param_name.c_str(), restart_sim);
    utils::AssessParamOK(param_name, restart_sim, do_thorough_checks);

    param_name = "input.get_box_layout_nodes";
    pp.query(param_name.c_str(), get_box_layout_nodes);
    validity = (utils::ErrorState)(utils::IsPowerOfTwo(get_box_layout_nodes) ||
                                   get_box_layout_nodes == 0);
    error_msg = param_name + " needs to be a power of 2!";
    warning_msg = "";
    utils::AssessParam(validity, param_name, get_box_layout_nodes, error_msg,
                       warning_msg, nerrors, do_thorough_checks);

    param_name = "amr.nghost";
    pp.query(param_name.c_str(), nghost);
    validity = utils::ErrorState::OK;
    for (int lev = 0; lev < max_level; ++lev) {
        if (nghost < 0 || nghost >= blocking_factor[lev][0]) {
            validity = utils::ErrorState::ERROR;
            break;
        }
    }
    error_msg = param_name + " needs to be >= 0 and < amr.blocking_factor!";
    warning_msg = "";
    utils::AssessParam(validity, param_name, nghost, error_msg, warning_msg,
                       nerrors, do_thorough_checks);

    param_name = "amr.tagging_on_gpu";
    pp.query(param_name.c_str(), tagging_on_gpu);
    utils::AssessParamOK(param_name, tagging_on_gpu, do_thorough_checks);

    param_name = "amr.coarse_level_grid_size";
    pp.get(param_name.c_str(), coarse_level_grid_size);
    validity = (utils::ErrorState)utils::IsPowerOfTwo(coarse_level_grid_size);
    error_msg = param_name + " needs to be a power of 2!";
    warning_msg = "";
    utils::AssessParam(validity, param_name, coarse_level_grid_size, error_msg,
                       warning_msg, nerrors, do_thorough_checks);

    param_name = "amr.increase_coarse_level_resolution";
    pp.query(param_name.c_str(), increase_coarse_level_resolution);
    validity = increase_coarse_level_resolution ? utils::ErrorState::WARNING
                                                : utils::ErrorState::OK;
    warning_msg = "Will increase coarse level resolution at the beginning.";
    utils::AssessParam(validity, param_name, increase_coarse_level_resolution,
                       "", warning_msg, nerrors, do_thorough_checks);

    param_name = "sim.t_start";
    pp.get(param_name.c_str(), t_start);
    utils::AssessParamOK(param_name, t_start, do_thorough_checks);

    param_name = "sim.t_end";
    pp.get(param_name.c_str(), t_end);
    utils::AssessParamOK(param_name, t_end, do_thorough_checks);

    param_name = "sim.L";
    pp.get(param_name.c_str(), L);
    utils::AssessParamOK(param_name, L, do_thorough_checks);

    param_name = "sim.cfl";
    pp.get(param_name.c_str(), cfl);
    utils::AssessParamOK(param_name, cfl, do_thorough_checks);

    param_name = "sim.gravitational_waves";
    pp.query(param_name.c_str(), with_gravitational_waves);
    utils::AssessParamOK(param_name, with_gravitational_waves,
                         do_thorough_checks);
}

/** @brief Parses all input parameters related to the individual scalar fields.
 */
void Sledgehamr::ParseInputScalars() {
    // Truncation error threshold.
    te_crit.resize(scalar_fields.size());
    double te_crit_default = DBL_MAX;

    amrex::ParmParse pp("");
    std::string param_name = "amr.te_crit";
    pp.query(param_name.c_str(), te_crit_default);
    utils::ErrorState validity = (utils::ErrorState)(te_crit_default > 0);
    std::string error_msg = param_name + " needs to be > 0!";
    std::string warning_msg = "";
    utils::AssessParam(validity, param_name, te_crit_default, error_msg,
                       warning_msg, nerrors, do_thorough_checks);

    shadow_hierarchy = false;
    for (int n = 0; n < scalar_fields.size(); ++n) {
        te_crit[n] = te_crit_default;

        param_name = "amr.te_crit_" + scalar_fields[n]->name;
        pp.query(param_name.c_str(), te_crit[n]);
        validity = (utils::ErrorState)(te_crit[n] > 0);
        error_msg = param_name + " needs to be > 0!";
        warning_msg = "";
        utils::AssessParam(validity, param_name, te_crit[n], error_msg,
                           warning_msg, nerrors, do_thorough_checks);

        if (te_crit[n] != DBL_MAX)
            shadow_hierarchy = true;
    }

    // Dissipation order and strength.
    dissipation_strength.resize(scalar_fields.size());
    double dissipation_default = 0;

    param_name = "sim.dissipation_strength";
    pp.query(param_name.c_str(), dissipation_default);
    validity = (utils::ErrorState)(dissipation_default >= 0);
    error_msg = param_name + " needs to be >= 0!";
    warning_msg = "";
    utils::AssessParam(validity, param_name, dissipation_default, error_msg,
                       warning_msg, nerrors, do_thorough_checks);

    with_dissipation = false;
    for (int n = 0; n < scalar_fields.size(); ++n) {
        dissipation_strength[n] = dissipation_default;

        param_name = "sim.dissipation_strength_" + scalar_fields[n]->name;
        pp.query(param_name.c_str(), dissipation_strength[n]);
        validity = (utils::ErrorState)(dissipation_strength[n] >= 0);
        error_msg = param_name + " needs to be >= 0!";
        warning_msg = "";
        utils::AssessParam(validity, param_name, dissipation_strength[n],
                           error_msg, warning_msg, nerrors, do_thorough_checks);

        if (dissipation_strength[n] > 0)
            with_dissipation = true;
    }

    if (with_dissipation) {
        dissipation_order = nghost;

        param_name = "sim.dissipation_order";
        pp.query(param_name.c_str(), dissipation_order);
        validity = (utils::ErrorState)(dissipation_order == 2 ||
                                       dissipation_order == 3);
        error_msg = "Currently only " + param_name + " = 2 or 3 supported!";
        warning_msg = "";
        utils::AssessParam(validity, param_name, dissipation_order, error_msg,
                           warning_msg, nerrors, do_thorough_checks);
    }
}

/** @brief Loads the spectrum binning from file. Aborts if they are not found.
 * @param   reload  We reload the binning from file. Needed e.g. when we change
 *                  the coarse level resolution.
 */
void Sledgehamr::ReadSpectrumKs(bool reload) {
    if (!spectrum_ks.empty() && !reload)
        return;

    if (reload) {
        spectrum_ks.clear();
        index_to_k.clear();
    }

    std::string filename = SLEDGEHAMR_DATA_PATH;
    filename += "spectra_ks.hdf5";
    std::string sdimN = std::to_string(dimN[0]);

    std::string msg =
        "Sledgehamr::ReadSpectrumKs: Could not find precomputed "
        "spectrum binning!\n Either the path to sledgehamr was set "
        "wrongly during compilation\n (currently set to " SLEDGEHAMR_DATA_PATH
        ")\n or data for a " +
        std::to_string(coarse_level_grid_size) +
        "^3 grid has"
        "not yet been added to the file (github repo only comes "
        "with binnings for a grid up to 512^3).\n Spectrum binning "
        "for a larger grid can generated by running the Jupyter "
        "notebook sledgehamr/notebooks/AddSpectrumBins.ipynb.";

    std::vector<int> nks(1);
    if (!utils::hdf5::Read(filename, {sdimN + "_nks"}, &(nks[0]))) {
        amrex::Abort(msg);
    }

    spectrum_ks.resize(nks[0]);
    if (!utils::hdf5::Read(filename, {sdimN + "_bins"}, &(spectrum_ks[0]))) {
        amrex::Abort(msg);
    }

    if (gravitational_waves == nullptr) {
        return;
    }

    std::string projection =
        std::to_string(gravitational_waves->GetProjectionType());
    index_to_k.resize(dimN[0]);
    if (!utils::hdf5::Read(filename, {sdimN + "_k" + projection},
                           &(index_to_k[0]))) {
        amrex::Abort(msg);
    }
}

}; // namespace sledgehamr
