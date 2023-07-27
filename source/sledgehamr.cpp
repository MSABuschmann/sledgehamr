#include <type_traits>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

#include "sledgehamr.h"
#include "sledgehamr_utils.h"

namespace sledgehamr {

Sledgehamr::Sledgehamr () {
    amrex::Print() << "Starting sledgehamr..." << std::endl;

    ParseInput();

    // Initialize modules.
    time_stepper = new TimeStepper(this);
    io_module = new IOModule(this);

    // Fill various level vectors.
    grid_new.resize(max_level+1);
    grid_old.resize(max_level+1);

    for (int lev=0; lev<=max_level; ++lev) {
        dimN.push_back( coarse_level_grid_size * pow(2,lev) );
        dx.push_back( L/(double)dimN[lev] );
        dt.push_back( dx[lev] * cfl );
    }
}

Sledgehamr::~Sledgehamr() {
    delete level_synchronizer;
    delete time_stepper;
    delete io_module;
}

void Sledgehamr::InitSledgehamr() {
    if (with_gravitational_waves)
        gravitational_waves = new GravitationalWaves(this);

    // Initialize here and not in the SledgeHAMR constructor such that it knows
    // about the number of scalar fields during construction. Necessary so it
    // can initialize boundary conditions.
    level_synchronizer = new LevelSynchronizer(this);
    performance_monitor = new PerformanceMonitor(this);

    ParseInputScalars();

    performance_monitor->Start(performance_monitor->idx_read_input);
    if (restart_sim) {
        io_module->RestartSim();
    } else {
        InitFromScratch( t_start );
    }
    performance_monitor->Stop(performance_monitor->idx_read_input);

    // Initialize project
    Init();
}

void Sledgehamr::Evolve() {
    // Main loop over time.
    while (grid_new[0].t < t_end) {
        // Advance all levels starting at lev=0. This performs an entire
        // shadow/coarse level time step.
        amrex::Print() << std::endl;
        time_stepper->Advance(0);

        // Write any output if requested.
        amrex::Print() << std::endl;
        io_module->Write();
    }

    // Force write at the end of simulation.
    io_module->Write(true);

    amrex::Print() << "Finished!" << std::endl;
}

void Sledgehamr::MakeNewLevelFromScratch(int lev, amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) {
    const int ncomp = scalar_fields.size();

    // Define lowest level from scratch.
    grid_new[lev].define(ba, dm, ncomp, nghost, time);
    grid_old[lev].define(ba, dm, ncomp, nghost);

    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // Fill current level lev with initial state data.
    io_module->FillLevelFromFile(lev);
}

void Sledgehamr::MakeNewLevelFromCoarse(int lev, amrex::Real time,
                                        const amrex::BoxArray& ba,
                                        const amrex::DistributionMapping& dm) {
    const int ncomp = grid_new[lev-1].nComp();
    const int nghost = grid_new[lev-1].nGrow();

    // Define a new level from scratch.
    grid_new[lev].define(ba, dm, ncomp, nghost, time);
    grid_old[lev].define(ba, dm, ncomp, nghost);

    // Fill new level with coarse level data.
    level_synchronizer->FillCoarsePatch(lev, time, grid_new[lev]);
}

void Sledgehamr::RemakeLevel(int lev, amrex::Real time,
                             const amrex::BoxArray& ba,
                             const amrex::DistributionMapping& dm) {
    const int ncomp = grid_new[lev].nComp();
    const int nghost = grid_new[lev].nGrow();

    // Remake new_grid and fill with data.
    LevelData new_state(ba, dm, ncomp, nghost, grid_new[lev].t);
    level_synchronizer->FillPatch(lev, time, new_state);
    std::swap(new_state, grid_new[lev]);
    new_state.clear();

    // Remake old_grid.
    grid_old[lev].clear();
    grid_old[lev].define(ba, dm, ncomp, nghost);
}

void Sledgehamr::ClearLevel(int lev) {
    grid_new[lev].clear();
    grid_old[lev].clear();
}

void Sledgehamr::ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time,
                          int ngrow) {
    // Skip regrid right at the beginning of the sim. Allowed to be overridden
    // if no truncation errors are used.
    if (time == t_start && shadow_hierarchy) return;

    performance_monitor->Start(performance_monitor->idx_tagging, lev);
    utils::sctp timer = utils::StartTimer();

    if (tagging_on_gpu)
        DoErrorEstGpu(lev, tags, time);
    else
        DoErrorEstCpu(lev, tags, time);

    amrex::Print() << "  Tagging took " << utils::DurationSeconds(timer)
                   << "s." << std::endl;
    performance_monitor->Stop(performance_monitor->idx_tagging, lev);
}

void Sledgehamr::CreateShadowLevel() {
    const int ncomp = scalar_fields.size();
    const double time = grid_old[0].t;
    amrex::BoxArray ba = grid_old[0].boxArray();
    ba.coarsen(2);
    amrex::DistributionMapping dm = grid_old[0].DistributionMap();

    shadow_level.define(ba, dm, ncomp, nghost);
    shadow_level_tmp.define(ba, dm, ncomp, nghost, time);

    shadow_level_geom = geom[0];
    shadow_level_geom.coarsen(amrex::IntVect(2,2,2));

    amrex::average_down(grid_old[0], shadow_level_tmp, geom[0],
                        shadow_level_geom, 0, ncomp, refRatio(0));

    time_stepper->integrator->Advance(-1);
}

void Sledgehamr::DoErrorEstCpu(int lev, amrex::TagBoxArray& tags, double time) {
    // Current state.
    const amrex::MultiFab& state = grid_new[lev];

    // State containing truncation errors if they have been calculated.
    const LevelData& state_te = grid_old[lev];

    // Initialize tag counters.
    int ntags_total = 0;
    int ntags_user = 0;
    std::vector<int> ntags_trunc(scalar_fields.size(), 0);

    std::vector<double> params_tag, params_mod;
    SetParamsTagCellForRefinement(params_tag);
    if (shadow_hierarchy)
        SetParamsTruncationModifier(params_mod);
 
    // Loop over boxes and cells.
#pragma omp parallel reduction(+: ntags_total) reduction(+: ntags_user) \
                     reduction(vec_int_plus : ntags_trunc)
    for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
        const amrex::Box& tilebox  = mfi.tilebox();
        const amrex::Array4<double const>& state_fab    = state.array(mfi);
        const amrex::Array4<double const>& state_fab_te = state_te.array(mfi);
        const amrex::Array4<char>& tag_arr = tags.array(mfi);

        // Tag with or without truncation errors.
        if (shadow_hierarchy && state_te.contains_truncation_errors) {
            TagWithTruncationCpu(state_fab, state_fab_te, tag_arr, tilebox,
                                 time, lev, &ntags_total, &ntags_user,
                                 &(ntags_trunc[0]), params_tag, params_mod);
        } else if (shadow_hierarchy) {
            std::string msg = "Trying to tag using truncation errors but no ";
            msg += "truncation errors are computed!";
            amrex::Abort(msg.c_str());
        } else {
            TagWithoutTruncationCpu(state_fab, tag_arr, tilebox, time, lev,
                                    &ntags_total, params_tag);
        }
    }

    // Collect all tags across MPI ranks.
    amrex::ParallelDescriptor::ReduceIntSum(ntags_total, 0);

    if (shadow_hierarchy) {
        amrex::ParallelDescriptor::ReduceIntSum(ntags_user, 0);
        amrex::ParallelDescriptor::ReduceIntSum(&(ntags_trunc[0]),
                                                ntags_trunc.size(), 0);
    }

    // Print statistics.
    long ncells = CountCells(lev);
    double ftotal = (double)ntags_total/(double)ncells;
    double fuser  = (double)ntags_user /(double)ncells;
    amrex::Print()  << "  Tagged cells at level " << lev << ": " << ntags_total
                    << " of " << ncells << " (" << ftotal*100. << "\%)"
                    << std::endl;

    if (shadow_hierarchy) {
        amrex::Print() << "    User-defined tags: " << ntags_user << std::endl;

        for (int i=0; i<scalar_fields.size(); ++i) {
            amrex::Print()  << "    Truncation error tags on "
                            << scalar_fields[i]->name << ": " << ntags_trunc[i]
                            << std::endl;
        }
    }
}

void Sledgehamr::DoErrorEstGpu(int lev, amrex::TagBoxArray& tags, double time) {
    // Current state.
    const amrex::MultiFab& state = grid_new[lev];

    // State containing truncation errors if they have been calculated.
    const amrex::MultiFab& state_te = grid_old[lev];

    std::vector<double> params_tag, params_mod;
    SetParamsTagCellForRefinement(params_tag);
    if (shadow_hierarchy)
        SetParamsTruncationModifier(params_mod);

    // Loop over boxes and cells.
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(state, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& tilebox  = mfi.tilebox();
        const amrex::Array4<double const>& state_fab    = state.array(mfi);
        const amrex::Array4<double const>& state_fab_te = state_te.array(mfi);
        const amrex::Array4<char>& tag_arr = tags.array(mfi);

        // Tag with or without truncation errors.
        if (shadow_hierarchy) {
            TagWithTruncationGpu(state_fab, state_fab_te, tag_arr, tilebox,
                                 time, lev, params_tag, params_mod);
        } else {
            TagWithoutTruncationGpu(state_fab, tag_arr, tilebox, time, lev,
                                    params_tag);
        }
    }

    amrex::Print()  << "  Tagged cells at level " << lev << "." << std::endl;
}

void Sledgehamr::ParseInput() {
    amrex::ParmParse pp_amr("amr");
    pp_amr.query("nghost", nghost);
    pp_amr.query("coarse_level_grid_size", coarse_level_grid_size);
    pp_amr.query("tagging_on_gpu", tagging_on_gpu);

    amrex::ParmParse pp_sim("sim");
    pp_sim.get("t_start", t_start);
    pp_sim.get("t_end", t_end);
    pp_sim.query("restart", restart_sim);
    pp_sim.get("L", L);
    pp_sim.get("cfl", cfl);
    pp_sim.query("gravitational_waves", with_gravitational_waves);
}

void Sledgehamr::ParseInputScalars() {
    // Truncation error threshold.
    te_crit.resize( scalar_fields.size() );
    double te_crit_default = DBL_MAX;

    amrex::ParmParse pp_amr("amr");
    pp_amr.query("te_crit", te_crit_default);

    shadow_hierarchy = false;
    for (int n=0; n<scalar_fields.size(); ++n) {
        te_crit[n] = te_crit_default;
        std::string ident = "te_crit_" + scalar_fields[n]->name;
        pp_amr.query(ident.c_str(), te_crit[n]);

        if (te_crit[n] != DBL_MAX)
            shadow_hierarchy = true;
    }

    // Dissipation order and strength.
    dissipation_strength.resize( scalar_fields.size() );
    double dissipation_default = 0;

    amrex::ParmParse pp_sim("sim");
    pp_sim.query("dissipation_strength", dissipation_default);
    with_dissipation = false;
    for (int n=0; n<scalar_fields.size(); ++n) {
        dissipation_strength[n] = dissipation_default;
        std::string ident = "dissipation_strength_" + scalar_fields[n]->name;
        pp_sim.query(ident.c_str(), dissipation_strength[n]);

        if (dissipation_strength[n] > 0)
            with_dissipation = true;
    }

    if (with_dissipation) {
        dissipation_order = nghost;
        pp_sim.query("dissipation_order", dissipation_order);

        if (dissipation_order < 2) {
            amrex::Print() << "#WARNING: Dissipation order "
                           << dissipation_order
                           << " not supported! Dissipation disabled."
                           << std::endl;
        }

        if (dissipation_order > 3) {
            amrex::Print() << "#WARNING: Dissipation order "
                           << dissipation_order
                           << " not supported! Dissipation defaults to order=3."
                           << std::endl;
            dissipation_order = 3;
        }
    }
}

void Sledgehamr::ReadSpectrumKs() {
    if (!spectrum_ks.empty())
        return;

    std::string filename = SLEDGEHAMR_DATA_PATH;
    filename += "spectra_ks.hdf5";
    std::string sdimN = std::to_string(dimN[0]);

    std::vector<int> nks(1);
    IOModule::ReadFromHDF5(filename, {sdimN+"_nks"}, &(nks[0]));

    spectrum_ks.resize(nks[0]);
    IOModule::ReadFromHDF5(filename, {sdimN}, &(spectrum_ks[0]));
}

}; // namespace sledgehamr
