#ifndef SLEDGEHAMR_SLEDGEHAMR_H_
#define SLEDGEHAMR_SLEDGEHAMR_H_

#include <AMReX_AmrCore.H>

#include "macros.h"
#include "level_data.h"
#include "level_synchronizer.h"
#include "time_stepper.h"
#include "integrator.h"
#include "io_module.h"
#include "scalars.h"
#include "kernels.h"

namespace sledgehamr {

template <typename T> class LevelSynchronizer;
template <typename T> class TimeStepper;
template <typename T> class Integrator;
template <typename T> class IOModule;

/** @brief Base class for all derived projects. Combines all the ingredients to
 *         make this code work.
 */
template <typename T>
class Sledgehamr : public amrex::AmrCore {
    // Give submodules access to data.
    friend class LevelSynchronizer<T>;
    friend class TimeStepper<T>;
    friend class Integrator<T>;
    friend class IOModule<T>;

  public:
    /** @brief Creates instances of submodules and reads input parameters.
     */
    Sledgehamr();

    virtual ~Sledgehamr();

    /** @brief Initalizes data from scratch or from checkpoint file.
     */
    void Init();

    /** @brief Starts the evolution
     */
    void Evolve();

    /** @brief Virtual function that loops over a state to fill the RHS. Has to
     *         be defined by derived project class, though this will be hidden
     *         from the user. Reason for this is so that the RHS function,
     *         declared also by the project class, can be inlined within this
     *         function, which it otherwise couldn't as it is not clear at
     *         compile time which project will be run.
     * @param   rhs_mf      Empty MultiFab to be filled with RHS.
     * @param   state_mf    State from which the RHS is to be computed.
     * @param   time        Current time.
     * @param   geom        Geometry of the current level.
     * @param   lev         Currently level.
     */
    AMREX_FORCE_INLINE
    void FillRhs(amrex::MultiFab& rhs_mf,
                         const amrex::MultiFab& state_mf, const double time,
                         const amrex::Geometry& geom, int lev);

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    virtual void Rhs(const amrex::Array4<double>& rhs,
                     const amrex::Array4<const double>& state,
                     const int i, const int j, const int k, const int lev,
                     const double time, const double dt, const double dx) = 0;


  protected:
    /** @brief Make a new level from scratch using provided BoxArray and
     *         DistributionMapping. Only used during initialization. Overrides
     *         the pure virtual function in amrex::AmrCore.
     * @param   lev     Level to be created.
     * @param   time    Time of new grid.
     * @param   ba      New amrex::BoxArray.
     * @param   dm      New amrex::DistributionMapping.
     */
    virtual void MakeNewLevelFromScratch(int lev, amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm)
                                         override;

    /** @brief Make a new level using provided BoxArray and DistributionMapping,
     *         and fills it with interpolated coarse level data. Overrides the
     *         pure virtual function in amrex::AmrCore.
     * @param   lev     Level to be created.
     * @param   time    Time of new grid.
     * @param   ba      New amrex::BoxArray.
     * @param   dm      New amrex::DistributionMapping.
     */
    virtual void MakeNewLevelFromCoarse(int lev, amrex::Real time,
                                        const amrex::BoxArray& ba,
                                        const amrex::DistributionMapping& dm)
                                        override;

    /** @brief Remake a new level using provided BoxArray and
     *         DistributionMapping, and fills it with interpolated coarse level
     *         data. Overrides the pure virtual function in amrex::AmrCore.
     * @param   lev     Level to be remade.
     * @param   time    Time of new grid.
     * @param   ba      New amrex::BoxArray.
     * @param   dm      New amrex::DistributionMapping.
     */
    virtual void RemakeLevel(int lev, amrex::Real time,
                             const amrex::BoxArray& ba,
                             const amrex::DistributionMapping& dm) override;

    /** @brief Delete level data. Overrides the pure virtual function in
     *         amrex::AmrCore.
     * @param   lev Level to be deleted.
     */
    virtual void ClearLevel(int lev) override;

    /** @brief Tag cells for refinement. Overrides the pure virtual function in
     *         amrex::AmrCore.
     * @param   lev         Level on which cells are tagged.
     * @param   time        Time of said level.
     * @param   ngrow       Grid growth factor.
     * @param   ntags_user  Counts number of user-defined tags.
     */
    virtual void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time,
                          int ngrow) override;

    /** @brief Virtual function that loop over a state to tags cells for
     *         refinement. Includes user-defined tags and truncation error tags.
     *         Will automatically be overriden by the project class.
     * @param   state_fab       Data.
     * @param   state_fab_te    Truncation errors.
     * @param   tagarr          Tag status.
     * @param   tilebox         Current (tile)box.
     * @param   time            Current time.
     * @param   lev             Current level.
     * @param   ntags_user      Counts number of user-defined tags.
     */
    virtual void TagWithTruncationCpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<const double>& state_fab_te,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev, int* ntags_total, int* ntags_user,
            int* ntags_trunc) {} ;

    virtual void TagWithTruncationGpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<const double>& state_fab_te,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev) {};

    /** @brief Same as TagWithTruncation but does not include
     *         truncation error tags.
     */
    virtual void TagWithoutTruncationCpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev, int* ntags_total) {};

    virtual void TagWithoutTruncationGpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev) {};

    /** @brief Function to check whether a given level is allowed to be created
     *         by a given time. Can be overridden by the project class. Will
     *         always allow a level to be created by default.
     * @param   lev     Level to be created.
     * @param   time    Time before which the level would be created.
     * @return  If the given level is allowed to be created.
     */
    virtual bool CreateLevelIf(int lev, double time) {
        return true;
    };

    /** @brief Instance to perform operations between two levels.
     */
    LevelSynchronizer<T>* level_synchronizer;

    /** @brief Instance to perform the sub-cycling in time and that handles
     *         regrid calls.
     */
    TimeStepper<T>* time_stepper;

    /** @brief Module that handles all IO operations (with the exception of
     *         parsing the inputs file.
     */
    IOModule<T>* io_module;

    /** @brief Holds the actual simulation data for all levels at two different
     *         states in time.
     */
    std::vector<LevelData> grid_old, grid_new;

    /** @brief Flag whether simulation should use a shadow hierarchy. To be set
     *         in inputs file.
     */
    bool shadow_hierarchy = false;

    /** @brief Holds pointers to all simulated scalar fields.
     */
    std::vector<ScalarField*> scalar_fields;

    /** @brief Number of ghost cells.
     */
    int nghost = 0;

    /** @brief Start and end times of the simulation.
     */
    double t_start, t_end;

    /** @brief Time step size and grid spacing at each level.
     */
    std::vector<double> dt, dx;

    /** @brief CFL criteria.
     */
    double cfl;

    /** @brief Box length.
     */
    double L;

    /** @brief Number of cells in each direction for each level.
     */
    std::vector<int> dimN;

    /** @brief Number of coarse level cells in each direction.
     */
    int coarse_level_grid_size;

    /** @brief Truncation error thresholds.
     */
    std::vector<double> te_crit;

private:

    /* @brief TODO
     */
    void DoErrorEstCpu(int lev, amrex::TagBoxArray& tags, double time);
    void DoErrorEstGpu(int lev, amrex::TagBoxArray& tags, double time);

    /** @brief Parse various input parameters from inputs file.
     */
    void ParseInput ();

    /** @brief Parse various scalar dependent input parameters from inputs file.
     *         This cannot be done in ParseInput as the project class has not
     *         been constructed at this point, hence we have no knowledge about
     *         the scalars quite yet. This function will be called in Init()
     *         instead.
     */
    void ParseInputScalars();

    /* @brief Whether tagging should be performed on gpu if possible.
     */
    bool tagging_on_gpu = false;
};

template <typename T>
AMREX_FORCE_INLINE
void Sledgehamr<T>::FillRhs(amrex::MultiFab& rhs_mf,
                         const amrex::MultiFab& state_mf,
                         const double time, const amrex::Geometry& geom,
                         int lev) {
    double l_dt = dt[lev];
    double l_dx = dx[lev];

    void (*rhs)(const amrex::Array4<double>&,
                const amrex::Array4<const double>&, const int, const int,
                const int, const int, const double, const double, const double);
    T* psim = static_cast<T*>(this);
    rhs = (void*)&T::Rhs;

#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const amrex::Array4<double>& rhs_fab = rhs_mf.array(mfi);
        const amrex::Array4<double const>& state_fab = state_mf.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rhs(rhs_fab, state_fab, i, j, k, lev, time, l_dt, l_dx);
            //static_cast<T*>(this)->Rhs(rhs_fab, state_fab, i, j, k, lev, time, l_dt, l_dx);
        });
    }
};

template <typename T>
Sledgehamr<T>::Sledgehamr () {
    amrex::Print() << "Starting sledgehamr..." << std::endl;

    ParseInput();

    // Initialize modules.
    time_stepper = new TimeStepper(this);
    io_module = new IOModule(this);

    // Fill various level vectors.
    grid_new.resize(max_level+1);
    grid_old.resize(max_level+1);

    for (int lev=0; lev<=max_level; ++lev) {
        dimN.push_back( coarse_level_grid_size * pow(2,lev-shadow_hierarchy) );
        dx.push_back( L/(double)dimN[lev] );
        dt.push_back( dx[lev] * cfl );
    }
}

template <typename T>
Sledgehamr<T>::~Sledgehamr() {
    delete level_synchronizer;
    delete time_stepper;
    delete io_module;
}

template <typename T>
void Sledgehamr<T>::Init() {
    // Initialize here and not in the SledgeHAMR constructor such that it knows
    // about the number of scalar fields during construction. Necessary so it
    // can initialize boundary conditions.
    level_synchronizer = new LevelSynchronizer(this);

    ParseInputScalars();

    // TODO: Check for checkpoint file etc.
    InitFromScratch( t_start );
}

template <typename T>
void Sledgehamr<T>::Evolve() {
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
}

template <typename T>
void Sledgehamr<T>::MakeNewLevelFromScratch(int lev, amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) {
    const int ncomp = scalar_fields.size();

    // Define lowest level from scratch.
    grid_new[lev].define(ba, dm, ncomp, nghost, time);
    grid_old[lev].define(ba, dm, ncomp, nghost);

    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // If shadow hierarchy is used, above level is the shadow level. We now need
    // to also make the coarse level.
    if (shadow_hierarchy && lev == 0) {
        ++lev;
        ++finest_level;

        // Create local copy of ba since ba is const.
        amrex::BoxArray rba = ba;
        rba.refine(2);

        grid_new[lev].define(rba, dm, ncomp, nghost, time);
        grid_old[lev].define(rba, dm, ncomp, nghost);

        // These are already set for shadow level.
        SetBoxArray(lev, rba);
        SetDistributionMap(lev, dm);
    }

    // Fill current level lev with initial state data.
    io_module->FillLevelFromFile(lev);

    // Fill shadow level with data from coarse level.
    if (shadow_hierarchy)
        level_synchronizer->AverageDownTo(0);
}

template <typename T>
void Sledgehamr<T>::MakeNewLevelFromCoarse(int lev, amrex::Real time,
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

template <typename T>
void Sledgehamr<T>::RemakeLevel(int lev, amrex::Real time,
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

template <typename T>
void Sledgehamr<T>::ClearLevel(int lev) {
    grid_new[lev].clear();
    grid_old[lev].clear();
}

template <typename T>
void Sledgehamr<T>::ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time,
                          int ngrow) {
    // Skip regrid right at the beginning of the sim. Allowed to be overridden
    // if no truncation errors are used (TODO).
    if (time == t_start) return;

    if (tagging_on_gpu)
        DoErrorEstGpu(lev, tags, time);
    else
        DoErrorEstCpu(lev, tags, time);
}

template <typename T>
void Sledgehamr<T>::DoErrorEstCpu(int lev, amrex::TagBoxArray& tags, double time) {
    // Current state.
    const amrex::MultiFab& state = grid_new[lev];

    // State containing truncation errors if they have been calculated.
    const amrex::MultiFab& state_te = grid_old[lev];

    // Initialize tag counters.
    int ntags_total = 0;
    int ntags_user = 0;
    std::vector<int> ntags_trunc(scalar_fields.size(), 0);

    // Loop over boxes and cells.
#pragma omp parallel reduction(+: ntags_total) reduction(+: ntags_user)\
                     reduction(vec_int_plus : ntags_trunc)
    for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
        const amrex::Box& tilebox  = mfi.tilebox();
        const amrex::Array4<double const>& state_fab    = state.array(mfi);
        const amrex::Array4<double const>& state_fab_te = state_te.array(mfi);
        const amrex::Array4<char>& tag_arr = tags.array(mfi);

        // Tag with or without truncation errors.
        if (shadow_hierarchy) {
            TagWithTruncationCpu(state_fab, state_fab_te, tag_arr, tilebox,
                                 time, lev, &ntags_total, &ntags_user,
                                 &(ntags_trunc[0]));
        } else {
            TagWithoutTruncationCpu(state_fab, tag_arr, tilebox, time, lev,
                                    &ntags_total);
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

template <typename T>
void Sledgehamr<T>::DoErrorEstGpu(int lev, amrex::TagBoxArray& tags, double time) {
    // Current state.
    const amrex::MultiFab& state = grid_new[lev];

    // State containing truncation errors if they have been calculated.
    const amrex::MultiFab& state_te = grid_old[lev];

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
                                 time, lev);
        } else {
            TagWithoutTruncationGpu(state_fab, tag_arr, tilebox, time, lev);
        }
    }

    amrex::Print()  << "  Tagged cells at level " << lev << "." << std::endl;
}

template <typename T>
void Sledgehamr<T>::ParseInput() {
    amrex::ParmParse pp_amr("amr");
    pp_amr.query("nghost", nghost);
    pp_amr.query("shadow_hierarchy", shadow_hierarchy);
    pp_amr.query("coarse_level_grid_size", coarse_level_grid_size);
    pp_amr.query("tagging_on_gpu", tagging_on_gpu);

    amrex::ParmParse pp_sim("sim");
    pp_sim.get("t_start", t_start);
    pp_sim.get("t_end", t_end);
    pp_sim.get("L", L);
    pp_sim.get("cfl", cfl);
}

template <typename T>
void Sledgehamr<T>::ParseInputScalars() {
    te_crit.resize( scalar_fields.size() );
    double te_crit_default = 1e99;

    amrex::ParmParse pp("amr");
    pp.query("te_crit", te_crit_default);
    for (int n=0; n<scalar_fields.size(); ++n) {
        te_crit[n] = te_crit_default;
        std::string ident = "te_crit_" + scalar_fields[n]->name;
        pp.query(ident.c_str(), te_crit[n]);
    }
}

}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMR_H_
