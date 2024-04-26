#ifndef SLEDGEHAMR_SLEDGEHAMR_H_
#define SLEDGEHAMR_SLEDGEHAMR_H_

#include <AMReX_AmrCore.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "kernels.h"
#include "macros.h"

#include "gravitational_waves.h"
#include "io_module.h"
#include "level_data.h"
#include "level_synchronizer.h"
#include "local_regrid/local_regrid.h"
#include "output_types/checkpoint.h"
#include "performance_monitor.h"
#include "projection.h"
#include "scalars.h"
#include "spectrum.h"
#include "time_stepper.h"

namespace sledgehamr {

class LevelSynchronizer;
class LocalRegrid;
class TimeStepper;
class Integrator;
class IOModule;
class Projection;
class Spectrum;
class GravitationalWaves;
class Checkpoint;
class PerformanceMonitor;

/** @brief Abstract base class for all derived projects. Combines all the
 *         ingredients to make this code work.
 */
class Sledgehamr : public amrex::AmrCore {
    friend class LevelSynchronizer;
    friend class LocalRegrid;
    friend class TimeStepper;
    friend class Integrator;
    friend class IOModule;
    friend class Projection;
    friend class Spectrum;
    friend class GravitationalWaves;
    friend class Checkpoint;
    friend class PerformanceMonitor;

  public:
    Sledgehamr();
    void InitSledgehamr();
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
     * @param   lev         Currently level.
     * @param   dt          Time step size.
     * @param   dx          Grid spacing.
     */
    virtual void FillRhs(amrex::MultiFab &rhs_mf,
                         const amrex::MultiFab &state_mf, const double time,
                         const int lev, const double dt, const double dx) = 0;

    /** @brief Like FillRhs but the result will be added to rhs_mf rather than
     *         set.
     * @param   rhs_mf      Empty MultiFab to be filled with RHS.
     * @param   state_mf    State from which the RHS is to be computed.
     * @param   time        Current time.
     * @param   lev         Currently level.
     * @param   dt          Time step size.
     * @param   dx          Grid spacing.
     * @param   weight      rhs = weight*rhs + ...
     */
    virtual void FillAddRhs(amrex::MultiFab &rhs_mf,
                            const amrex::MultiFab &state_mf, const double time,
                            const int lev, const double dt, const double dx,
                            const double weight) = 0;

    /** @brief Returns the box length L.
     */
    double GetL() const { return L; };

    /** @brief Returns the grid spacing dx at a given level.
     * @param   lev Level.
     */
    double GetDx(const int lev) const { return dx[lev]; }

    /** @brief Returns the time step size dt at a given level.
     * @param   lev Level.
     */
    double GetDt(const int lev) const { return dt[lev]; }

    /** @brief Returns the number of (potential) cells along an axis at a given
     *         level.
     * @param   lev Level.
     */
    int GetDimN(const int lev) const { return dimN[lev]; }

    /** @brief Returns the maximum number of allowed levels.
     */
    int GetMaxLevel() const { return max_level; }

    /** @brief Returns the number of current levels.
     */
    int GetFinestLevel() const { return finest_level; }

    /** @brief Returns the vector of amrex::Geometry's.
     */
    amrex::Vector<amrex::Geometry> &GetGeometry() { return geom; }

    /** @brief Returns reference to current data at a given level.
     * @param   lev Level.
     */
    LevelData &GetLevelData(const int lev) { return grid_new[lev]; }

    /** @brief Returns reference to old data at a given level.
     * @param   lev Level.
     */
    LevelData &GetOldLevelData(const int lev) { return grid_old[lev]; }

    /** @brief Returns the blocking factor at a given level.
     * @param   lev Level.
     */
    int GetBlockingFactor(const int lev) const {
        return blocking_factor[lev][0];
    }

    /** @brief Returns the name of a scalar field.
     * @param   comp    Scalar field component.
     */
    std::string GetScalarFieldName(const int comp) const {
        return scalar_fields[comp]->name;
    }

    /** @brief Pointer to synchronization module.
     */
    std::unique_ptr<LevelSynchronizer> level_synchronizer;

    /** @brief Pointer to I/O module.
     */
    std::unique_ptr<IOModule> io_module;

    /** @brief Number of ghost cells.
     */
    int nghost = 0;

    /** @brief Whether we are running with gravitional waves.
     */
    bool with_gravitational_waves = false;

    /** @brief Whether we are carefully checking all input parameters. If 'true'
     *         we will not start the actual simulation.
     */
    bool do_thorough_checks = false;

    /** @brief  If we are doing thorough checks of parameters the node layout
     *          can be different than what we want to use during the simulation
     *          run. This variable contains the number of nodes we want to
     *          ultimately use during simulation.
     */
    int check_mpi_ranks = 0;

    /** @brief Number of parameter errors we encountered during initialization.
     */
    int nerrors = 0;

  protected:
    virtual void
    MakeNewLevelFromScratch(int lev, amrex::Real time,
                            const amrex::BoxArray &ba,
                            const amrex::DistributionMapping &dm) override;
    virtual void
    MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray &ba,
                           const amrex::DistributionMapping &dm) override;
    virtual void RemakeLevel(int lev, amrex::Real time,
                             const amrex::BoxArray &ba,
                             const amrex::DistributionMapping &dm) override;

    virtual void ClearLevel(int lev) override;
    virtual void ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time,
                          int ngrow) override;

    /** @brief Virtual function that loop over a state to tags cells for
     *         refinement. Includes user-defined tags and truncation error tags.
     *         Will automatically be overriden by the project class.
     *         Work is performed on CPUs.
     * @param   state_fab       Data.
     * @param   state_fab_te    Truncation errors.
     * @param   tagarr          Tag status.
     * @param   tilebox         Current (tile)box.
     * @param   time            Current time.
     * @param   lev             Current level.
     * @param   ntags_user      Counts number of user-defined tags.
     * @param   params_mod      User-defined parameters.
     */
    virtual void
    TagWithTruncationCpu(const amrex::Array4<const double> &state_fab,
                         const amrex::Array4<const double> &state_fab_te,
                         const amrex::Array4<char> &tagarr,
                         const amrex::Box &tilebox, double time, int lev,
                         long *ntags_total, long *ntags_user, long *ntags_trunc,
                         const std::vector<double> &params_tag,
                         const std::vector<double> &params_mod) = 0;

    /** @brief Same as TagWithTruncationCpu but performs work on GPUs.
     */
    virtual void
    TagWithTruncationGpu(const amrex::Array4<const double> &state_fab,
                         const amrex::Array4<const double> &state_fab_te,
                         const amrex::Array4<char> &tagarr,
                         const amrex::Box &tilebox, double time, int lev,
                         const std::vector<double> &params_tag,
                         const std::vector<double> &params_mod) = 0;

    /** @brief Same as TagWithTruncationCpu but does not include
     *         truncation error tags.
     */
    virtual void
    TagWithoutTruncationCpu(const amrex::Array4<const double> &state_fab,
                            const amrex::Array4<char> &tagarr,
                            const amrex::Box &tilebox, double time, int lev,
                            long *ntags_total,
                            const std::vector<double> &params) = 0;

    /** @brief Same as TagWithTruncationGpu but does not include
     *         truncation error tags.
     */
    virtual void
    TagWithoutTruncationGpu(const amrex::Array4<const double> &state_fab,
                            const amrex::Array4<char> &tagarr,
                            const amrex::Box &tilebox, double time, int lev,
                            const std::vector<double> &params) = 0;

    /** @brief Initialize project specific details. To be overriden by each
     *         project if required.
     */
    virtual void Init(){};

    /** @brief Function to check whether a given level is allowed to be created
     *         by a given time. Will always allow a level to be created by
     *         default.
     * @param   lev     Level to be created.
     * @param   time    Time before which the level would be created.
     * @return  If the given level is allowed to be created.
     */
    bool DoCreateLevelIf(const int lev, const double time) {
        // We always need to have a shadow or coarse level.
        return lev <= 0 ? true : CreateLevelIf(lev, time);
    }

    /** @brief Same as DoCreateLevelIf but can be overriden by the project
               class.
     */
    virtual bool CreateLevelIf(const int lev, const double time) {
        return true;
    };

    /** @brief Any work that should be done before a coarse level time step.
     *         To be overriden by the project class.
     */
    virtual void BeforeTimestep(const double time){};

    /** @brief Stopping condition of the simulation. Can be overriden by the
     *         project class.
     * @param   time    Current time.
     */
    virtual bool StopRunning(const double time) { return time >= t_end; }

    /** @brief Function to be overriden by the project class to set user-defined
     *         parameters that will be passed to the Rhs computation.
     * @param   params  Parameters to be set.
     * @param   time    Current time.
     * @param   lev     Current level.
     */
    virtual void SetParamsRhs(std::vector<double> &params, const double time,
                              const int lev){};

    /** @brief Function to be overriden by the project class to set user-defined
     *         parameters that will be passed to the Rhs computation of
     *         gravitational waves.
     * @param   params  Parameters to be set.
     * @param   time    Current time.
     * @param   lev     Current level.
     */
    virtual void SetParamsGravitationalWaveRhs(std::vector<double> &params,
                                               const double time,
                                               const int lev){};

    /** @brief Function to be overriden by the project class to set user-defined
     *         parameters that will be passed to the tagging procedure.
     * @param   params  Parameters to be set.
     * @param   time    Current time.
     * @param   lev     Current level.
     */
    virtual void SetParamsTagCellForRefinement(std::vector<double> &params,
                                               const double time,
                                               const int lev){};

    /** @brief Function to be overriden by the project class to set user-defined
     *         parameters that will be available when modifying the truncation
     *         error threshold.
     * @param   params  Parameters to be set.
     * @param   time    Current time.
     * @param   lev     Current level.
     */
    virtual void SetParamsTruncationModifier(std::vector<double> &params,
                                             const double time,
                                             const int lev){};

    /** @brief Function to be overriden by the project class to set user-defined
     *         parameters that will be passed to the spectrum computation.
     * @param   params  Parameters to be set.
     * @param   time    Current time.
     * @param   lev     Current level.
     */
    virtual void SetParamsSpectra(std::vector<double> &params,
                                  const double time){};

    /** @brief Function to be overriden by the project class to set user-defined
     *         parameters that will be passed to the computation of projections.
     * @param   params  Parameters to be set.
     * @param   time    Current time.
     * @param   lev     Current level.
     */
    virtual void SetParamsProjections(std::vector<double> &params,
                                      const double time){};

    void CreateShadowLevel();

    /** @brief Pointer to the time stepper submodule.
     */
    std::unique_ptr<TimeStepper> time_stepper;

    /** @brief Pointer to the gravitational wave submodule.
     */
    std::unique_ptr<GravitationalWaves> gravitational_waves;

    /** @brief Pointer to the performance monitor.
     */
    std::unique_ptr<PerformanceMonitor> performance_monitor;

    /** @brief Holds the most recent simulation data for all levels.
     */
    std::vector<LevelData> grid_new;

    /** @brief Holds old simulation data from the previous time step for all
     *         levels.
     */
    std::vector<LevelData> grid_old;

    /** @brief Shadow level in case we need it.
     */
    bool shadow_hierarchy = false;

    /** @brief Actual shadow level data when we need it.
     */
    LevelData shadow_level, shadow_level_tmp;

    /** @brief Geometry of shadow level.
     */
    amrex::Geometry shadow_level_geom;

    /** @brief Holds pointers to all simulated scalar fields. We assume
     *         ownership of all objects within this vector - one day we upgrade
     *         to a smart pointer.
     */
    std::vector<ScalarField *> scalar_fields;

    /** @brief Start times of the simulation.
     */
    double t_start;

    /** @brief End times of the simulation.
     */
    double t_end;

    /** @brief Time step size at each level.
     */
    std::vector<double> dt;

    /** @brief Grid spacing at each level.
     */
    std::vector<double> dx;

    /** @brief CFL criteria.
     */
    double cfl;

    /** @brief Box length.
     */
    double L;

    /** @brief Number of potential cells in each direction for each level.
     */
    std::vector<int> dimN;

    /** @brief Number of coarse level cells in each direction.
     */
    int coarse_level_grid_size;

    /** @brief Vector of truncation error thresholds.
     */
    std::vector<double> te_crit;

    /** @brief Unique bins for spectrum calculation.
     */
    std::vector<int> spectrum_ks;

    /** @brief Pre-computed conversion of index to k for gravitational wave
     * spectrum calculation.
     */
    std::vector<double> index_to_k;

    /** @brief  Vector of respective dissipation strenths.
     */
    std::vector<double> dissipation_strength;

    /** @brief  Whether we want to simulate with Kreiss-Oliger dissipation.
     */
    bool with_dissipation = false;

    /** @brief Order of Kreiss-Oliger dissipation. Needs to be 2 or 3 for
     *         the dissipation to be enabled.
     */
    int dissipation_order = 0;

    /** @brief Whether we are restarting the sim or from scratch.
     */
    bool restart_sim = false;

  private:
    void DoErrorEstCpu(int lev, amrex::TagBoxArray &tags, double time);
    void DoErrorEstGpu(int lev, amrex::TagBoxArray &tags, double time);

    void ParseInput();
    void ParseInputScalars();
    void ReadSpectrumKs(bool reload = false);
    void DoPrerunChecks();
    void DetermineBoxLayout();

    /** @brief Whether tagging should be performed on gpu if possible.
     */
    bool tagging_on_gpu = false;

    /** @brief Whether we actually want to perform a simulation or just checking
     *         parameters or other things.
     */
    bool no_simulation = false;

    /** @brief If this is larger zero we want to save the coarse level box
     *         layout if we were to run with this amount of nodes.
     */
    int get_box_layout_nodes = 0;

    /** @brief Whether we want to increase the coarse level resolution once.
     */
    bool increase_coarse_level_resolution = false;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMR_H_
