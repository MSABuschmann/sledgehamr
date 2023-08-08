#ifndef SLEDGEHAMR_SLEDGEHAMR_H_
#define SLEDGEHAMR_SLEDGEHAMR_H_

#include <AMReX_AmrCore.H>

#include "kernels.h"
#include "macros.h"
#include "level_data.h"
#include "level_synchronizer.h"
#include "local_regrid.h"
#include "time_stepper.h"
#include "io_module.h"
#include "scalars.h"
#include "projection.h"
#include "spectrum.h"
#include "gravitational_waves.h"
#include "checkpoint.h"
#include "performance_monitor.h"

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

/** @brief Base class for all derived projects. Combines all the ingredients to
 *         make this code work.
 */
class Sledgehamr : public amrex::AmrCore {
    // Give submodules access to data.
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
    /** @brief Creates instances of submodules and reads input parameters.
     */
    Sledgehamr();

    virtual ~Sledgehamr();

    /** @brief Initalizes data from scratch or from checkpoint file.
     */
    void InitSledgehamr();

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
     * @param   lev         Currently level.
     * @param   dt          Time step size.
     * @param   dx          Grid spacing.
     */
    virtual void FillRhs(amrex::MultiFab& rhs_mf,
                         const amrex::MultiFab& state_mf, const double time,
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
    virtual void FillAddRhs(amrex::MultiFab& rhs_mf,
                            const amrex::MultiFab& state_mf, const double time,
                            const int lev, const double dt, const double dx,
                            const double weight) = 0;

    /** @brief Pointer to synchronization module.
     */
    LevelSynchronizer* level_synchronizer;
    IOModule* io_module;

    /** @brief Number of ghost cells.
     */
    int nghost = 0;

    bool with_gravitational_waves = false;

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
            int* ntags_trunc, const std::vector<double>& params_tag,
            const std::vector<double>& params_mod) = 0;

    virtual void TagWithTruncationGpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<const double>& state_fab_te,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev, const std::vector<double>& params_tag,
            const std::vector<double>& params_mod) = 0;

    /** @brief Same as TagWithTruncation but does not include
     *         truncation error tags.
     */
    virtual void TagWithoutTruncationCpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev, int* ntags_total,
            const std::vector<double>& params) = 0;

    virtual void TagWithoutTruncationGpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev, const std::vector<double>& params) = 0;

    /** @brief Initialize project specific details.
     */
    virtual void Init() {};

    /** @brief Function to check whether a given level is allowed to be created
     *         by a given time. Can be overridden by the project class. Will
     *         always allow a level to be created by default.
     * @param   lev     Level to be created.
     * @param   time    Time before which the level would be created.
     * @return  If the given level is allowed to be created.
     */
    virtual bool CreateLevelIf(const int lev, const double time) {
        return true;
    };

    virtual void SetParamsRhs(std::vector<double>& params) {};
    virtual void SetParamsGravitationalWaveRhs(std::vector<double>& params) {};
    virtual void SetParamsTagCellForRefinement(std::vector<double>& params) {};
    virtual void SetParamsTruncationModifier(std::vector<double>& params) {};
    virtual void SetParamsSpectra(std::vector<double>& params) {};
    virtual void SetParamsProjections(std::vector<double>& params) {};

    /** @brief Creates a shadow level and evolves it by one time step. Needed
     *         to compute truncation errors on the coarse level.
     */
    void CreateShadowLevel();

    virtual void BeforeTimestep(const double time) {};

    /** @brief Pointer to sub-modules.
     */
    TimeStepper* time_stepper;
    GravitationalWaves* gravitational_waves;
    PerformanceMonitor* performance_monitor;

    /** @brief Holds the actual simulation data for all levels at two different
     *         states in time.
     */
    std::vector<LevelData> grid_old, grid_new;

    /** @brief Shadow level in case we need it.
     */
    bool shadow_hierarchy = false;
    LevelData shadow_level, shadow_level_tmp;
    amrex::Geometry shadow_level_geom;

    /** @brief Holds pointers to all simulated scalar fields.
     */
    std::vector<ScalarField*> scalar_fields;

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

    /** @brief TODO
     */
    std::vector<int> spectrum_ks;

    std::vector<double> dissipation_strength;
    bool with_dissipation = false;
    int dissipation_order = 0;

    bool restart_sim = false;

private:

    /* @brief Do ErrorEst on either CPU or GPU.
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

    /** @brief TODO
     */
    void ReadSpectrumKs();

    /* @brief Whether tagging should be performed on gpu if possible.
     */
    bool tagging_on_gpu = false;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMR_H_
