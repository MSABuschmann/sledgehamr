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

class LevelSynchronizer;
class TimeStepper;
class Integrator;
class IOModule;

/** @brief Base class for all derived projects. Combines all the ingredients to
 *         make this code work.
 */
class Sledgehamr : public amrex::AmrCore {
    // Give submodules access to data.
    friend class LevelSynchronizer;
    friend class TimeStepper;
    friend class Integrator;
    friend class IOModule;

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
            int* ntags_trunc) = 0;

    virtual void TagWithTruncationGpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<const double>& state_fab_te,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev) = 0;

    /** @brief Same as TagWithTruncation but does not include
     *         truncation error tags.
     */
    virtual void TagWithoutTruncationCpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev, int* ntags_total) = 0;

    virtual void TagWithoutTruncationGpu(
            const amrex::Array4<const double>& state_fab,
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox,
            double time, int lev) = 0;

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
    virtual void FillRhs(amrex::MultiFab& rhs_mf,
                         const amrex::MultiFab& state_mf, const double time,
                         const amrex::Geometry& geom, int lev) = 0;

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
    LevelSynchronizer* level_synchronizer;

    /** @brief Instance to perform the sub-cycling in time and that handles
     *         regrid calls.
     */
    TimeStepper* time_stepper;

    /** @brief Module that handles all IO operations (with the exception of
     *         parsing the inputs file.
     */
    IOModule* io_module;

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

}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMR_H_
