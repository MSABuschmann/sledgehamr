#ifndef SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_
#define SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_

#include <AMReX_BCUtil.H>
#include <AMReX_Interpolater.H>

#include "sledgehamr.h"
#include "level_data.h"

namespace sledgehamr {

/** @brief Struct with overloaded operator to handle boundary conditions. Empty
 *         because we do not have boundary conditions beyond periodic.
 */
struct NullFill {
    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& /*iv*/,
                     amrex::Array4<amrex::Real> const& /*dest*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/,
                     const amrex::Real /*time*/, const amrex::BCRec* /*bcr*/,
                     const int /*bcomp*/, const int /*orig_comp*/) const {}
};

class Sledgehamr;

/** @brief This class handles all operations between two levels such as
 *         averaging down, interpolation to fine, filling of ghost cells, etc.
 *         Class is friend of SledgeHAMR.
 */
class LevelSynchronizer {
  public:
    LevelSynchronizer(Sledgehamr* owner);

    /** @brief Fills LevelData with information from a coarse level. This is
     *         used e.g. when a new level of refinement is added.
     * @param   lev     New level to be filled with data from lev=1.
     * @param   time    Time of new level.
     * @param   ld      New level data.
     */
    void FillCoarsePatch(const int lev, const double time, amrex::MultiFab& ld);

    /** @brief Fills MultiFab with data from valid regions and fills ghost
     *         cells. Works for single level and 2-level cases (interpolating
     *         from coarse).
     * @param   lev     (Fine) Level to be filled.
     * @param   time    Time of level.
     * @param   mf      MultiFab to be filled.
     * @param   scomp   Starting component of source.
     * @param   scomp   Starting component of destination.
     * @param   ncomp   Total number of components.
     */
    void FillPatch(const int lev, const double time, amrex::MultiFab& mf,
                   const int scomp = 0, const int dcomp = 0, int ncomp = -1);

    /** @brief Fills MultiFab with data from valid regions and fills ghost cells
     *         during intermediate time steps. Works for single level and
     *         2-level cases (interpolating from coarse).
     * @param   lev     (Fine) Level to be filled.
     * @param   time    Time of level.
     * @param   mf      MultiFab to be filled.
     * @param   scomp   Starting component of source.
     * @param   scomp   Starting component of destination.
     * @param   ncomp   Total number of components.
     */
    void FillIntermediatePatch(const int lev, const double time,
                               amrex::MultiFab& mf, const int scomp = 0,
                               const int dcomp = 0, const int ncomp = -1);

    /** @brief Average down fine level (lev+1) onto coarse level (lev).
     * @param   lev Coarse Level onto which to be averaged down.
     */
    void AverageDownTo(const int lev);

    /** @brief Compute truncation errors for level lev and saves them in
     *         sim->grid_old[lev]. Also averages down lev onto lev-1 at the
     *         same time.
     * @param   lev Level for which truncation errors are to be computed.
     */
    void ComputeTruncationErrors(const int lev);

  private:
    /** @brief Fetches level data at a given level and time. Needs to be
     *         amrex::Vector not std::vector.
     * @param   lev     Level at which data is to be fetched.
     * @param   time    Time at which data is to be fetched. If time does not
     *                  align with t_old or t_new both states will be returned.
     * @return  Vector with pointer to fetched data.
     */
    amrex::Vector<amrex::MultiFab*> GetLevelData(const int lev, 
                                                 const double time);

    /** @brief Integer array containing the type of boundary condition at each
     *         boundary edge. Needs to be amrex::Vector not std::vector.
     */
    amrex::Vector<amrex::BCRec> bcs;

    /** @brief Pointer to AMReX interpolator to be used between levels.
     */
    amrex::Interpolater* mapper = nullptr;

    /** @brief Pointer to owner on whose data this class operates.
     */
    Sledgehamr* sim;

    /** @brief enum for the various interpolation types.
     */
    enum InterpType {
            PCInterp = 0,
            NodeBilinear,
            CellConservativeLinear,
            CellBilinear,
            CellQuadratic,
            CellConservativeProtected,
            CellConservativeQuartic
    };
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_
