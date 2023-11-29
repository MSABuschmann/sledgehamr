#ifndef SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_
#define SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_

#include <AMReX_BCUtil.H>
#include <AMReX_Interpolater.H>

#include "sledgehamr.h"
#include "level_data.h"

namespace sledgehamr {

/** @brief Struct with overloaded operator to handle boundary conditions. Empty
 *         because we do not have boundary conditions beyond periodic currently.
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
 *         Class is friend of sledgehamr.
 */
class LevelSynchronizer {
  public:
    LevelSynchronizer(Sledgehamr* owner);

    void FillCoarsePatch(const int lev, const double time, amrex::MultiFab& ld);
    void FillPatch(const int lev, const double time, amrex::MultiFab& mf,
                   const int scomp = 0, const int dcomp = 0, int ncomp = -1);
    void FillIntermediatePatch(const int lev, const double time,
                               amrex::MultiFab& mf, const int scomp = 0,
                               const int dcomp = 0, const int ncomp = -1);

    void AverageDownTo(const int lev);
    void ComputeTruncationErrors(const int lev);

    void IncreaseCoarseLevelResolution();
    void ChangeNGhost(int new_nghost);
    void RegridCoarse();

    /** @brief Integer array containing the type of boundary condition at each
     *         boundary edge. Needs to be amrex::Vector not std::vector.
     */
    amrex::Vector<amrex::BCRec> bcs;

  private:
    amrex::Vector<amrex::MultiFab*> GetLevelData(const int lev,
                                                 const double time);

    /** @brief Pointer to AMReX interpolator to be used between levels.
     */
    amrex::Interpolater* mapper = nullptr;

    /** @brief Pointer to the simulation.
     */
    Sledgehamr* sim;

    /** @brief enum for the various interpolation types.
     */
    enum InterpType {
            PCInterp = 0,
            CellConservativeLinear = 1,
            CellQuadratic = 2,
            CellConservativeQuartic = 4,
    };
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_
