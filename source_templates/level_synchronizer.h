#ifndef SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_
#define SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_

#include <AMReX_BCUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include "sledgehamr.h"
#include "level_data.h"
#include "kernels.h"

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

template <typename> class Sledgehamr;

/** @brief This class handles all operations between two levels such as
 *         averaging down, interpolation to fine, filling of ghost cells, etc.
 *         Class is friend of SledgeHAMR.
 */
template <typename T>
class LevelSynchronizer {
  public:
    LevelSynchronizer (Sledgehamr<T>* owner);

    /** @brief Fills LevelData with information from a coarse level. This is
     *         used e.g. when a new level of refinement is added.
     * @param   lev     New level to be filled with data from lev=1.
     * @param   time    Time of new level.
     * @param   ld      New level data.
     */
    void FillCoarsePatch(int lev, double time, amrex::MultiFab& ld);

    /** @brief Fills MultiFab with data from valid regions and fills ghost
     *         cells. Works for single level and 2-level cases (interpolating
     *         from coarse).
     * @param   lev     (Fine) Level to be filled.
     * @param   time    Time of level.
     * @param   mf      MultiFab to be filled.
     */
    void FillPatch(int lev, double time, amrex::MultiFab& mf);

    /** @brief Fills MultiFab with data from valid regions and fills ghost cells
     *         during intermediate time steps. Works for single level and
     *         2-level cases (interpolating from coarse).
     * @param   lev     (Fine) Level to be filled.
     * @param   time    Time of level.
     * @param   mf      MultiFab to be filled.
     */
    void FillIntermediatePatch(int lev, double time, amrex::MultiFab& mf);

    /** @brief Average down fine level (lev+1) onto coarse level (lev).
     * @param   lev Coarse Level onto which to be averaged down.
     */
    void AverageDownTo(int lev);

    /** @brief Compute truncation errors for level lev and saves them in
     *         sim->grid_old[lev]. Also averages down lev onto lev-1 at the
     *         same time.
     * @param   lev Level for which truncation errors are to be computed.
     */
    void ComputeTruncationErrors(int lev);

  private:
    /** @brief Fetches level data at a given level and time. Needs to be
     *         amrex::Vector not std::vector.
     * @param   lev     Level at which data is to be fetched.
     * @param   time    Time at which data is to be fetched. If time does not
     *                  align with t_old or t_new both states will be returned.
     * @return  Vector with pointer to fetched data.
     */
    amrex::Vector<amrex::MultiFab*> GetLevelData(int lev, double time);

    /** @brief Integer array containing the type of boundary condition at each
     *         boundary edge. Needs to be amrex::Vector not std::vector.
     */
    amrex::Vector<amrex::BCRec> bcs;

    /** @brief Pointer to AMReX interpolator to be used between levels.
     */
    amrex::Interpolater* mapper = nullptr;

    /** @brief Pointer to owner on whose data this class operates.
     */
    Sledgehamr<T>* sim;

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

template <typename T>
LevelSynchronizer<T>::LevelSynchronizer(Sledgehamr<T>* owner) {
    sim = owner;

    const int ncomp = sim->scalar_fields.size();

    // Boundary conditions
    bcs.resize(ncomp);
    for (int n = 0; n < ncomp; ++n) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            bcs[n].setLo(i, amrex::BCType::int_dir);
            bcs[n].setHi(i, amrex::BCType::int_dir);
        }
    }

    // Set interpolation type between levels
    amrex::ParmParse pp("amr");
    int interpolation_type = InterpType::CellConservativeQuartic;
    pp.query("interpolation_type", interpolation_type);

    if (interpolation_type == InterpType::CellConservativeLinear) {
        mapper = &amrex::cell_cons_interp;
    } else if (interpolation_type == InterpType::CellConservativeQuartic) {
        mapper = &amrex::quartic_interp;
    } else if (interpolation_type == InterpType::PCInterp) {
        mapper = &amrex::pc_interp;
    } else {
        amrex::Error("Unsupported interpolation type");
    }
}

template <typename T>
void LevelSynchronizer<T>::FillCoarsePatch(int lev, double time,
                                        amrex::MultiFab& mf) {
    // Get lev-1 data.
    std::vector<amrex::MultiFab*> cmf = GetLevelData(lev-1, time);

#ifdef AMREX_USE_GPU
    // Boundary conditions.
    amrex::GpuBndryFuncFab<NullFill> gpu_bndry_func(NullFill{});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > cphysbc(
            sim->geom[lev-1], bcs, gpu_bndry_func);
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > fphysbc(
            sim->geom[lev  ], bcs, gpu_bndry_func);
#else
    // Boundary conditions.
    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1], bcs,
                                                       bndry_func);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev  ], bcs,
                                                       bndry_func);
#endif

    // Interpolate lev from lev-1.
    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, 0, mf.nComp(),
                                 sim->geom[lev-1], sim->geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 sim->refRatio(lev-1), mapper, bcs, 0);
}

template <typename T>
void LevelSynchronizer<T>::FillPatch(int lev, double time, amrex::MultiFab& mf) {
    // Get data and boundary conditions for level lev.
    amrex::Vector<amrex::MultiFab*> fmfs = GetLevelData(lev, time);
    amrex::Vector<double> ftime = LevelData::getTimes(fmfs);

#ifdef AMREX_USE_GPU
    amrex::GpuBndryFuncFab<NullFill> gpu_bndry_func(NullFill{});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > fphysbc(
            sim->geom[lev], bcs, gpu_bndry_func);
#else
    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev], bcs,
                                                       bndry_func);
#endif

    if (lev == 0) {
        // Call FillPatchSingleLevel for the coarse level.
        amrex::FillPatchSingleLevel(mf, time, fmfs, ftime, 0, 0, mf.nComp(),
                                    sim->geom[lev], fphysbc, 0);
    } else {
        // Call FillPatchTwoLevels with data from fine (lev) and coarse (lev-1)
        // level.
        amrex::Vector<amrex::MultiFab*> cmfs = GetLevelData(lev-1, time);
        amrex::Vector<double> ctime = LevelData::getTimes(cmfs);

#ifdef AMREX_USE_GPU
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > cphysbc(
                sim->geom[lev-1], bcs, gpu_bndry_func);
#else
        amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1],
                                                           bcs, bndry_func);
#endif

        amrex::FillPatchTwoLevels(mf, time, cmfs, ctime, fmfs, ftime, 0, 0,
                                  mf.nComp(), sim->geom[lev-1], sim->geom[lev],
                                  cphysbc, 0, fphysbc, 0, sim->refRatio(lev-1),
                                  mapper, bcs, 0);
    }
}

template <typename T>
void LevelSynchronizer<T>::FillIntermediatePatch(int lev, double time,
                                              amrex::MultiFab& mf) {
    FillPatch(lev, time, mf);
}

template <typename T>
void LevelSynchronizer<T>::AverageDownTo(int lev) {
    amrex::average_down(sim->grid_new[lev+1], sim->grid_new[lev],
                        sim->geom[lev+1], sim->geom[lev], 0,
                        sim->grid_new[lev].nComp(), sim->refRatio(lev));

    // Since we averaged down we do not have truncation errors available at
    // lev+1.
    sim->grid_old[lev+1].contains_truncation_errors = false;
}

template <typename T>
void LevelSynchronizer<T>::ComputeTruncationErrors(int lev) {
    amrex::MultiFab& S_crse = sim->grid_new[lev-1];
    amrex::MultiFab& S_fine = sim->grid_new[lev];
    amrex::MultiFab& S_te   = sim->grid_old[lev];
    const int ncomp = sim->scalar_fields.size();

    // Coarsen() the fine stuff on processors owning the fine data.
    amrex::BoxArray crse_S_fine_BA = S_fine.boxArray();
    crse_S_fine_BA.coarsen(2);

    if (crse_S_fine_BA == S_crse.boxArray() &&
        S_fine.DistributionMap() == S_crse.DistributionMap()) {
#ifdef AMREX_USE_GPU
        if (amrex::Gpu::inLaunchRegion() && S_crse.isFusingCandidate()) {
            auto const& crsema = S_crse.arrays();
            auto const& finema = S_fine.const_arrays();
            auto const& tema = S_te.arrays();
            ParallelFor(S_crse, amrex::IntVect(0), S_crse.nComp(),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n)
                noexcept {
                sledgehamr::kernels::AverageDownWithTruncationError(i, j, k,
                        ncomp, crsema[box_no], finema[box_no], tema[box_no]);
            });

            if (!amrex::Gpu::inNoSyncRegion()) {
                amrex::Gpu::streamSynchronize();
            }
        } else
#endif
        {
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
            for (amrex::MFIter mfi(S_crse, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<double> const& crsearr = S_crse.array(mfi);
                amrex::Array4<double const> const& finearr =
                        S_fine.const_array(mfi);
                amrex::Array4<double> const& tearr = S_te.array(mfi);

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    sledgehamr::kernels::AverageDownWithTruncationError(i, j, k,
                            ncomp, crsearr, finearr, tearr);
                });
            }
        }
    } else {
        amrex::MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(),
                                    S_fine.nComp(), 0, amrex::MFInfo(),
                                    amrex::FArrayBoxFactory());
        crse_S_fine.ParallelCopy(S_crse, 0, 0, S_crse.nComp());

#ifdef AMREX_USE_GPU
        if (amrex::Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate()) {
            auto const& crsema = crse_S_fine.arrays();
            auto const& finema = S_fine.const_arrays();
            auto const& tema = S_te.arrays();
            ParallelFor(crse_S_fine, amrex::IntVect(0), ncomp,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n)
                noexcept {
                    sledgehamr::kernels::AverageDownWithTruncationError(i, j, k,
                            ncomp, crsema[box_no], finema[box_no],
                            tema[box_no]);
            });

            if (!amrex::Gpu::inNoSyncRegion()) {
                amrex::Gpu::streamSynchronize();
            }
        } else
#endif
        {
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
            for (amrex::MFIter mfi(crse_S_fine, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi ){
                const amrex::Box& bx = mfi.tilebox();
                amrex::Array4<double> const& crsearr = crse_S_fine.array(mfi);
                amrex::Array4<double const> const& finearr =
                        S_fine.const_array(mfi);
                amrex::Array4<double> const& tearr = S_te.array(mfi);

                // We copy from component scomp of the fine fab into
                // component 0 of the crse fab because the crse fab is a
                // temporary which was made starting at comp 0, it is not part
                // of the actual crse multifab which came in.
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    sledgehamr::kernels::AverageDownWithTruncationError(i, j, k,
                            ncomp, crsearr, finearr, tearr);
                });
            }
        }

        S_crse.ParallelCopy(crse_S_fine, 0, 0, crse_S_fine.nComp());
    }

    // We now have saved truncation errors.
    sim->grid_old[lev].contains_truncation_errors = true;
}

template <typename T>
amrex::Vector<amrex::MultiFab*> LevelSynchronizer<T>::GetLevelData(int lev,
                                                                double time) {
    amrex::Vector<amrex::MultiFab*> mfs;
    LevelData * New = &sim->grid_new[lev];
    LevelData * Old = &sim->grid_old[lev];

    double teps = fabs(New->t - Old->t)*1.e-3;

    // Add either new, old or both states.
    if (time > New->t - teps && time < New->t + teps) {
        mfs.push_back(New);
    } else if (time > Old->t - teps && time < Old->t + teps) {
        mfs.push_back(Old);
    } else {
        mfs.push_back(Old);
        mfs.push_back(New);
    }

    return mfs;
}

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LEVEL_SYNCHRONIZER_H_
