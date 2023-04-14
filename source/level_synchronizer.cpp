#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>

#include "level_synchronizer.h"
#include "kernels.h"

namespace sledgehamr{

LevelSynchronizer::LevelSynchronizer(Sledgehamr* owner) {
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

void LevelSynchronizer::FillCoarsePatch(int lev, double time,
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

    // Interpolate lev from lev-1.
    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, 0, mf.nComp(),
                                 sim->geom[lev-1], sim->geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 sim->refRatio(lev-1), mapper, bcs, 0);
#else
    // Boundary conditions.
    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1], bcs,
                                                       bndry_func);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev  ], bcs,
                                                       bndry_func);

    // Interpolate lev from lev-1.
    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, 0, mf.nComp(),
                                 sim->geom[lev-1], sim->geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 sim->refRatio(lev-1), mapper, bcs, 0);
#endif
}

void LevelSynchronizer::FillPatch(int lev, double time, amrex::MultiFab& mf) {
    // Get data and boundary conditions for level lev.
    amrex::Vector<amrex::MultiFab*> fmfs = GetLevelData(lev, time);
    amrex::Vector<double> ftime = LevelData::getTimes(fmfs);

#ifdef AMREX_USE_GPU
    amrex::GpuBndryFuncFab<NullFill> gpu_bndry_func(NullFill{});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > fphysbc(
            sim->geom[lev], bcs, gpu_bndry_func);

    if (lev == 0) {
        // Call FillPatchSingleLevel for the coarse level.
        amrex::FillPatchSingleLevel(mf, time, fmfs, ftime, 0, 0, mf.nComp(),
                                    sim->geom[lev], fphysbc, 0);
    } else {
        // Call FillPatchTwoLevels with data from fine (lev) and coarse (lev-1)
        // level.
        amrex::Vector<amrex::MultiFab*> cmfs = GetLevelData(lev-1, time);
        amrex::Vector<double> ctime = LevelData::getTimes(cmfs);

        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > cphysbc(
                sim->geom[lev-1], bcs, gpu_bndry_func);

        amrex::FillPatchTwoLevels(mf, time, cmfs, ctime, fmfs, ftime, 0, 0,
                                  mf.nComp(), sim->geom[lev-1], sim->geom[lev],
                                  cphysbc, 0, fphysbc, 0, sim->refRatio(lev-1),
                                  mapper, bcs, 0);
    }
#else
    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(sim->geom[lev], bcs,
                                                       bndry_func);

    if (lev == 0) {
        // Call FillPatchSingleLevel for the coarse level.
        amrex::FillPatchSingleLevel(mf, time, fmfs, ftime, 0, 0, mf.nComp(),
                                    sim->geom[lev], fphysbc, 0);
    } else {
        // Call FillPatchTwoLevels with data from fine (lev) and coarse (lev-1)
        // level.
        amrex::Vector<amrex::MultiFab*> cmfs = GetLevelData(lev-1, time);
        amrex::Vector<double> ctime = LevelData::getTimes(cmfs);

        amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1],
                                                           bcs, bndry_func);

        amrex::FillPatchTwoLevels(mf, time, cmfs, ctime, fmfs, ftime, 0, 0,
                                  mf.nComp(), sim->geom[lev-1], sim->geom[lev],
                                  cphysbc, 0, fphysbc, 0, sim->refRatio(lev-1),
                                  mapper, bcs, 0);
    }
#endif
}

void LevelSynchronizer::FillIntermediatePatch(int lev, double time,
                                              amrex::MultiFab& mf) {
    FillPatch(lev, time, mf);
}

void LevelSynchronizer::AverageDownTo(int lev) {
    amrex::average_down(sim->grid_new[lev+1], sim->grid_new[lev],
                        sim->geom[lev+1], sim->geom[lev], 0,
                        sim->grid_new[lev].nComp(), sim->refRatio(lev));

    // Since we averaged down we do not have truncation errors available at
    // lev+1.
    sim->grid_old[lev+1].contains_truncation_errors = false;
}

void LevelSynchronizer::ComputeTruncationErrors(int lev) {
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
        ifi (amrex::Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate()) {
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

amrex::Vector<amrex::MultiFab*> LevelSynchronizer::GetLevelData(int lev,
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
