#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>

#include "level_synchronizer.h"
#include "kernels.h"
#include "sledgehamr_utils.h"

namespace sledgehamr{

LevelSynchronizer::LevelSynchronizer(Sledgehamr* owner) : sim(owner) {
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
    int interpolation_type = InterpType::PCInterp;
    pp.query("interpolation_type", interpolation_type);

    switch (interpolation_type) {
        case InterpType::PCInterp:
            mapper = &amrex::pc_interp;
            break;
        case InterpType::CellConservativeLinear:
            mapper = &amrex::cell_cons_interp;
            break;
        case InterpType::CellQuadratic:
            mapper = &amrex::quadratic_interp;
            break;
        case InterpType::CellConservativeQuartic:
            mapper = &amrex::quartic_interp;
            break;
        default:
            amrex::Error("Unsupported interpolation type");
    }
}

void LevelSynchronizer::FillCoarsePatch(const int lev, const double time,
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
    amrex::InterpFromCoarseLevel(mf, amrex::IntVect(sim->nghost), time, *cmf[0],
                                 0, 0, mf.nComp(),
                                 sim->geom[lev-1], sim->geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 sim->refRatio(lev-1), mapper, bcs, 0);
}

void LevelSynchronizer::FillPatch(const int lev, const double time,
                                  amrex::MultiFab& mf, const int scomp,
                                  const int dcomp, int ncomp) {
    sim->performance_monitor->Start(
            sim->performance_monitor->idx_fill_patch, lev);

    if (ncomp == -1 ) {
        ncomp = mf.nComp();
    }

    // Get data and boundary conditions for level lev.
    amrex::Vector<amrex::MultiFab*> fmfs = GetLevelData(lev, time);
    amrex::Vector<double> ftime = LevelData::getTimes(fmfs);
    amrex::Geometry& geom = lev < 0 ? sim->shadow_level_geom : sim->geom[lev];

#ifdef AMREX_USE_GPU
    amrex::GpuBndryFuncFab<NullFill> gpu_bndry_func(NullFill{});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > fphysbc(
            geom, bcs, gpu_bndry_func);
#else
    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(geom, bcs, bndry_func);
#endif

    if (lev <= 0) {
        // Call FillPatchSingleLevel for the coarse level.
        amrex::FillPatchSingleLevel(mf, time, fmfs, ftime, scomp, dcomp, ncomp,
                                    geom, fphysbc, 0);
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

        amrex::FillPatchTwoLevels(mf, time, cmfs, ctime, fmfs, ftime, scomp,
                                  dcomp, ncomp, sim->geom[lev-1],
                                  geom, cphysbc, 0, fphysbc, 0,
                                  sim->refRatio(lev-1), mapper, bcs, 0);
    }

    sim->performance_monitor->Stop(
            sim->performance_monitor->idx_fill_patch, lev);
}

void LevelSynchronizer::FillIntermediatePatch(const int lev, const double time,
        amrex::MultiFab& mf, const int scomp, const int dcomp, int ncomp) {
    sim->performance_monitor->Start(
            sim->performance_monitor->idx_fill_intermediate_patch, lev);

    // Figure out which components to fill.
    if (ncomp == -1 ) {
        ncomp = mf.nComp();
    }

    // Get data and boundary conditions for level lev.
    amrex::Vector<amrex::MultiFab*> fmfs{&mf};
    amrex::Vector<double> ftime{time};
    amrex::Geometry& geom = lev < 0 ? sim->shadow_level_geom : sim->geom[lev];

#ifdef AMREX_USE_GPU
    amrex::GpuBndryFuncFab<NullFill> gpu_bndry_func(NullFill{});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > fphysbc(
            geom, bcs, gpu_bndry_func);
#else
    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(geom, bcs, bndry_func);
#endif

    if (lev <= 0) {
        // Call FillPatchSingleLevel for the coarse level.
        amrex::FillPatchSingleLevel(mf, time, fmfs, ftime, scomp, dcomp, ncomp,
                                    geom, fphysbc, 0);
    } else {
        // Call FillPatchTwoLevels with data from fine (lev) and coarse (lev-1)
        // level.
        amrex::Vector<amrex::MultiFab*> cmfs = GetLevelData(lev-1, time);
        amrex::Vector<double> ctime = LevelData::getTimes(cmfs);

        // TODO check if copy is really needed.
        amrex::MultiFab mf_tmp(mf.boxArray(), mf.DistributionMap(), mf.nComp(),
                         sim->nghost);

#ifdef AMREX_USE_GPU
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<NullFill> > cphysbc(
                sim->geom[lev-1], bcs, gpu_bndry_func);
#else
        amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(sim->geom[lev-1],
                                                           bcs, bndry_func);
#endif

        amrex::FillPatchTwoLevels(mf_tmp, time, cmfs, ctime, fmfs, ftime, scomp,
                                  dcomp, ncomp, sim->geom[lev-1],
                                  geom, cphysbc, 0, fphysbc, 0,
                                  sim->refRatio(lev-1), mapper, bcs, 0);

        std::swap(mf, mf_tmp);
    }

    sim->performance_monitor->Stop(
            sim->performance_monitor->idx_fill_intermediate_patch, lev);
}

void LevelSynchronizer::AverageDownTo(const int lev) {
    sim->performance_monitor->Start(
            sim->performance_monitor->idx_average_down, lev);

    amrex::average_down(sim->grid_new[lev+1], sim->grid_new[lev],
                        sim->geom[lev+1], sim->geom[lev], 0,
                        sim->grid_new[lev].nComp(), sim->refRatio(lev));

    // Since we averaged down we do not have truncation errors available at
    // lev+1.
    sim->grid_old[lev+1].contains_truncation_errors = false;

    sim->performance_monitor->Stop(
            sim->performance_monitor->idx_average_down, lev);
}

void LevelSynchronizer::ComputeTruncationErrors(int lev) {
    sim->performance_monitor->Start(
            sim->performance_monitor->idx_truncation_error, lev);

    // Sanity check shadow level was created at the right time and is properly
    // sync'd with the coarse level.
    if ( lev == 0 &&
         !utils::ApproxEqual(sim->shadow_level.t, sim->grid_new[lev].t) ) {
        std::string msg = "Shadow level not sync'd with coarse level! "
                        + std::to_string(sim->shadow_level.t) + " (shadow) vs "
                        + std::to_string(sim->grid_new[lev].t) + " (coarse)";
        amrex::Abort(msg);
    }

    amrex::MultiFab& S_crse = lev == 0 ? sim->shadow_level
                                       : sim->grid_new[lev-1];
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
            ParallelFor(S_crse, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
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
            ParallelFor(crse_S_fine, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
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

    if (lev == 0)
        sim->shadow_level.clear();

    sim->performance_monitor->Stop(
            sim->performance_monitor->idx_truncation_error, lev);
}

amrex::Vector<amrex::MultiFab*> LevelSynchronizer::GetLevelData(const int lev,
        const double time) {
    amrex::Vector<amrex::MultiFab*> mfs;
    LevelData * New = lev < 0 ? &sim->shadow_level     : &sim->grid_new[lev];
    LevelData * Old = lev < 0 ? &sim->shadow_level_tmp : &sim->grid_old[lev];

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

void LevelSynchronizer::IncreaseCoarseLevelResolution() {
    amrex::Print() << "Increase Coarse Level resolution!" << std::endl;

    if (sim->finest_level > 0) {
        std::string msg = std::string("Increasing coarse level resolution ")
                        + "is currently only supported if grid has not been "
                        + "refined yet.";
        amrex::Abort(msg);
    }

    // Initalize new level.
    const int lev = 0;
    const int ncomp = sim->grid_new[lev].nComp();
    const int nghost = sim->grid_new[lev].nGrow();
    const double time = sim->grid_new[lev].t;

    amrex::Geometry old_geom = sim->geom[lev];
    amrex::Geometry new_geom = amrex::refine(sim->geom[lev],
                                             amrex::IntVect(2,2,2));
    amrex::BoxArray ba = sim->grid_new[lev].boxArray();
    ba.refine(2);

    sim->geom[lev] = new_geom;
    sim->ChopGrids(lev, ba, amrex::ParallelDescriptor::NProcs());
    sim->geom[lev] = old_geom;

    amrex::DistributionMapping dm =
            amrex::DistributionMapping(ba, amrex::ParallelDescriptor::NProcs());
    LevelData ld(ba, dm, ncomp, nghost, time);

    // Interpolate new from old coarse level.
    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(old_geom,bcs,bndry_func);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(new_geom,bcs,bndry_func);

    amrex::Vector<amrex::MultiFab*> cmf{&sim->grid_new[lev]};
    amrex::Vector<amrex::MultiFab*> fmf{&ld};
    amrex::Vector<double> ctime{time};
    amrex::Vector<double> ftime{time};

    amrex::InterpFromCoarseLevel(ld, time, sim->grid_new[lev], 0, 0, ncomp,
                                 old_geom, new_geom, cphysbc, 0, fphysbc, 0,
                                 amrex::IntVect(2,2,2), mapper, bcs, 0);

    // Replace old coarse level.
    sim->grid_new[lev].clear();
    sim->grid_old[lev].clear();

    std::swap(sim->grid_new[lev], ld);
    sim->grid_old[lev].define(ba, dm, ncomp, nghost, time);

    // Update metadata.
    sim->SetBoxArray(lev, ba);
    sim->SetDistributionMap(lev, dm);
    sim->geom[lev] = new_geom;

    sim->max_level -= 1;
    sim->coarse_level_grid_size *= 2;
    sim->dimN.erase(sim->dimN.begin());
    sim->dx.erase(sim->dx.begin());
    sim->dt.erase(sim->dt.begin());

    sim->ReadSpectrumKs(true);
}

void LevelSynchronizer::ChangeNGhost(int new_nghost) {
    // Create new LevelData.
    const int lev                        = 0;
    LevelData& ld_old                    = sim->grid_new[lev];
    const amrex::BoxArray& ba            = ld_old.boxArray();
    const amrex::DistributionMapping& dm = ld_old.DistributionMap();
    const int ncomp                      = ld_old.nComp();
    const double time                    = ld_old.t;
    const amrex::Geometry& geom          = sim->geom[lev];

    // Allocate and fill.
    LevelData ld_new(ba, dm, ncomp, new_nghost, time);

    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> physbc(
            geom, sim->level_synchronizer->bcs, bndry_func);

    amrex::Vector<amrex::MultiFab*> smf{static_cast<amrex::MultiFab*>(&ld_old)};
    amrex::Vector<double> stime{time};
    amrex::MultiFab& mf = ld_new;

    amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp, geom,
                                physbc, 0);

    // Swap and update.
    std::swap(sim->grid_new[lev], ld_new);
    sim->nghost = new_nghost;
}

void LevelSynchronizer::RegridCoarse() {
    const int lev               = 0;
    LevelData& ld_old           = sim->grid_new[lev];
    const int ncomp             = ld_old.nComp();
    const double time           = ld_old.t;
    const amrex::Geometry& geom = sim->geom[lev];

    sim->grid_old[lev].clear();

    // New layout.
    amrex::BoxArray ba = ld_old.boxArray();
    amrex::Box bx = ba.minimalBox();
    ba = amrex::BoxArray(bx);
    sim->ChopGrids(lev, ba, amrex::ParallelDescriptor::NProcs());
    amrex::DistributionMapping dm(ba, amrex::ParallelDescriptor::NProcs());

    // Allocate and fill.
    LevelData ld_new(ba, dm, ncomp, sim->nghost, time);

    amrex::CpuBndryFuncFab bndry_func(nullptr);
    amrex::PhysBCFunct<amrex::CpuBndryFuncFab> physbc(
            geom, sim->level_synchronizer->bcs, bndry_func);

    amrex::Vector<amrex::MultiFab*> smf{static_cast<amrex::MultiFab*>(&ld_old)};
    amrex::Vector<double> stime{time};
    amrex::MultiFab& mf = ld_new;

    amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp, geom,
                                physbc, 0);

    std::swap(sim->grid_new[lev], ld_new);
    sim->grid_old[lev] = LevelData(ba, dm, ncomp, sim->nghost, time);
    sim->SetBoxArray(lev, ba);
    sim->SetDistributionMap(lev, dm);
}

}; // namespace sledgehamr
