#include "integrator.h"
#include "level_data.h"

namespace sledgehamr{

Integrator::Integrator(Sledgehamr* owner) {
    sim = owner;
}

void Integrator::Advance(const int lev) {
    if (lev >= 0) {
        std::swap(sim->grid_old[lev], sim->grid_new[lev]);
    }

    const double dt   = lev < 0 ? sim->dt[0]*2.         : sim->dt[lev];
    const double dx   = lev < 0 ? sim->dx[0]*2.         : sim->dx[lev];
    LevelData& mf_old = lev < 0 ? sim->shadow_level_tmp : sim->grid_old[lev];
    LevelData& mf_new = lev < 0 ? sim->shadow_level     : sim->grid_new[lev];
    sim->level_synchronizer->FillPatch(lev, mf_old.t, mf_old);

    Integrate(mf_old, mf_new, lev, dt, dx);

    mf_new.t = mf_old.t + dt;
    mf_new.istep = mf_old.istep + 1;

    if (lev < 0) {
        sim->shadow_level_tmp.clear();
    }
}

std::string Integrator::Name(IntegratorType type) {
    switch (type) {
        case AmrexRkButcherTableau:
            return "User-defined RK Butcher Tableau.";
        case AmrexForwardEuler:
            return "Forward Euler.";
        case AmrexTrapezoid:
            return "Trapezoid Method.";
        case AmrexSsprk3:
            return "SSPRK3 (AMReX implementation).";
        case AmrexRk4:
            return "RK4.";
        case Lsssprk3:
            return "SSPRK3 (Low-storage sledgehamr implementation).";
        default:
            return "Unkown!";
    }
}

void Integrator::DebugMessage(amrex::MultiFab& mf, std::string msg) {
#ifndef AMREX_USE_GPU
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid();
            ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto& state = mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            noexcept {
            if (i==0 && j==0 && k==0 ) {
                double Psi1 = state(i, j, k, 0);
                double Psi2 = state(i, j, k, 1);
                double Pi1  = state(i, j, k, 2);
                double Pi2  = state(i, j, k, 3);

                amrex::AllPrint() << msg << ": Psi1" << " " << Psi1 << " " << Psi2 << " "
                                  << Pi1 << " " << Pi2 << " | "
                                  << state(i-1, j-1, k-1, 0) << " "
                                  << std::endl;
            }
        });
    }
#endif
}

};  // namespace sledgehamr
