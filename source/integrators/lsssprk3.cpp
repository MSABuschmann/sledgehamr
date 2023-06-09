#include <AMReX_TimeIntegrator.H>

#include "lsssprk3.h"

namespace sledgehamr {

void IntegratorLSSSPRK3::Integrate(LevelData& mf_old, LevelData& mf_new,
        const int lev, const double dt, const double dx) {
    const int ncomp = mf_old.nComp();
    const double t0 = mf_old.t;
    const double t1 = t0 + dt;
    amrex::MultiFab k1(mf_old.boxArray(), mf_old.DistributionMap(), ncomp,
                       sim->nghost);

    sim->FillRhs(k1, mf_old, t0, lev, dt, dx);
    amrex::MultiFab::LinComb(mf_new, 1, mf_old, 0, dt, k1, 0, 0, ncomp, 0);
    sim->level_synchronizer->FillIntermediatePatch(lev, t1, mf_new);
    sim->FillAddRhs(k1, mf_new, t1, lev, dt, dx, 1.);
    amrex::MultiFab::LinComb(mf_new, 1, mf_old, 0, dt/4., k1, 0, 0, ncomp, 0);
    sim->level_synchronizer->FillIntermediatePatch(lev, t0 + dt/2., mf_new);
    sim->FillAddRhs(k1, mf_new, t0 + dt/2., lev, dt, dx, 0.25);
    amrex::MultiFab::LinComb(mf_new, 1, mf_old, 0, dt*2./3., k1, 0, 0, ncomp,0);
    sim->level_synchronizer->FillIntermediatePatch(lev, t1, mf_new);
}

}; // namespace sledgehamr
