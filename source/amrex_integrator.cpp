#include "amrex_integrator.h"

virtual void AMReXIntegrator::Integrate(LevelData& mf_old, LevelData& mf_new, amrex::Geometry& geom, const int lev, const double dt) {
    // TODO Make 4th order time interpolation available.
    amrex::TimeIntegrator<amrex::MultiFab> integrator(mf_old);

    auto source_fun = [&](amrex::MultiFab& rhs, const amrex::MultiFab& state,
                          const double time) {
        sim->FillRhs(rhs, state, time, geom, lev);
    };

    auto post_update_fun = [&](amrex::MultiFab& S_data, const double time) {
        sim->level_synchronizer->FillIntermediatePatch(lev, time, S_data);
    };

    integrator.set_rhs(source_fun);
    integrator.set_post_update(post_update_fun);
    integrator.advance(mf_old, mf_new, mf_old.t, dt);
}
