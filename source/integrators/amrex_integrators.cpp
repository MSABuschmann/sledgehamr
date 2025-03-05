#include <AMReX_TimeIntegrator.H>

#include "amrex_integrators.h"

namespace sledgehamr {

/** @brief Advances one level by one time step using the Runge-Kutta
 *         integration scheme.
 * @param   mf_old  Current state.
 * @param   mf_new  New state after advancement.
 * @param   lev     Current level.
 * @param   dt      Time step size.
 * @param   dx      Grid spacing.
 */
void IntegratorAMReX::Integrate(LevelData &mf_old, LevelData &mf_new,
                                const int lev, const double dt,
                                const double dx) {

    amrex::Abort("AMReX Integrator interface changed. Sledgehamr needs to be "
                 "adapted to reflect those changes first.");

    // TODO Make 4th order time interpolation available.
    amrex::TimeIntegrator<amrex::MultiFab> integrator(mf_old);
    /*
        auto source_fun = [&](amrex::MultiFab &rhs, const amrex::MultiFab
       &state, const double time) {
            sim->level_synchronizer->FillIntermediatePatch(lev, time, state);
            sim->FillRhs(rhs, state, time, lev, dt, dx);
        };

        //   auto post_update_fun = [&](amrex::MultiFab& S_data, const double
       time)
        //   {
        //       sim->level_synchronizer->FillIntermediatePatch(lev, time,
       S_data);
        //   };

        integrator.set_rhs(source_fun);
        */
    //    integrator.set_post_update(post_update_fun);
    integrator.advance(mf_old, mf_new, mf_old.t, dt);
}

}; // namespace sledgehamr
