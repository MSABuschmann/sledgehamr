#ifndef SLEDGEHAMR_INTEGRATOR_H_
#define SLEDGEHAMR_INTEGRATOR_H_

#include <AMReX_TimeIntegrator.H>

#include "sledgehamr.h"
#include "level_data.h"

namespace sledgehamr {

template <typename> class Sledgehamr;

/** @brief Pure virtual base class that handles the time integration for a
 *         single level.
 *         TODO: Currently not pure virtual, needs to be expanded for other
 *               integrators.
 */
template <typename T>
class Integrator {
  public:
    Integrator(Sledgehamr<T>* owner);

    /** @brief Advances a single level.
     * @param   lev Level to be advanced.
     */
    virtual void Advance(int lev);

  private:

    /** @brief Pointer to owner on whose data this class operates.
     */
    Sledgehamr<T>* sim;
};

template <typename T>
Integrator<T>::Integrator(Sledgehamr<T>* owner) {
    sim = owner;
}

template <typename T>
void Integrator<T>::Advance(int lev) {
    std::swap(sim->grid_old[lev], sim->grid_new[lev]);
    LevelData& S_new = sim->grid_new[lev];

    // State with ghost cells as the integrator initial condition
    LevelData Sborder(sim->grids[lev], sim->dmap[lev], S_new.nComp(),
                      sim->nghost, sim->grid_old[lev].t);

    sim->level_synchronizer->FillPatch(lev, sim->grid_old[lev].t, Sborder);

    // TODO Only create once if we want to use 4th order time interpolation.
    amrex::TimeIntegrator<amrex::MultiFab> integrator(Sborder);
    const auto geom_lev = sim->geom[lev];

    auto source_fun = [&](amrex::MultiFab& rhs, const amrex::MultiFab& state,
                          const double time) {
        static_cast<T*>(sim)->FillRhs(rhs, state, time, geom_lev, lev);
    };

    auto post_update_fun = [&](amrex::MultiFab& S_data, const double time) {
        sim->level_synchronizer->FillIntermediatePatch(lev, time, S_data);
    };

    integrator.set_rhs(source_fun);
    integrator.set_post_update(post_update_fun);
    integrator.advance(Sborder, S_new, sim->grid_old[lev].t, sim->dt[lev]);

    sim->grid_new[lev].t = sim->grid_old[lev].t + sim->dt[lev];
    sim->grid_new[lev].istep = sim->grid_old[lev].istep + 1;
}

}; // namespace sledgehamr

#endif // SLEDGEHAMR_INTEGRATOR_H_
