#include "integrator.h"
#include "level_data.h"

namespace sledgehamr{

Integrator::Integrator(Sledgehamr* owner) {
    sim = owner;
}

void Integrator::Advance(const int lev) {
    if (lev >= 0) {
        sim->grid_old[lev].contains_truncation_errors = false;
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
            return "User-defined RK Butcher Tableau";
        case AmrexForwardEuler:
            return "Forward Euler";
        case AmrexTrapezoid:
            return "Trapezoid Method";
        case AmrexSsprk3:
            return "SSPRK3 (AMReX implementation)";
        case AmrexRk4:
            return "RK4";
        case Lsssprk3:
            return "SSPRK3 (Low-storage sledgehamr implementation)";
        case RknButcherTableau:
            return "User-defined RKN Butcher Tableau";
        case Rkn4:
            return "4th order RKN";
        case Rkn5:
            return "5th order RKN";
        default:
            return "Unkown!";
    }
}

void Integrator::DebugMessage(amrex::MultiFab& mf, std::string msg) { }

};  // namespace sledgehamr
