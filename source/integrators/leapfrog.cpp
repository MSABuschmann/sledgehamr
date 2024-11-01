#include "leapfrog.h"

namespace sledgehamr {

/** @brief Advances one level by one time step using leap-frog algorithm in the
 * kick-drift-kick form.
 * @param   mf_old  Current state.
 * @param   mf_new  New state after advancement.
 * @param   lev     Current level.
 * @param   dt      Time step size.
 * @param   dx      Grid spacing.
 */
void IntegratorLeapfrog::Integrate(LevelData &mf_old, LevelData &mf_new,
                                   const int lev, const double dt,
                                   const double dx) {
    const double t0 = mf_old.t;
    const double t1 = t0 + dt;

    // Total number of fields, gravitational field components and
    // user-defined fields (includes conjugate momenta).
    const int N = mf_old.nComp();
    const int Ngrav = sim->with_gravitational_waves ? 12 : 0;
    const int Nf = N - Ngrav;

    // Start and end points of field pairings for user fields and
    // gravitational wave components.
    const int uN = Nf / 2;
    const int uN0 = 0;
    const int uN1 = uN0 + uN;

    const int gN = Ngrav / 2;
    const int gN0 = Nf;
    const int gN1 = gN0 + gN;

    // temp states.
    amrex::MultiFab a(mf_old.boxArray(), mf_old.DistributionMap(), N,
                      sim->nghost);
    amrex::MultiFab vh(mf_old.boxArray(), mf_old.DistributionMap(), N,
                       sim->nghost);

    // integrate.
    sim->FillRhs(a, mf_old, t0, lev, dt, dx);
    amrex::MultiFab::LinComb(vh, 1, mf_old, 0, dt / 2., a, 0, 0, N, 0);
    amrex::MultiFab::LinComb(mf_new, 1, mf_old, uN0, dt, vh, uN1, uN0, uN, 0);
    if (Ngrav > 0) {
        amrex::MultiFab::LinComb(mf_new, 1, mf_old, gN0, dt, vh, gN1, gN0, gN,
                                 0);
    }
    sim->level_synchronizer->FillIntermediatePatch(lev, t1, mf_new);
    sim->FillRhs(a, mf_new, t1, lev, dt, dx);
    amrex::MultiFab::LinComb(mf_new, 1, vh, uN1, dt / 2., a, uN1, uN1, uN, 0);
    if (Ngrav > 0) {
        amrex::MultiFab::LinComb(mf_new, 1, vh, gN1, dt / 2., a, gN1, gN1, gN,
                                 0);
    }
    sim->level_synchronizer->FillIntermediatePatch(lev, t1, mf_new);
}

void IntegratorLeapfrog::DebugPrint(amrex::MultiFab &mf, const char *msg) {
#pragma omp parallel
    for (amrex::MFIter mfi(mf, true); mfi.isValid(); ++mfi) {
        const amrex::Array4<double> &state_arr = mf.array(mfi);
        const amrex::Box &bx = mfi.tilebox();

        const amrex::Dim3 lo = amrex::lbound(bx);
        const amrex::Dim3 hi = amrex::ubound(bx);

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (i == 0 && j == 0 && k == 0) {
                        amrex::AllPrint()
                            << msg << ": " << state_arr(i, j, k, 0) << " "
                            << state_arr(i, j, k, 1) << " | "
                            << state_arr(i, j, k, 2) << " "
                            << state_arr(i, j, k, 3) << std::flush << std::endl;
                    }
                }
            }
        }
    }
}

}; // namespace sledgehamr
