#include "rkn.h"
#include "sledgehamr_utils.h"

namespace sledgehamr {

/** @brief Set up integrator in advance.
 * @param   owner   Pointer to simulation.
 * @oaram   id      Integrator type. Must be either RknButcherTableau, Rkn4, or
 *                  Rkn5.
 */
IntegratorRkn::IntegratorRkn(Sledgehamr* owner, const IntegratorType id)
  : Integrator{owner},
    integrator_type(id) {
    SetButcherTableau();
}

/** @brief Advances one level by one time step using the Runge-Kutta-Nystroem
 *         integration scheme.
 * @param   mf_old  Current state.
 * @param   mf_new  New state after advancement.
 * @param   lev     Current level.
 * @param   dt      Time step size. 
 * @param   dx      Grid spacing.
 */
void IntegratorRkn::Integrate(LevelData& mf_old, LevelData& mf_new,
        const int lev, const double dt, const double dx) {
        const double t = mf_old.t;
        const int nghost = sim->nghost;

        // Total number of fields, gravitational field components and
        // user-defined fields (includes conjugate momenta).
        const int N = mf_old.nComp();
        const int Ngrav = sim->with_gravitational_waves ? 12 : 0;
        const int Nf = N - Ngrav;

        // Start and end points of field pairings for user fields and
        // gravitational wave components.
        const int uN = Nf/2;
        const int uN0 = 0;
        const int uN1 = uN0 + uN;

        const int gN = Ngrav/2;
        const int gN0 = Nf;
        const int gN1 = gN0 + gN;

        std::vector<std::unique_ptr<amrex::MultiFab> > F_nodes;
        for (int i = 0; i < number_nodes; ++i) {
            F_nodes.emplace_back( std::make_unique<amrex::MultiFab>(
                    mf_old.boxArray(), mf_old.DistributionMap(), N, nghost) );
        }

        for (int i = 0; i < number_nodes; ++i) {
            // stage_time = x_0 + c_i*h
            double stage_time = t + dt * nodes[i];

            // mf_new = (y_0, y'_0)
            amrex::MultiFab::Copy(mf_new, mf_old, 0, 0, N, nghost);

            // mf_new = (y_0 + c_i * h * y'_0, y'_0)
            amrex::MultiFab::Saxpy(mf_new, dt * nodes[i], mf_old,
                                   uN1, uN0, uN, nghost);
            if (Ngrav > 0)
                amrex::MultiFab::Saxpy(mf_new, dt * nodes[i], mf_old,
                                       gN1, gN0 , gN, nghost);

            if (i > 0) {
                // mf_new = (y_0 + c_i * h * y'_0 + sum_j( h^2 * a_ij * k'_j ),
                //           y'_0)
                for (int j = 0; j < i; ++j) {
                    amrex::MultiFab::Saxpy(
                            mf_new, dt*dt * tableau[i][j], *F_nodes[j],
                            uN1, uN0 , uN, nghost);

                    if (Ngrav == 0) continue;

                    amrex::MultiFab::Saxpy(
                            mf_new, dt*dt * tableau[i][j], *F_nodes[j],
                            gN1, gN0 , gN, nghost);
                }

            }

            sim->level_synchronizer->FillIntermediatePatch(
                    lev, stage_time, mf_new);

            // F_nodes[i] = (null, f(stage_time, mf_new)
            //            = (null, k'_i)
            sim->FillRhs(*F_nodes[i], mf_new, stage_time, lev, dt, dx);
        }

        // mf_new = (y_0, y'_0)
        amrex::MultiFab::Copy(mf_new, mf_old, 0, 0, N, nghost);

        // mf_new = (y_0 + h * y'_0, y'_0)
        amrex::MultiFab::Saxpy(mf_new, dt, mf_new, uN1, uN0, uN, nghost);
        if( Ngrav > 0 )
             amrex::MultiFab::Saxpy(mf_new, dt, mf_new, gN1, gN0, gN, nghost);

        // mf_new = (y_0 + h * y'_0 + sum_i(h^2*\bar[b_i]*k'_i),
        //           y'_0 + sum_i(b_i*k'_i))
        for (int i = 0; i < number_nodes; ++i) {
            amrex::MultiFab::Saxpy(mf_new, dt * dt * weights_bar_b[i],
                                   *F_nodes[i],
                                   uN1, uN0, uN, nghost);
            amrex::MultiFab::Saxpy(mf_new, dt * weights_b[i], *F_nodes[i],
                                   uN1, uN1, uN, nghost);

            if (Ngrav == 0) continue;

            amrex::MultiFab::Saxpy(mf_new, dt * dt * weights_bar_b[i],
                                   *F_nodes[i],
                                   gN1, gN0, gN, nghost);
            amrex::MultiFab::Saxpy(mf_new, dt * weights_b[i], *F_nodes[i],
                                   gN1, gN1, gN, nghost);
        }

        // Call the post-update hook for S_new
        sim->level_synchronizer->FillIntermediatePatch(
                        lev, t + dt, mf_new);
}

/** @brief Sets up the correct Butcher Tableau given an integration scheme.
 */
void IntegratorRkn::SetButcherTableau() {
    switch (integrator_type) {
        case RknButcherTableau:
            ReadUserDefinedButcherTableau();
            break;

        case Rkn4:
            nodes = {0.0,
                     0.5,
                     1.0};
            tableau = {{0.0},
                       {1./8., 0.0},
                       {0.0,   0.5, 0.0}};
            weights_bar_b = {1./6., 1./3., 0.0};
            weights_b     = {1./6., 4./6., 1./6.};
            break;

        case Rkn5:
            nodes = {0.0,
                     1./5.,
                     2./3.,
                     1.0};
            tableau = {{0.0},
                       {1./50., 0.0},
                       {-1./27., 7./27., 0.0},
                       {3./10., -2./35., 9./35., 0.0}};
            weights_bar_b = {14./336., 100./336., 54./336., 0.0};
            weights_b = {14./336., 125./336., 162./336., 35./336.};
            break;
    }

    number_nodes = weights_b.size();
}

/** @brief Reads in a user-defined Butcher Tableau.
 */
void IntegratorRkn::ReadUserDefinedButcherTableau() {
    amrex::ParmParse pp("integrator.rkn");

    pp.getarr("weights_bar_b", weights_bar_b);
    pp.getarr("weights_b", weights_b);
    pp.getarr("nodes", nodes);

    amrex::Vector<amrex::Real> btable; // flattened into row major format
    pp.getarr("tableau", btable);

    if (weights_bar_b.size() != nodes.size() ||
        weights_b.size() != nodes.size()) {
        std::string msg = "integrator.weights should be the same length as ";
        msg += "integrator.nodes";
        amrex::Error(msg);
    } else {
        number_nodes = weights_b.size();
        const int nTableau = (number_nodes * (number_nodes + 1)) / 2;
        if (btable.size() != nTableau) {
            std::string msg = "integrator.tableau incorrect length - ";
            msg += "should include the Butcher Tableau diagonal.";
            amrex::Error(msg);
        }
    }

    // Fill tableau from the flattened entries.
    int k = 0;
    for (int i = 0; i < number_nodes; ++i) {
        std::vector<double> stage_row;
        for (int j = 0; j <= i; ++j) {
            stage_row.push_back(btable[k]);
            ++k;
        }

        tableau.push_back(stage_row);
    }

    // Check that this is an explicit method.
    for (const auto& astage : tableau) {
        if (astage.back() != 0.0) {
            std::string msg = "RKN integrator currently only supports ";
            msg += "explicit Butcher tableaus.";
            amrex::Error(msg);
        }
    }
}

}; // namespace sledgehamr
