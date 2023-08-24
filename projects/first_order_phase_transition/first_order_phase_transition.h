#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_

#include <sledgehamr.h>
#include <sledgehamr_utils.h>

#include "bubbles.h"

namespace first_order_phase_transition {

ADD_SCALARS(Phi)
ADD_CONJUGATE_MOMENTA(dPhi)

/** @brief Function that calculates the RHS of the EOM at a single cell.
 * @param   rhs     Container to be filled with RHS.
 * @param   state   Data from which to calculate RHS (current state).
 * @param   i       i-th cell index.
 * @param   j       j-th cell index.
 * @param   k       k-th cell index.
 * @param   lev     Current level.
 * @param   time    Current time.
 * @param   dt      Time step size.
 * @param   dx      Grid spacing.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Rhs(const amrex::Array4<double>& rhs,
         const amrex::Array4<const double>& state,
         const int i, const int j, const int k, const int lev,
         const double time, const double dt, const double dx,
         const double* params) {
    double quadratic = params[0];
    double cubic     = params[1];
    double quartic   = params[2];
    double Phi       = state(i, j, k, Scalar::Phi);
    double potential = quadratic*Phi + cubic*Phi*Phi + quartic*Phi*Phi*Phi;

    constexpr int order = 1;
    double laplacian_Phi = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Scalar::Phi, dx*dx);

    rhs(i, j, k, Scalar::Phi)  = state(i, j, k, Scalar::dPhi);
    rhs(i, j, k, Scalar::dPhi) = laplacian_Phi + potential;
}

template<> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void GravitationalWavesRhs<true>(const amrex::Array4<double>& rhs,
        const amrex::Array4<const double>& state, const int i, const int j,
        const int k, const int lev, const double time, const double dt,
        const double dx, const double* params) {
    // Compute Laplacians.
    double dx2 = dx * dx;
    constexpr int order = 1;
    double laplacian_u_xx = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_xx, dx2);
    double laplacian_u_yy = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_yy, dx2);
    double laplacian_u_zz = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_zz, dx2);
    double laplacian_u_xy = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_xy, dx2);
    double laplacian_u_xz = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_xz, dx2);
    double laplacian_u_yz = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Gw::u_yz, dx2);

    // Compute gradients.
    double grad_x_Phi = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Phi, dx, 'x');
    double grad_y_Phi = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Phi, dx, 'y');
    double grad_z_Phi = sledgehamr::utils::Gradient<order>(
            state, i, j, k, Scalar::Phi, dx, 'z');

    // Compute EOM.
    rhs(i, j, k, Gw::u_xx) = state(i, j, k, Gw::du_xx);
    rhs(i, j, k, Gw::u_yy) = state(i, j, k, Gw::du_yy);
    rhs(i, j, k, Gw::u_zz) = state(i, j, k, Gw::du_zz);
    rhs(i, j, k, Gw::u_xy) = state(i, j, k, Gw::du_xy);
    rhs(i, j, k, Gw::u_xz) = state(i, j, k, Gw::du_xz);
    rhs(i, j, k, Gw::u_yz) = state(i, j, k, Gw::du_yz);
    rhs(i, j, k, Gw::du_xx) = laplacian_u_xx + grad_x_Phi*grad_x_Phi;
    rhs(i, j, k, Gw::du_yy) = laplacian_u_yy + grad_y_Phi*grad_y_Phi;
    rhs(i, j, k, Gw::du_zz) = laplacian_u_zz + grad_z_Phi*grad_z_Phi;
    rhs(i, j, k, Gw::du_xy) = laplacian_u_xy + grad_x_Phi*grad_y_Phi;
    rhs(i, j, k, Gw::du_xz) = laplacian_u_xz + grad_x_Phi*grad_z_Phi;
    rhs(i, j, k, Gw::du_yz) = laplacian_u_yz + grad_y_Phi*grad_z_Phi;
}

/** @brief Modifies the truncation error criteria for Pi1 and Pi2 from its
 *         default \tau > \tau_{crit} to \tau * \Delta t_{\ell} > \tau_{crit}.
 * @param   truncation_error    \tau
 * @return f(\tau) for criteria f(\tau) > \tau_{crit}.
 */
template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Scalar::dPhi>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_xx>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_yy>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_zz>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_xy>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_xz>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double TruncationModifier<Gw::du_yz>(const amrex::Array4<const double>& state,
        const int i, const int j, const int k, const int lev, const double time,
        const double dt, const double dx, const double truncation_error,
        const double* params) {
    return truncation_error * dt;
}

/** @brief TODO
 */
AMREX_FORCE_INLINE
double dPhi2(amrex::Array4<double const> const& state, const int i,
        const int j, const int k, const int lev, const double time,
        const double dt, const double dx,
        const std::vector<double>& params) {
    double dPhi = state(i, j, k, Scalar::dPhi);
    return dPhi * dPhi;
}

AMREX_FORCE_INLINE
double Distance(const double a, const double b, const double L) {
    const double dist = fabs(a-b);
    return dist < L/2. ? dist : L - dist;
}

AMREX_FORCE_INLINE
void AddBubble(const int i, const int j, const int k, const double dx,
               const double L, amrex::Array4<amrex::Real> const& fab,
               const Bubble& bubble) {
    const double Dx = Distance(i*dx,bubble.x, L);
    const double Dy = Distance(j*dx,bubble.y, L);
    const double Dz = Distance(k*dx,bubble.z, L);

    double D = std::sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
    double pos = bubble.GetPos(D);
    if( pos == -1 ) {
        return;
    }
    int ind = pos;
    double frac = pos - (double)ind;

    for (int n = 0; n < fab.nComp(); ++n) {
        int n0 = n;
        if (n == Scalar::Phi) n0 = 0;
        else if (n == Scalar::dPhi) n0 = 1;
        else continue;

        double val = bubble.GetVal(n0, ind, frac);
        fab(i, j, k, n) += val;
    }
}

FINISH_SLEDGEHAMR_SETUP

/** @brief Class to simulate a first order phase transition.
 */
class first_order_phase_transition : public sledgehamr::Sledgehamr {
  public:
    START_PROJECT(first_order_phase_transition)

    void Init() override;
    void SetParamsRhs(std::vector<double>& params) override;
    void BeforeTimestep(const double time) override;

  private:
    void ParseVariables();
    void ParseBubbles();
    void ComputeParameters();
    void SetProjections();

    void InjectBubbles(const double time);
    void InjectBubbleLevels(std::vector<int> ab);
    void FillBubbleLayout(const int lev, std::vector<int> ab);
    void AddBubbleValues(std::vector<int> ab);

    double lambda_bar;
    double quadratic, cubic, quartic;

    std::vector<Bubble> bubbles;
    int next_bubble = 0;
    int idx_perfmon_add_bubbles;
    std::vector<int> bubbles_to_inject;
};

}; // namespace first_order_phase_transition

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_H_
