#ifndef PROJECTS_AXION_STRINGS_H_
#define PROJECTS_AXION_STRINGS_H_

#include <sledgehamr.h>
#include <sledgehamr_utils.h>

namespace axion_strings{

ADD_SCALARS(Psi1, Psi2, Pi1, Pi2);

/** @brief Function that calculates the RHS of the EOM at a single cell.
 * @param   i           i-th cell index.
 * @param   j           j-th cell index.
 * @param   k           k-th cell index.
 * @param   time        Current time.
 * @param   lev         Current level.
 * @param   dx          Grid spacing.
 * @param   rhs_fab     Container to be filled with RHS.
 * @param   state_fab   Data from which to calculate RHS.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Rhs(const int i, const int j, const int k, const double time,
         const int lev, const double dx, amrex::Array4<double> const& rhs_fab,
         amrex::Array4<double const> const& state_fab) {
    // Fetch field values.
    double Psi1 = state_fab(i, j, k, Scalar::Psi1);
    double Psi2 = state_fab(i, j, k, Scalar::Psi2);
    double Pi1  = state_fab(i, j, k, Scalar::Pi1);
    double Pi2  = state_fab(i, j, k, Scalar::Pi2);

    double eta = time;

    // Compute Laplacians.
    double dx2 = dx * dx;
    constexpr int order = 1;
    double LaplacianPsi1 = sledgehamr::utils::Laplacian<order>(
            state_fab, i, j, k, Scalar::Psi1, dx2);
    double LaplacianPsi2 = sledgehamr::utils::Laplacian<order>(
            state_fab, i, j, k, Scalar::Psi2, dx2);

    // Compute EOM.
    double cross_term = eta*eta*( Psi1*Psi1 + Psi2*Psi2 - 1. ) + 0.56233;

    rhs_fab(i, j, k, Scalar::Psi1) =  Pi1;
    rhs_fab(i, j, k, Scalar::Psi2) =  Pi2;
    rhs_fab(i, j, k, Scalar::Pi1)  = -Pi1*2./eta + LaplacianPsi1
                                    - Psi1 * cross_term;
    rhs_fab(i, j, k, Scalar::Pi2)  = -Pi2*2./eta + LaplacianPsi2
                                    - Psi2 * cross_term;
}

/** @brief Checks for zero-crossings between two points in the complex scalar
 *         field.
 * @param   Psi1_1  \Psi_1 of first point.
 * @param   Psi2_1  \Psi_2 of first point.
 * @param   Psi1_2  \Psi_1 of second point.
 * @param   Psi2_2  \Psi_2 of second point.
 * @return Sign of slope of zero-crossing. 0 if no crossing.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int ZeroXing(double Psi1_1, double Psi2_1, double Psi1_2, double Psi2_2) {
    if (Psi2_1 * Psi2_2 >= 0) return 0;
    if (Psi2_1 * Psi1_2 - Psi1_1 * Psi2_2 > 0) return 1;
    return -1;
}

/** @brief Computes the winding factor along a given axis. Will be non-zero if
 *         plaquette is pierced by a string.
 * @param   i           i-th cell index.
 * @param   j           j-th cell index.
 * @param   k           k-th cell index.
 * @param   state_fab   Data.
 * @return  Winding factor.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int WindingAxis1(int i, int j, int k,
                 amrex::Array4<double const> const& state_fab) {
    return ZeroXing(state_fab(i  ,j  ,k  ,Scalar::Psi1),
                    state_fab(i  ,j  ,k  ,Scalar::Psi2),
                    state_fab(i+1,j  ,k  ,Scalar::Psi1),
                    state_fab(i+1,j  ,k  ,Scalar::Psi2))
         + ZeroXing(state_fab(i+1,j  ,k  ,Scalar::Psi1),
                    state_fab(i+1,j  ,k  ,Scalar::Psi2),
                    state_fab(i+1,j+1,k  ,Scalar::Psi1),
                    state_fab(i+1,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state_fab(i+1,j+1,k  ,Scalar::Psi1),
                    state_fab(i+1,j+1,k  ,Scalar::Psi2),
                    state_fab(i  ,j+1,k  ,Scalar::Psi1),
                    state_fab(i  ,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state_fab(i  ,j+1,k  ,Scalar::Psi1),
                    state_fab(i  ,j+1,k  ,Scalar::Psi2),
                    state_fab(i  ,j  ,k  ,Scalar::Psi1),
                    state_fab(i  ,j  ,k  ,Scalar::Psi2));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int WindingAxis2(int i, int j, int k,
                 amrex::Array4<double const> const& state_fab) {
    return ZeroXing(state_fab(i  ,j  ,k  ,Scalar::Psi1),
                    state_fab(i  ,j  ,k  ,Scalar::Psi2),
                    state_fab(i+1,j  ,k  ,Scalar::Psi1),
                    state_fab(i+1,j  ,k  ,Scalar::Psi2))
         + ZeroXing(state_fab(i+1,j  ,k  ,Scalar::Psi1),
                    state_fab(i+1,j  ,k  ,Scalar::Psi2),
                    state_fab(i+1,j  ,k+1,Scalar::Psi1),
                    state_fab(i+1,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state_fab(i+1,j  ,k+1,Scalar::Psi1),
                    state_fab(i+1,j  ,k+1,Scalar::Psi2),
                    state_fab(i  ,j  ,k+1,Scalar::Psi1),
                    state_fab(i  ,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state_fab(i  ,j  ,k+1,Scalar::Psi1),
                    state_fab(i  ,j  ,k+1,Scalar::Psi2),
                    state_fab(i  ,j  ,k  ,Scalar::Psi1),
                    state_fab(i  ,j  ,k  ,Scalar::Psi2));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int WindingAxis3(int i, int j, int k,
                 amrex::Array4<double const> const& state_fab) {
    return ZeroXing(state_fab(i  ,j  ,k  ,Scalar::Psi1),
                    state_fab(i  ,j  ,k  ,Scalar::Psi2),
                    state_fab(i  ,j+1,k  ,Scalar::Psi1),
                    state_fab(i  ,j+1,k  ,Scalar::Psi2))
         + ZeroXing(state_fab(i  ,j+1,k  ,Scalar::Psi1),
                    state_fab(i  ,j+1,k  ,Scalar::Psi2),
                    state_fab(i  ,j+1,k+1,Scalar::Psi1),
                    state_fab(i  ,j+1,k+1,Scalar::Psi2))
         + ZeroXing(state_fab(i  ,j+1,k+1,Scalar::Psi1),
                    state_fab(i  ,j+1,k+1,Scalar::Psi2),
                    state_fab(i  ,j  ,k+1,Scalar::Psi1),
                    state_fab(i  ,j  ,k+1,Scalar::Psi2))
         + ZeroXing(state_fab(i  ,j  ,k+1,Scalar::Psi1),
                    state_fab(i  ,j  ,k+1,Scalar::Psi2),
                    state_fab(i  ,j  ,k  ,Scalar::Psi1),
                    state_fab(i  ,j  ,k  ,Scalar::Psi2));
}

/** @brief Function that tags individual cells for refinement.
 * @param   i           i-th cell index.
 * @param   j           j-th cell index.
 * @param   k           k-th cell index.
 * @param   time        Current time.
 * @param   lev         Current level.
 * @param   state_fab   Data.
 * @return  Boolean value as to whether cell should be refined or not.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool TagCellForRefinement(const int i, const int j, const int k,
                          const double time, const int lev,
                          amrex::Array4<double const> const& state_fab) {
    // Check all three plaquettes (in positive index direction) for string
    // piercings.
    if (WindingAxis1(i, j, k, state_fab) != 0) return true;
    if (WindingAxis2(i, j, k, state_fab) != 0) return true;
    if (WindingAxis3(i, j, k, state_fab) != 0) return true;
    return false;
}

/** @brief Class to simulate axion strings.
 */
class axion_strings : public sledgehamr::Sledgehamr {
  public:
    START_PROJECT(axion_strings)
};

}; // namespace axion_strings

#endif // PROJECTS_AXION_STRINGS_H_
