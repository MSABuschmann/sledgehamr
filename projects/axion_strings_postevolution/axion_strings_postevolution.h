#ifndef PROJECTS_AXION_STRINGS_POSTEVOLUTION_H_
#define PROJECTS_AXION_STRINGS_POSTEVOLUTION_H_

#include <sledgehamr.h>
#include <sledgehamr_utils.h>
#include "../axion_strings/cosmology.h"
#include "../axion_strings/axion_strings.h"

namespace axion_strings_postevolution {
using namespace axion_strings;

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
    // Fetch field values.
    double Psi1 = state(i, j, k, Scalar::Psi1);
    double Psi2 = state(i, j, k, Scalar::Psi2);
    double Pi1  = state(i, j, k, Scalar::Pi1);
    double Pi2  = state(i, j, k, Scalar::Pi2);

    double eta = time;

    // Compute Laplacians.
    constexpr int order = 2;
    double laplacian_Psi1 = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Scalar::Psi1, dx*dx);
    double laplacian_Psi2 = sledgehamr::utils::Laplacian<order>(
            state, i, j, k, Scalar::Psi2, dx*dx);

    // Compute EOM.
    const double eta_0 = params[0];
    const double eta_pre = params[1];
    const double eta_sq = eta_pre*eta_pre;
    const double lambda = eta_0*eta_0;
    const double drag = eta_0;

    double potential_pre = lambda/eta_sq*( Psi1*Psi1 + Psi2*Psi2 - 1. );
    double rhs_Psi1_pre = -Pi1*3./eta_pre + laplacian_Psi1/eta_sq/drag
                                 - Psi1*potential_pre;
    double rhs_Psi2_pre = -Pi2*3./eta_pre + laplacian_Psi2/eta_sq/drag
                                 - Psi2*potential_pre;

    double potential_post = eta*eta*( Psi1*Psi1 + Psi2*Psi2 - 1. ) + 0.56233;
    double rhs_Psi1_post = -Pi1*2./eta + laplacian_Psi1 - Psi1*potential_post;
    double rhs_Psi2_post = -Pi2*2./eta + laplacian_Psi2 - Psi2*potential_post;

    double frac = params[2];
    rhs(i, j, k, Scalar::Psi1) =  Pi1;
    rhs(i, j, k, Scalar::Psi2) =  Pi2;
    rhs(i, j, k, Scalar::Pi1)  = (1.-frac)*rhs_Psi1_pre + frac*rhs_Psi1_post;
    rhs(i, j, k, Scalar::Pi2)  = (1.-frac)*rhs_Psi2_pre + frac*rhs_Psi2_post;
}

using axion_strings::GravitationalWavesRhs;
using axion_strings::TruncationModifier;
using axion_strings::TagCellForRefinement;

FINISH_SLEDGEHAMR_SETUP

/** @brief Class to simulate axion strings.
 */
class axion_strings_postevolution : public sledgehamr::Sledgehamr {
  public:
    START_PROJECT(axion_strings_postevolution)

    void SetParamsRhs(std::vector<double>& params, const double time,
                      const int lev) override;

    void Init() override;

    bool CreateLevelIf(const int lev, const double time) override {
        return cosmo.CreateLevelIf(lev, time);
    };

  private:
    void ParseConstants();
    void GetPreevolutionTime();
    void ReinterpretInitialState();

    double log_0 = 2;
    double eta_0 = 2.3;
    double eta_pre_0 = -1;
    double eta_transition = 2.8;
    double f_transition = 10;
 
    Cosmology cosmo;
};

}; // namespace axion_strings_postevolution

#endif // PROJECTS_AXION_STRINGS_POSTEVOLUTION_H_
