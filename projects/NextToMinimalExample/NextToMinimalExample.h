#pragma once

#include <sledgehamr.h>

namespace NextToMinimalExample {

SLEDGEHAMR_ADD_SCALARS(Psi1, Psi2)
SLEDGEHAMR_ADD_CONJUGATE_MOMENTA(Pi1, Pi2)

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Rhs(const amrex::Array4<double>& rhs,
         const amrex::Array4<const double>& state,
         const int i, const int j, const int k, const int lev,
         const double time, const double dt, const double dx,
         const double* params) {
    // Retrieve parameter.
    double lambda = params[0];

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

    // Compute potential.
    double potential = eta*eta*( Psi1*Psi1 + Psi2*Psi2 - 1. ) + 0.56233;
    potential *= lambda;

    // Full EOM.
    rhs(i, j, k, Scalar::Psi1) =  Pi1;
    rhs(i, j, k, Scalar::Psi2) =  Pi2;
    rhs(i, j, k, Scalar::Pi1)  = -Pi1*2./eta + laplacian_Psi1 - Psi1*potential;
    rhs(i, j, k, Scalar::Pi2)  = -Pi2*2./eta + laplacian_Psi2 - Psi2*potential;
}


// Extra kernel function to compute \dot{a}^2 = (Psi1*Pi2 - Psi2-Pi1)^2 / r^4
// with r^2 = Psi1*Psi2 + Psi2*Psi2
AMREX_FORCE_INLINE
double a_dot_sq(
        amrex::Array4<amrex::Real const> const& state,
        const int i, const int j, const int k,
        const int lev, const double time, const double dt, const double dx,
        const std::vector<double>& params) {
    double Psi1    = state(i, j, k, Scalar::Psi1);
    double Psi2    = state(i, j, k, Scalar::Psi2);
    double Pi1     = state(i, j, k, Scalar::Pi1);
    double Pi2     = state(i, j, k, Scalar::Pi2);
    double r_sq    = Psi1*Psi1 + Psi2*Psi2;
    double dot_a = (Psi1*Pi2 - Psi2*Pi1)/r_sq;
    return dot_a*dot_a;
}

SLEDGEHAMR_FINISH_SETUP

class NextToMinimalExample : public sledgehamr::Sledgehamr {
  public:
    SLEDGEHAMR_INITIALIZE_PROJECT(NextToMinimalExample)

    void Init() override {
        // Read in parameter from inputs file
        amrex::ParmParse pp("project");
        pp.get("lambda", lambda);

        // Let us output a projection of \dot{a}^2, the power spectrum of it,
        // and also compute the average of the radial mode sqrt(Psi1^2 + Psi2^2)
        io_module->projections.emplace_back(a_dot_sq, "a_dot_sq");
        io_module->spectra.emplace_back(a_dot_sq, "a_dot_sq");
        io_module->output.emplace_back(
                "avg", OUTPUT_FCT(NextToMinimalExample::WriteAvg));
    }

    void SetParamsRhs(std::vector<double>& params, const double time,
                      const int lev) override {
        params.push_back(lambda);
     }

    bool WriteAvg(double time, std::string prefix);

  private:
    double lambda;
};

}; // namespace NextToMinimalExample
