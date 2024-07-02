#include "AxionStringsPreevolution.h"

namespace AxionStringsPreevolution {

void AxionStringsPreevolution::SetParamsRhs(std::vector<double> &params,
                                            const double time, const int lev) {
    params.push_back(eta_0);
}

void AxionStringsPreevolution::Init() {
    cosmo.Init(this);
    ParseConstants();

    if (random_state > 0)
        SetRandomState();
}

bool AxionStringsPreevolution::StopRunning(const double time) {
    double xi = cosmo.Xi(finest_level, eta_0);
    amrex::Print() << "String length: " << xi << ", target: " << xi_0
                   << std::endl;
    return (xi <= xi_0 && time >= min_eta);
}

void AxionStringsPreevolution::SetRandomState() {
    sledgehamr::LevelData &state = GetLevelData(0);
    amrex::InitRandom(random_state + amrex::ParallelDescriptor::MyProc());
    amrex::Print() << "Set Random Initial Conditions" << std::endl;

    // #pragma omp parallel
    for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
        const amrex::Array4<double> &state_arr = state.array(mfi);
        const amrex::Box &bx = mfi.tilebox();

        const amrex::Dim3 lo = amrex::lbound(bx);
        const amrex::Dim3 hi = amrex::ubound(bx);

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    state_arr(i, j, k, Scalar::Psi1) =
                        amrex::Random() * 2. - 1.;
                    state_arr(i, j, k, Scalar::Psi2) =
                        amrex::Random() * 2. - 1.;
                }
            }
        }
    }
}

void AxionStringsPreevolution::ParseConstants() {
    amrex::ParmParse pp("project");
    pp.get("random", random_state);
    pp.get("starting_log", log_0);
    pp.get("starting_xi", xi_0);
    pp.get("min_eta", min_eta);
    eta_0 = std::sqrt(std::exp(log_0) / std::sqrt(2.));
}

} // namespace AxionStringsPreevolution
