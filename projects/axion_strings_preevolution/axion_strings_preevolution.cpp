#include "axion_strings_preevolution.h"

namespace axion_strings_preevolution {

void axion_strings_preevolution::SetParamsRhs(std::vector<double>& params) {
    params.push_back(eta_0);
}

void axion_strings_preevolution::Init() {
    cosmo.Init(this);
    ParseConstants();
}

bool axion_strings_preevolution::StopRunning(const double time) {
    const int lev = 0;
    int string_tags = GetStringTags(lev);
    double xi = cosmo.Xi(string_tags, lev, eta_0);
    return (xi <= xi_0 && time >= min_eta);    
}

void axion_strings_preevolution::ParseConstants() {
    amrex::ParmParse pp("project");
    pp.get("starting_log", log_0);
    pp.get("starting_xi", xi_0);
    pp.get("min_eta", min_eta);
    eta_0 = std::sqrt(std::exp(log_0)/std::sqrt(2.));
}

int axion_strings_preevolution::GetStringTags(const int lev) {
    const sledgehamr::LevelData& state = grid_new[lev];
    int ntags = 0;

    // Lets just do this on CPU even if GPUs available. This section is not
    // performace critical and it is a lot simpler this way.
#pragma omp parallel reduction(+: ntags)
    for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
        const amrex::Array4<double const>& state_fab = state.array(mfi);
        const amrex::Box& tilebox  = mfi.tilebox();
        const amrex::Dim3 lo = amrex::lbound(tilebox);
        const amrex::Dim3 hi = amrex::ubound(tilebox);

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    ntags += axion_strings::TagCellForRefinement<true>(
                                    state_fab, i, j, k, lev, state.t, dt[lev],
                                    dx[lev], NULL);
                }
            }
        }
    } 

    amrex::ParallelDescriptor::ReduceIntSum(ntags);
    return ntags; 
}

} // namespace axion_strings_preevolution
