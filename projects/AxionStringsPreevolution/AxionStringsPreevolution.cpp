#include "AxionStringsPreevolution.h"

namespace AxionStringsPreevolution {

void AxionStringsPreevolution::SetParamsRhs(
        std::vector<double>& params, const double time, const int lev) {
    params.push_back(eta_0);
}

void AxionStringsPreevolution::Init() {
    cosmo.Init(this);
    ParseConstants();
}

bool AxionStringsPreevolution::StopRunning(const double time) {
    double xi = cosmo.Xi(finest_level, eta_0);
    amrex::Print() << "String length: " << xi << ", target: " << xi_0
                   << std::endl;
    return (xi <= xi_0 && time >= min_eta);
}

void AxionStringsPreevolution::ParseConstants() {
    amrex::ParmParse pp("project");
    pp.get("starting_log", log_0);
    pp.get("starting_xi", xi_0);
    pp.get("min_eta", min_eta);
    eta_0 = std::sqrt(std::exp(log_0)/std::sqrt(2.));
}

} // namespace AxionStringsPreevolution
