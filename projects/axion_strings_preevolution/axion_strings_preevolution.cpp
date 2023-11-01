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
    double xi = cosmo.Xi(finest_level, eta_0);
    amrex::Print() << "String length: " << xi << ", target: " << xi_0
                   << std::endl;
    return (xi <= xi_0 && time >= min_eta);
}

void axion_strings_preevolution::ParseConstants() {
    amrex::ParmParse pp("project");
    pp.get("starting_log", log_0);
    pp.get("starting_xi", xi_0);
    pp.get("min_eta", min_eta);
    eta_0 = std::sqrt(std::exp(log_0)/std::sqrt(2.));
}

} // namespace axion_strings_preevolution
