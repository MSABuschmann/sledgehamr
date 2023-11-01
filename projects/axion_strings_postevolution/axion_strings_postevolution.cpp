#include "axion_strings_postevolution.h"

namespace axion_strings_postevolution {

void axion_strings_postevolution::SetParamsRhs(std::vector<double>& params) {
    params.push_back(eta_0);
}

void axion_strings_postevolution::Init() {
    cosmo.Init(this);
    ParseConstants();
    ReinterpretInitialState();
}

void axion_strings_postevolution::ParseConstants() {
    amrex::ParmParse pp("project");
    pp.get("starting_log", log_0);
    eta_0 = std::sqrt(std::exp(log_0)/std::sqrt(2.));
}

void axion_strings_postevolution::ReinterpretInitialState() {
    // Set starting time if we are starting a new sim.
    if (!restart_sim) {
        grid_new[0].t = eta_0;
    }
}

} // namespace axion_strings_postevolution
