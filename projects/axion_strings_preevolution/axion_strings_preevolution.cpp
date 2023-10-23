#include "axion_strings_preevolution.h"

namespace axion_strings_preevolution {

void axion_strings_preevolution::SetParamsRhs(std::vector<double>& params) {
    params.push_back(eta_0);
}

void axion_strings_preevolution::Init() {
    cosmo.Init(this);
    ParseConstants();
}

void axion_strings_preevolution::ParseConstants() {
    amrex::ParmParse pp("project");
    pp.get("StartingLog", Log_0);
    pp.get("StartingXi", Xi_0);
    eta_0 = sqrt(exp(Log_0)/sqrt(2.));
}

} // namespace axion_strings_preevolution
