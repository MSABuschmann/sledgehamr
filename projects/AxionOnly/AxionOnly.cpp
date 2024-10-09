#include "AxionOnly.h"

namespace AxionOnly {

void AxionOnly::Init() {
    ParseVariables();
    SetProjections();
}

void AxionOnly::ParseVariables() {
    amrex::ParmParse pp_prj("project");
    pp_prj.get("n", n);
    pp_prj.get("eta_c", eta_c);
}

void AxionOnly::SetProjections() {
    io_module->projections.emplace_back(theta_prime2, "theta_prime2");
}

void AxionOnly::SetParamsRhs(std::vector<double> &params, const double time,
                             const int lev) {
    const double eta = time;
    const double eta_growth = std::pow(std::min(eta, eta_c), n);
    params.push_back(eta_growth);
}

}; // namespace AxionOnly
