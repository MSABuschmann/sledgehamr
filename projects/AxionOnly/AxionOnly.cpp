#include "AxionOnly.h"

namespace AxionOnly {

void AxionOnly::Init() {
    ParseVariables();
    SetProjections();
    SetSpectrum();
}

void AxionOnly::ParseVariables() {
    amrex::ParmParse pp_prj("project");
    pp_prj.get("n", n);
    pp_prj.get("eta_c", eta_c);
    pp_prj.get("eta_star", eta_star);
    pp_prj.get("N_QCD", N_QCD);
}

void AxionOnly::SetProjections() {
    io_module->projections.emplace_back(dtheta_prime2, "dtheta_prime2");
}

void AxionOnly::SetSpectrum() {
    io_module->spectra.emplace_back(theta_spectrum, "theta");
    io_module->spectra.emplace_back(dtheta_spectrum, "dtheta");
}

void AxionOnly::SetParamsRhs(std::vector<double> &params, const double time,
                             const int lev) {
    const double eta = time;
    const double ma_sq = std::pow(std::min(eta, eta_c) / eta_star, n) / N_QCD;
    params.push_back(ma_sq);
}

}; // namespace AxionOnly
