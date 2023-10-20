#include "cosmology.h"
#include <sledgehamr_utils.h>

namespace axion_strings {

void Cosmology::Init(sledgehamr::Sledgehamr* owner) {
    sim = owner;
    ParseVariables();
    PrintRefinementTimes();
    SetProjections();
    SetSpectra();
}

bool Cosmology::CreateLevelIf(const int lev, const double time) {
    if (lev <= 0) return true;
    return StringWidth(lev-1, time) <= string_width_threshold;
}

void Cosmology::ParseVariables() {
    amrex::ParmParse pp_prj("project");
    pp_prj.get("string_width_threshold", string_width_threshold);

    amrex::ParmParse pp_out("output");
    pp_out.get("spectra_log_min", spectra_log_min);
}

void Cosmology::PrintRefinementTimes() {
    for (int lev = 1; lev <= sim->GetMaxLevel(); ++lev) {
        amrex::Print() << "Level " << lev << " ("
                       << sledgehamr::utils::LevelName(lev)
                       << ") will be introduced at eta = "
                       << RefinementTime(lev-1) << std::endl;
    }
}

void Cosmology::SetProjections() {
    // Add projections.
    sledgehamr::Projection proj1(a_prime2, "a_prime2");
    sledgehamr::Projection proj2(r_prime2, "r_prime2");
    sim->io_module->projections.push_back(proj1);
    sim->io_module->projections.push_back(proj2);
}

void Cosmology::SetSpectra() {
    // Add spectra and change time interval to log(m_r/H).
    sledgehamr::Spectrum spec1(a_prime_screened, "a_prime_screened");
    sim->io_module->spectra.push_back(spec1);
    sim->io_module->output[sim->io_module->idx_spectra].SetTimeFunction(
            TIME_FCT(Cosmology::LogTruncated));
}

}; // namespace axion_strings
