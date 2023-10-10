#include "axion_strings.h"

namespace axion_strings{

void axion_strings::Init() {
    ParseVariables();
    PrintRefinementTimes();
    SetProjections();
    SetSpectra();
}

bool axion_strings::CreateLevelIf(const int lev, const double time) {
    if (lev <= 0) return true;
    return StringWidth(lev-1, time) <= string_width_threshold;
}

void axion_strings::ParseVariables() {
    amrex::ParmParse pp_prj("project");
    pp_prj.get("string_width_threshold", string_width_threshold);

    amrex::ParmParse pp_out("output");
    pp_out.get("spectra_log_min", spectra_log_min);
}

void axion_strings::PrintRefinementTimes() {
    for (int lev = 1; lev <= max_level; ++lev) {
        amrex::Print() << "Level " << lev << " ("
                       << sledgehamr::utils::LevelName(lev)
                       << ") will be introduced at eta = "
                       << RefinementTime(lev-1) << std::endl;
    }
}

void axion_strings::SetProjections() {
    // Add projections.
    sledgehamr::Projection proj1(a_prime2, "a_prime2");
    sledgehamr::Projection proj2(r_prime2, "r_prime2");
    io_module->projections.push_back(proj1);
    io_module->projections.push_back(proj2);
}

void axion_strings::SetSpectra() {
    // Add spectra and change time interval to log(m_r/H).
    sledgehamr::Spectrum spec1(a_prime_screened, "a_prime_screened");
    io_module->spectra.push_back(spec1);
    io_module->output[io_module->idx_spectra].SetTimeFunction(
            TIME_FCT(axion_strings::LogTruncated));
}

}; // namespace axion_strings
