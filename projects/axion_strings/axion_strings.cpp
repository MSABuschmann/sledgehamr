#include "axion_strings.h"

namespace axion_strings{

void axion_strings::Init() {
    ParseVariables();
    PrintRefinementTimes();
    SetProjections();
}

bool axion_strings::CreateLevelIf(const int lev, const double time) {
    return StringWidth(lev, time) <= string_width_threshold;
}

void axion_strings::ParseVariables() {
    amrex::ParmParse pp("project");
    pp.get("string_width_threshold", string_width_threshold);
}

void axion_strings::PrintRefinementTimes() {
    for (int lev = 1; lev <= max_level; ++lev) {
        amrex::Print() << "Level " << lev << " ("
                       << sledgehamr::utils::LevelName(lev)
                       << ") will be introduced at eta = "
                       << RefinementTime(lev) << std::endl;
    }
}

void axion_strings::SetProjections() {
    sledgehamr::Projection proj1(a_prime2, "a_prime2");
    sledgehamr::Projection proj2(r_prime2, "r_prime2");
    io_module->projections.push_back(proj1);
    io_module->projections.push_back(proj2);

    sledgehamr::Spectrum spec1(a_prime2_screened, "a_prime2_screened");
    io_module->spectra.push_back(spec1);
}

double axion_strings::StringWidth(const int lev, const double time) {
    const double eta = time;
    const double m_r = sqrt(2.*lambda) * eta; 
    return 1./(m_r * dx[lev]); 
}

double axion_strings::RefinementTime(const int lev) {
    return dimN[lev] / (sqrt(2.*lambda) * string_width_threshold * L);
}

}; // namespace axion_strings
