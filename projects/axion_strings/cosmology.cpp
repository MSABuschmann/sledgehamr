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

double Cosmology::Xi(const int string_tags, const int lev, const double eta) {
    const double T1    = 5.4954174414835757e+17;
    const double mpl   = 1.22e19;
    const double gStar = 106;

    double string_length_sim_units = static_cast<double>(string_tags) 
                                        * 2./3. * sim->GetDx(lev);
    double physical_string_length  = BoxToPhysical(string_length_sim_units, eta,
                                                   T1, mpl, gStar);
    double physical_box_size       = BoxToPhysical(sim->GetL(), eta, T1, mpl,
                                                   gStar);

    double T    = XiTemp(eta, T1);
    double time = XiTime(T, mpl, gStar);
    double xi   = physical_string_length * std::pow(time, 2) 
                    / std::pow(physical_box_size, 3);
    return xi;
}

}; // namespace axion_strings
