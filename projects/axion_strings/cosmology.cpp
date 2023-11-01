#include "cosmology.h"
#include <sledgehamr_utils.h>

namespace axion_strings {

void Cosmology::Init(sledgehamr::Sledgehamr* owner) {
    sim = owner;
    ParseVariables();
    PrintRefinementTimes();
    SetProjections();
    SetSpectra();
    SetXiMeasurement();
}

void Cosmology::ParseVariables() {
    amrex::ParmParse pp_prj("project");
    pp_prj.get("string_width_threshold", string_width_threshold);

    amrex::ParmParse pp_out("output");
    pp_out.get("spectra_log_min", spectra_log_min);
    pp_out.get("interval_xi_log", interval_xi_log);
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

void Cosmology::SetXiMeasurement() {
    sim->io_module->output.emplace_back(sim->io_module->output_folder, "xi",
            OUTPUT_FCT(Cosmology::WriteXi), interval_xi_log);
    sim->io_module->output.back().SetTimeFunction(TIME_FCT(Cosmology::Log));
}

bool Cosmology::WriteXi(double time, std::string prefix) {
    int lev    = sim->GetFinestLevel();
    double xi  = Xi(lev, time);
    double log = Log(time);
    constexpr int nsize = 4;
    double data[nsize] = {static_cast<double>(lev), time, log, xi};

    amrex::Print() << "Write xi: " << prefix << ", xi=" << xi << std::endl;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/xi.h5";
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
        sledgehamr::IOModule::WriteToHDF5(file_id, "data", data, nsize);
        H5Fclose(file_id); 
    }

    return true; 
}

double Cosmology::Xi(const int lev, const double eta) {
    int string_tags = GetStringTags(lev);

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

int Cosmology::GetStringTags(const int lev) {
    const sledgehamr::LevelData& state = sim->GetLevelData(lev);
    double dx = sim->GetDx(lev);
    double dt = sim->GetDt(lev);
    int ntags = 0;

    // Lets just do this on CPU even if GPUs available. This section is not
    // performace critical and it is a lot simpler this way.
#pragma omp parallel reduction(+: ntags)
    for (amrex::MFIter mfi(state, true); mfi.isValid(); ++mfi) {
        const amrex::Array4<double const>& state_fab = state.array(mfi);
        const amrex::Box& tilebox  = mfi.tilebox();
        const amrex::Dim3 lo = amrex::lbound(tilebox);
        const amrex::Dim3 hi = amrex::ubound(tilebox);

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    ntags += TagCellForRefinement<true>(state_fab, i, j, k, lev,
                                                        state.t, dt, dx, NULL);
                }
            }
        }
    } 

    amrex::ParallelDescriptor::ReduceIntSum(ntags);
    return ntags; 
}

}; // namespace axion_strings
