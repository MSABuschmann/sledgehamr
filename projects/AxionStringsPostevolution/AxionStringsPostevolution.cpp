#include "AxionStringsPostevolution.h"

namespace AxionStringsPostevolution {

void AxionStringsPostevolution::SetParamsRhs(
        std::vector<double>& params, const double time, const int lev) {
    double eta = time;
    double eta_pre = eta_pre_0 + (eta - eta_0);
    double frac = 1./(1. + std::exp(-f_transition*(eta-eta_transition)));

    params.push_back(eta_0);
    params.push_back(eta_pre);
    params.push_back(frac);
}

void AxionStringsPostevolution::Init() {
    cosmo.Init(this);
    ParseConstants();
    GetPreevolutionTime();
    ReinterpretInitialState();
}

void AxionStringsPostevolution::ParseConstants() {
    amrex::ParmParse pp("project");
    pp.get("starting_log", log_0);
    pp.get("eta_transition", eta_transition);
    pp.get("f_transition", f_transition);

    eta_0 = std::sqrt(std::exp(log_0)/std::sqrt(2.));
}

void AxionStringsPostevolution::ReinterpretInitialState() {
    // Set starting time if we are starting a new sim.
    if (!restart_sim) {
        grid_new[0].t = eta_0;
    }
}

void AxionStringsPostevolution::GetPreevolutionTime() {
    // Get end time of the preevolution from the initial checkpoint file.
    std::string folder;
    amrex::ParmParse pp("input");
    pp.get("initial_state", folder);

    const int nparams = 8;
    double header[nparams];
    std::string filename = folder + "/Meta.hdf5";
    sledgehamr::IOModule::ReadFromHDF5(filename, {"Header"}, header);
    eta_pre_0 = header[0];
}

} // namespace AxionStringsPostevolution
